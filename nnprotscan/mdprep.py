import os
import json
import subprocess
import numpy as np
from collections import OrderedDict

import openmm.unit as unit
import openmm.app as app

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB import Atom, Residue, Chain, Model

from .linkage_model import LinkageModel

ROSETTA_LIGAND_RES = "LIG"
MD_LIGAND_RES = "SEM"
LIGAND_CHAIN = "X"
LIGAND_PROXIMITY = 0.2 * unit.nanometer

FF_DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "ff")
VALID_FFS = ["pyRED_amber03"]

def protonate(pdb_file, out_pdb_file="", pH=7.0):
    """
    Protonates the protein using propka-3 with pdb2pqr
    """
        
    if not out_pdb_file:
        out_pdb_file = pdb_file
    rng = int(10000 * np.random.random())
    prefix = f"tmp_{rng}"
    
    cmd = [
        "pdb2pqr30",
        "--keep-chain", "--ff=AMBER", "--titration-state-method=propka", f"--with-ph={pH}", "--protonate-all", 
        "--whitespace", "--quiet",
        f"--pdb-output={out_pdb_file}", 
        out_pdb_file, prefix + ".pqr"
    ]
    
    try:
        subprocess.check_call(cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    except subprocess.CalledProcessError:
        print("pdb2pqr30 failed.")
    
    # get the total charge in the system
    with open(prefix + ".pqr", "r") as of:
        lines = of.readlines()
    charge = sum([float(line.strip().split()[-2]) for line in lines])
    
    for ext in [".pqr", ".log"]:
        if os.path.isfile(prefix + ext):
            os.remove(prefix + ext)
    
    return charge


class TopologyBuilder:
    def __init__(self, pdb_file, match_site, linkage="amide"):
        """
        Build acylated topology from PDB file containing the binder
        with or without other chains (such as target proteins).
        
        Args:
            pdb_file (str): PDB file containing the linker, binder and any other
            protein chains.
            
            match_site (str): Acylation site provided as single letter
            chain id and residue id (1 based). Eg. A12 means acylate
            the 12th residue of chain A.
            
            linkage (str, optional): The linkage model. Defaults to "amide".
        """
        
        self.pdb_file = pdb_file
        self.linkage_model = LinkageModel(linkage)
        self.match_site = {
            "chain": match_site[0], 
            "resid": int(match_site[1:])
        }
        
        self._modeller = None
        self._set_topology_from_pdb(pdb_file)
            
    @property
    def topology(self):
        if self._modeller is None:
            return None
        else:
            return self._modeller.topology
    
    @property
    def positions(self):
        if self._modeller is None:
            return None
        else:
            return self._modeller.positions
    
    @property
    def atom_and_res_map(self):
        atom_map = dict()
        res_map = dict()
        
        if self.topology is None:
            return atom_map, res_map
        
        for c in self.topology.chains():
            for r in c.residues():
                res_map[(c.id, int(r.id))] = r
                for a in r.atoms():
                    atom_map[(c.id, int(r.id), a.name)] = a
        return atom_map, res_map
    
    def _set_topology_from_pdb(self, input_pdb_file):
        p = app.PDBFile(input_pdb_file)
        self._modeller = app.Modeller(p.topology, p.positions)
    
    def _build_noncanonical_residue(self):
        # store ligand atom info
        ligatoms = OrderedDict()
        for k, a in self.atom_and_res_map[0].items():
            if k[0] == LIGAND_CHAIN and k[1] == 1:
                ligatoms[a] = self.positions[a.index]
        
        # sort chains starting with binder (exclude ligand chain)
        new_chainids = sorted([c.id for c in self.topology.chains()])
        new_chainids.remove(self.match_site["chain"])
        new_chainids.remove(LIGAND_CHAIN)
        new_chainids.insert(0, self.match_site["chain"])
        
        # build a new topology and add chains (except the new ligand chain)
        new_top = app.Topology()
        new_chains = {i: new_top.addChain(i) for i in new_chainids}
        
        # add residues and atoms of the binder chain upto the prospective non-canonical residue
        new_residues = dict()
        old2new_map = dict()
        residues_before = [v for k, v in self.atom_and_res_map[1].items()
                           if k[0] == self.match_site["chain"] and k[1] <= self.match_site["resid"]]
        for r in residues_before:
            new_resname = r.name if int(r.id) != self.match_site["resid"] else MD_LIGAND_RES
            new_residues[r.id] = new_top.addResidue(new_resname, new_chains[r.chain.id], r.id)
            for a in r.atoms():
                new_a = new_top.addAtom(a.name, a.element, new_residues[r.id])
                old2new_map[a] = new_a
        
        # now add the ligand atoms to the non-canonical residue
        nca_res = [r for r in new_top.residues() if r.name == MD_LIGAND_RES][0]
        for a in ligatoms:
            new_a = new_top.addAtom(a.name + "x", a.element, nca_res)
            old2new_map[a] = new_a
        
        # finally add all remaining residues and atoms
        residues_after = [v for k, v in self.atom_and_res_map[1].items()
                          if k[0] == self.match_site["chain"] and k[1] > self.match_site["resid"]]
        
        residues_after.extend([
            v for k, v in self.atom_and_res_map[1].items()
            if k[0] not in (self.match_site["chain"], LIGAND_CHAIN)
        ])
        
        for r in residues_after:
            new_residues[r.id] = new_top.addResidue(r.name, new_chains[r.chain.id], r.id)
            for a in r.atoms():
                new_a = new_top.addAtom(a.name, a.element, new_residues[r.id])
                old2new_map[a] = new_a
        
        # add all bonds
        for (a1, a2) in self.topology.bonds():
            new_top.addBond(old2new_map[a1], old2new_map[a2])
        
        # add positions
        new2old_map = {v: k for k, v in old2new_map.items()}
        new_pos = []
        for a in new_top.atoms():
            vec = self.positions[new2old_map[a].index]._value
            new_pos.append(vec)
        new_pos *= unit.nanometer
        
        self._modeller = app.Modeller(new_top, new_pos)
    
    def _link_noncanonical_sidechain(self):
        # remove extra atoms
        atom_map = self.atom_and_res_map[0]
        delatoms = []
        for a in self.linkage_model.extra_atoms["protein"]:
            key = (self.match_site["chain"], self.match_site["resid"], a)
            if key in atom_map:
                delatoms.append(atom_map[key])
                
        for a in self.linkage_model.extra_atoms["protractor"]:
            key = (self.match_site["chain"], self.match_site["resid"], a + "x")
            if key in atom_map:
                delatoms.append(atom_map[key])
        
        if delatoms:
            self._modeller.delete(delatoms)
        
        # connect the sidechain to the linker
        atom_map = self.atom_and_res_map[0] # update atom map!
        for bondpair in self.linkage_model.bond_atoms:
            key_a = (self.match_site["chain"], self.match_site["resid"], bondpair["protein"])
            a = atom_map[key_a]
            
            key_b = (self.match_site["chain"], self.match_site["resid"], bondpair["protractor"] + "x")
            b = atom_map[key_b]
            
            self._modeller.topology.addBond(a, b)
            self._link_coords(a, b)
        
    def _link_coords(self, a, b):
        # a is protein, b is protractor
        coords_a = self._modeller.positions[a.index]
        coords_b = self._modeller.positions[b.index]
        dist = np.linalg.norm(coords_a - coords_b)
        
        if dist <= LIGAND_PROXIMITY:
            return
        shift = dist - LIGAND_PROXIMITY
        
        atom_map = self.atom_and_res_map[0]
        ligand_indices = [a.index for k, a in atom_map.items()
                          if k[0] == LIGAND_CHAIN]
        for i in ligand_indices:
            self._modeller.positions[i] += shift * (coords_a - coords_b) / dist
    
    
    def protonate(self, pH=7.0):
        rng = int(10000 * np.random.random())
        prefix = f"tmp_{rng}"
        tmp_pdb_file = prefix + ".pdb"
        
        self.write(tmp_pdb_file)
        charge = protonate(tmp_pdb_file, pH=pH)
        self._set_topology_from_pdb(tmp_pdb_file)
        
        if os.path.isfile(tmp_pdb_file):
            os.remove(tmp_pdb_file)
        
        return charge
        
    def acylate(self):
        """
        Connects both the topologies as well as coordinates of the protein
        and ligand.
        """
        
        self._build_noncanonical_residue()
        self._link_noncanonical_sidechain()
    
    def write(self, out_pdb_file):
        """
        Write the current topology and coordinates to disk.
        
        Args:
            out_pdb_file (str): Ouptut PDB file.
        """
        with open(out_pdb_file, "w") as of:
            app.PDBFile.writeFile(self.topology, self.positions, file=of,
                                  keepIds=True)
    


def convert_to_gmx(pdb_file, chain_id, ff, out_pdb_file=""):
    """
    Convert a pdb file built using TopologyBuilder to a Gromacs (GMX) compatible format.
    This assumes that the PDB file already has the protractor added. 
    
    Args:
        pdb_file (str): PDB file with protractor linked covalently.
        
        chain_id (str): Chain ID of the (protracted) binder in the pdb_file.
        
        ff (str): Name of forcefield type for the protractor. This should match ffs kept in ../data/ff.
        For this function to work, the ff path must contain at least:
            - the template for the new residue ('new_residue.pdb')
            
            - a json file (residue_name_mapping.json) that maps 
            atom names of the protracted residue in given pdb file to those in the template.
        
        out_pdb_file (str, optional): Output PDB file. If nothing is provided, input will be over-written.
    """
    
    # helper function to map to the reference residue
    def _get_new_residue(r_old, ref_atoms, atom_map, offset):
        atomdict = {a.name : a for a in r_old.get_atoms()}
        inv_atom_map = {v: k for k,v in atom_map.items()}
        
        r_new = Residue.Residue(
            id=(" ", *r_old.id[1:]), # coerce this to be an ATOM even if originally a HETATM
            resname=MD_LIGAND_RES,
            segid=r_old.segid
        )
        
        for i, a in enumerate(ref_atoms):
            if a.name in inv_atom_map:
                name_old = inv_atom_map[a.name]
            else:
                name_old = a.name
            a_old = atomdict[name_old]
            
            a_new = Atom.Atom(
                # retain
                coord=a_old.coord,
                
                # replace
                name=a.name, serial_number=offset + i,
                bfactor=a.bfactor, occupancy=a.occupancy, altloc=a.altloc,
                fullname=a.fullname, element=a.fullname,
                pqr_charge=a.pqr_charge, radius=a.radius
            )
            
            r_new.add(a_new)
        
        return r_new
        
    
    if not out_pdb_file:
        out_pdb_file = pdb_file
    
    if ff not in VALID_FFS:
        valid_ff_str = ", ".join(VALID_FFS)
        msg = f"Protractor forcefield {ff} not found. Currently implemented forcefields are: {valid_ff_str}"
        raise TypeError(msg)
    
    if not os.path.isdir(os.path.join(FF_DATA_PATH, ff)):
        msg = f"No data found for protractor forcefield {ff}"
        raise ValueError(msg)
    
    # get all relevant ff data    
    mapping_file = os.path.join(FF_DATA_PATH, ff, "residue_name_mapping.json")
    with open(mapping_file, "r") as of:
        atom_map = json.load(of)
    
    ref_pdb_file = os.path.join(FF_DATA_PATH, ff,  "new_residue.pdb")
    ref_model = PDBParser(QUIET=True).get_structure("x", ref_pdb_file)[0]
    assert len(list(ref_model.get_chains())) == 1
    ref_atoms = list(ref_model.get_atoms())
    
    # parse given pdb file
    model = PDBParser(QUIET=True).get_structure("x", pdb_file)[0]
    
    # get offset
    offset = -1
    for a in model.get_atoms():
        r = a.get_parent()
        c = r.get_parent()
        if r.resname == MD_LIGAND_RES and c.id == chain_id:
            offset = a.serial_number
            break
    if offset < 0:
        raise ValueError("Cannot determine positive residue offset from given PDB file when converting to GMX")
    
    # build new model
    model_new = Model.Model(0)
    for chain in model.get_chains():
        chain_new = Chain.Chain(chain.id)
        for r in chain.get_residues():
            if chain.id == chain_id and r.resname == MD_LIGAND_RES:
                r_new = _get_new_residue(r, ref_atoms, atom_map, offset)
            else:
                r_new = r
            chain_new.add(r_new)
        model_new.add(chain_new)
    
    io = PDBIO()
    io.set_structure(model_new)
    io.save(out_pdb_file)
    