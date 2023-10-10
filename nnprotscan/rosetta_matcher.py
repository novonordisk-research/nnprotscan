import os
import glob
import subprocess

import pandas as pd
import numpy as np
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.SeqUtils import seq1 as seq3to1
from Bio.PDB import PDBParser

from .linkage_model import LinkageModel

LIGAND_RES = "LIG"
LIGAND_CHAIN = "X"
INTERFACE_CUTOFF = 10.0 # A

MATCH_CST_SCORE_CUTOFF = 20.0
INTERFACE_SCORE_CUTOFF = 0.0
NUM_CONFORMERS_CUTOFF = 1

PACKING_FLAGS = [
    "-ex1",
    "-ex2",
    "-ex1aro",
    "-ex2aro",
    "-no_optH false",
    "-flip_HNQ",
    "-use_input_sc",
    "-linmem_ig 10",
    "-no_his_his_pairE"
]

def _check_rosetta_binary(binary, relpath=""):
    msg = "Rosetta binary path not found. Add this path using the environment variable ROSETTA_HOME"
    
    if not "ROSETTA_HOME" in os.environ:
        raise IOError(msg)
    
    if not os.path.isdir(os.environ["ROSETTA_HOME"]):
        raise IOError(msg)
    
    full_binary_path = os.path.join(os.environ["ROSETTA_HOME"], relpath, binary)
    if not os.path.isfile(full_binary_path):
        raise IOError(msg)
    
    return full_binary_path
    
    
def _run_shell_cmd(cmd, from_path=""):
    curr_path = os.getcwd()
    
    if from_path:
        os.chdir(from_path)
    
    if not isinstance(cmd, list):
        cmd =[cmd]
    subprocess.check_call(" ".join(cmd), shell=True)
    
    #TODO: make the above work without the shell=True
    # and then make the following version work
    # p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = p.communicate()
    # print(stderr.decode("utf-8"))
    
    if from_path:
        os.chdir(curr_path)


def standardize_smiles(smiles):
    """
    Standardize an input SMILES string for a protractor. The SMILES string
    is converted to a rdkit mol, which is processed, deprotonated and 
    converted back into a SMILES without hydrogens.
    
    Args:
        smiles (str): SMILES string for the protractor (must include both
        fatty acid + linker)

    Returns:
        (str): Standardized SMILES string.
    """
    smiles = Chem.CanonSmiles(smiles)
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.RemoveHs(mol, sanitize=True)
    mol = Chem.AddHs(mol)
    
    # deprotonate
    for indices in mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O")):
        atom = mol.GetAtomWithIdx(indices[2])
        atom.SetFormalCharge(-1)
    
    mol = Chem.RemoveHs(mol, sanitize=True)
    return Chem.MolToSmiles(mol)


def generate_small_mol_protractor_conformations(smiles, out_sdf_file, n=1):
    """
    Generates conformers of a small molecule protractor (ususally a fatty acid
    + some linker) from a SMILES description.
    
    Args:
        smiles (str): SMILES string for the protractor (must include both
        fatty acid + linker)
        
        out_sdf_file (str): Conformers will be written to this SDF file.
        
        n (int, optional): Number of conformers to generate. Defaults to 1.
    """
    
    smiles = standardize_smiles(smiles)
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # generate conformers and write to sdf file
    # taken from here: https://gist.github.com/tdudgeon/b061dc67f9d879905b50118408c30aac
    ids = AllChem.EmbedMultipleConfs(
        mol, numConfs=n,
        maxAttempts=10000,
        pruneRmsThresh=1, 
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True,
        enforceChirality=True
    )
    
    writer = Chem.SDWriter(out_sdf_file)
    for i in range(n):
        writer.write(mol, confId=i)


def relax_protein(pdb_file, out_prefix=""):
    """
    Relax a protein according to the rosetta energy function.
    
    Args:
        pdb_file (str): Input PDB file.
        
        out_prefix (str, optional): Prefix for output PDB and score files.
        If not provided, the input PDB file is over-written and a score file
        with the same prefix as the PDB file is produced.
    """
    
    binary = _check_rosetta_binary(
        "relax.static.linuxgccrelease",
        relpath="bin"
    )
    
    if not out_prefix:
        out_prefix = pdb_file.split(".pdb")[0]
    
    cmd = [
        binary,
        
        f"-in:file:s {pdb_file}",
        f"-out:path:pdb {os.path.dirname(out_prefix)}",
        f"-out:file:scorefile {out_prefix}.score.sc",
        "-out:file:output_pose_energies_table false",
        "-out:file:output_pose_cache_data false",
        "-score:weights ref2015",
        "-nstruct 1",
         
        "-relax:constrain_relax_to_start_coords",
        "-relax:coord_constrain_sidechains",
        "-relax:fast",
        "-relax:min_type dfpmin_armijo_nonmonotone",
        "-relax:ramp_constraints false",
        
        *PACKING_FLAGS
    ]
    
    _run_shell_cmd(cmd)
    
    # rename relaxed pdb file
    src_pdb_file = os.path.join(
            os.path.dirname(out_prefix),
            os.path.basename(pdb_file).replace(".pdb", "_0001.pdb")
    )
    target_pdb_file = out_prefix + ".pdb"
    os.rename(src=src_pdb_file, dst=target_pdb_file)
    

def prepare_linker_params(sdf_file, out_prefix="", linkage="amide"):
    """
    Prepare small molecule rosetta params for the linker + biopolymer 
    (usually a fatty acid)
    
    Args:
        sdf_file (str): SDF file containing ligand conformers.
        
        out_prefix (str, optional): Output prefix for param and constraint
        files ultimately used in the matcher protocol. Absolute paths 
        can be used to produce output files in intended folders. If nothing
        is provided, files will be created in the same path as the sdf file.
        
        linkage (str, optional): Either "amide" (default) for lysine
        based acylation or "sulfide" for cysteine based acylation.
    """
    
    binary = _check_rosetta_binary(
        "molfile_to_params.py",
        relpath=os.path.join("scripts", "python", "public")
    )
    
    if not out_prefix:
        out_prefix = sdf_file.split(".sdf")[0]
    
    # write params file
    cmd = [
        binary,
        sdf_file,
        f"-n {LIGAND_RES}",
        f"-p {out_prefix}",
        "--recharge=0",
        "--clobber"
    ]
    
    _run_shell_cmd(cmd)
    
    # write constraint file
    lm = LinkageModel(linkage)
    src_file = lm.match_constraint_file
    target_file = out_prefix + ".cst"
    with open(src_file, "r") as of:
        s = of.read()
    with open(target_file, "w") as of:
        of.write(s)
        

def _run(input_files, output_path, rms_cutoff):
    
    binary = _check_rosetta_binary(
        "match.static.linuxgccrelease",
        relpath="bin"
    )
    
    pdb_file, params_file, constraint_file, position_file = input_files
    
    cmd = [
        binary,
        f"-s {pdb_file}",
        
        "-match:output_format CloudPDB",
        "-match:only_enumerate_non_match_redundant_ligand_rotamers true",
        "-match:consolidate_matches",
        "-out:file:output_pose_energies_table false",
        "-out:file:output_pose_cache_data false",
        
        f"-match:lig_name {LIGAND_RES}",
        f"-match:scaffold_active_site_residues {position_file}",
        f"-match:geometric_constraint_file {constraint_file}",
        f"-extra_res_fa {params_file}",
        
        "-match:match_grouper SameSequenceAndDSPositionGrouper",
        "-match:output_matches_per_group 100",
        f"-match:grouper_downstream_rmsd {rms_cutoff}",
        "-mute protocols.idealize",

        *PACKING_FLAGS
    ]
    
    _run_shell_cmd(cmd, from_path=output_path)
    
    
def run(pdb_file, params_file, constraint_file, output_path="matcher_outputs",
        chain="A", fixed_positions=[], rms_cutoff=3.0, n_cores=8):
    
    """
    Run the Rosetta Matcher protocol for a given number of residue positions.

    Args:
        pdb_file (str): (Relaxed) PDB file.
        
        params_file (str): File containing small molecule params for the linker
        + biopolymer (usually a fatty acid).
        
        constraint_file (str): File containing constraints defining the
        geometry of the binding site.
        
        output_path (str, optional): Output path. Will be created if it
        doesn't exist.
        
        chain (str, optional): Chain along which acylation positions were
        scanned. Defaults to "A".
        
        fixed_positions (list, optional): List of residue positions (in pdb numbering) to not scan.
        
        rms_cutoff (float, optional): Structural clustering cutoff for matcher solutions.
        Defaults to 3.0 A.
        
        n_cores (int, optional): Number of (cpu) cores for parallel processing.
        Defaults to 8.
    """
    
    # get a list of all positions to scan
    model = PDBParser(QUIET=True).get_structure("x", pdb_file)[0][chain]
    positions = sorted(
        [r.id[1] for r in model.get_residues()
        if (r.id[0] == " ") and (r.id[1] not in fixed_positions)]
    )
    
    # use absolute paths everywhere
    output_path = os.path.abspath(output_path)
    pdb_file = os.path.abspath(pdb_file)
    params_file = os.path.abspath(params_file)
    constraint_file = os.path.abspath(constraint_file)
    
    # create output path
    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)
        
    # divide the positions into chunks
    if len(positions) < n_cores:
        n_cores = len(positions)
        chunksize = len(positions)
        poslist = [positions]
    else:
        chunksize = int(len(positions) / n_cores)
        poslist = []
        for i in range(n_cores):
            start = i * chunksize
            stop = (i+1) * chunksize
            if i == n_cores-1:
                stop = len(positions)
            this_positions = [positions[j] for j in range(start, stop)]
            poslist.append(this_positions)
    
    print(f"Using {n_cores} non-redundant cores with {chunksize} positions per core")
    
    # create position files
    posfiles = []
    for i, positions in enumerate(poslist):
        s = " ".join([str(p) for p in positions])
        posfile = os.path.join(output_path, f"positions_chunk_{i}.txt")
        with open(posfile, "w") as of:
            of.write(s)
        posfiles.append(posfile)
    
    arglist = [
        ((pdb_file, params_file, constraint_file, posfile),         
          output_path,
          rms_cutoff) 
        for posfile in posfiles
    ]
    
    with Pool(n_cores) as pool:    
        pool.starmap(_run, arglist)



def _score(input_files, output_path, chunk_id):
    binary = _check_rosetta_binary(
        "rosetta_scripts.static.linuxgccrelease",
        relpath="bin"
    )
    
    print("Running chunk", chunk_id)
    
    # files and paths
    pdblist_file, params_file, xml_file = input_files
    
    # scorefile
    scorefile = os.path.join(output_path, f"score_{chunk_id}.sc") 
    
    # pdb output path
    output_pdb_path = os.path.join(output_path, "pdbs")
    if not os.path.isdir(output_pdb_path):
        os.makedirs(output_pdb_path, exist_ok=True)
        
    cmd = [
        binary,
        f"-l {pdblist_file}",
        "-parser_read_cloud_pdb true",
        "-obey_ENDMDL true",
        "-run:preserve_header true",
        "-overwrite true",
        
        f"-out:path:pdb {output_pdb_path}", 
        f"-out:file:scorefile {scorefile}",
        "-out:file:output_pose_energies_table false",
        "-out:file:output_pose_cache_data false",
        
        f"-extra_res_fa {params_file}",
        f"-parser:protocol {xml_file}",
        
        *PACKING_FLAGS
    ]
    
    # TODO: remove bare exception when this step becomes bug-free
    try:
        _run_shell_cmd(cmd)
    except Exception as e:
        return
    
    
def score(matched_pdb_path, params_file, constraint_file, 
          chain="A", output_path="scores", n_cores=8):
    """
    Score the matched pdbs for each acylation position according to the 
    Rosetta energy function. For each acylation position that has been
    succesfully matched, the score can be computed either as an avg. over
    upto the first 100 models in the cluster or the first model.
    
    Args:
        match_path (str): Path containing all poses produced by the Matcher
        protocol.
        
        params_file (str): File containing small molecule arams for the linker
        + biopolymer (usually a fatty acid)
        
        constraint_file (str): File containing constraints defining the
        geometry of the binding site.
        
        chain (str, optional): Chain along which acylation positions were
        scanned. Defaults to "A".
        
        output_path (str, optional): Output path containing scores for 
        each cluster. Defaults to "scores".
        
        n_cores (int, optional): Number of (cpu) cores for parallel processing.
        Defaults to 8.
    """
    
    xml_str ="""<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction
            name="hard_rep"
            weights="ref2015_cst" />
    </SCOREFXNS>
    
    <RESIDUE_SELECTORS>
        <Chain name="binder" chains="{binder_chain}" />
        <Chain name="ligand" chains="{ligand_chain}" />
        <Logic name="receptor" selector="not (binder or ligand)" />
        
        <InterfaceByVector
            name="iface_protein_selector"
            cb_dist_cut="{iface_cutoff}"
            vector_dist_cut="{iface_cutoff}"
            grp1_selector="binder"
            grp2_selector="receptor "/>
    </RESIDUE_SELECTORS>
    
    <TASKOPERATIONS>
        <DetectProteinLigandInterface
            name="iface_ligand_taskop"
            cut1="{iface_cutoff}"
            design="0" />
        
        <RestrictToRepacking
            name="repack_only" />
    </TASKOPERATIONS>
    
    <MOVERS>
        <PackRotamersMover
            name="repack"
            scorefxn="hard_rep"
            nloop="10"
            task_operations="repack_only,iface_ligand_taskop" />
        
        <MinMover
            name="min"
            scorefxn="hard_rep"
            bb="0"
            chi="1"
            chi_task_operations="iface_ligand_taskop"
            type="dfpmin_armijo_nonmonotone"
            tolerance="0.01"
            max_iter="1000" />
        
        <AddOrRemoveMatchCsts
            name="cstadd"
            cst_instruction="add_new"
            cstfile="{constraint_file}"
            keep_covalent="1" />
    </MOVERS>
    
    <FILTERS>
        <EnzScore
            name="allcst"
            score_type="cstE"
            scorefxn="hard_rep"
            whole_pose="1"
            energy_cutoff="10000"
            confidence="0" />
            
        <ScorePoseSegmentFromResidueSelectorFilter
            name="iface"
            residue_selector="iface_protein_selector"
            scorefxn="hard_rep"
            in_context="1"
            confidence="0" />
    </FILTERS>
    
    <PROTOCOLS>
        <Add mover_name="cstadd" />
        
        <Add mover_name="min" />
        <Add mover_name="repack" />
        <Add mover_name="min" />
        
        <Add filter_name="allcst" />
        <Add filter_name="iface" />
  </PROTOCOLS>
</ROSETTASCRIPTS>
    """
    
    # use absolute paths everywhere
    matched_pdb_path = os.path.abspath(matched_pdb_path)
    output_path = os.path.abspath(output_path)
    params_file = os.path.abspath(params_file)
    constraint_file = os.path.abspath(constraint_file)
    
    # create output path
    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)
    
    # write the rosetta script xml file
    xml_file = os.path.join(output_path, "min_cst.xml")
    if not os.path.isfile(xml_file):
        with open(xml_file, "w") as of:
            of.write(xml_str.format(
                constraint_file=constraint_file,
                binder_chain=chain,
                ligand_chain=LIGAND_CHAIN,
                iface_cutoff=INTERFACE_CUTOFF
            ))
    
    # divide the matched pdb files into chunks of pdblists for each chunk.
    # This may be slow.
    print("Preparing lists of matched pdb files")
    cloudpdb_files = glob.glob(os.path.join(matched_pdb_path, "*_1.pdb"))
    n_files = len(cloudpdb_files)
    pdblists = []
    
    if n_files < n_cores:
        n_cores = n_files
        chunksize = n_files
        pdblists = [cloudpdb_files]
    else:
        chunksize = int(n_files/ n_cores)
        for i in range(n_cores):
            start = i * chunksize
            stop = (i+1) * chunksize
            if i == n_cores-1:
                stop = n_files
            pdblist = cloudpdb_files[start : stop]
            pdblists.append(pdblist)
    
    pdblist_files = []
    for i, pdblist in enumerate(pdblists):
        f = os.path.join(output_path, f"pdblist_{i}.txt")
        with open(f, "w") as of:
            of.write("\n".join(pdblist))
        pdblist_files.append(f)
        
    print(f"Using {n_cores} non-redundant cores, scoring {chunksize} matched pdb files per core")
    
    arglist = [
        ((pdblist_file, params_file, xml_file),
         output_path, i) 
        for i, pdblist_file in enumerate(pdblist_files)
    ]
    
    with Pool(n_cores) as pool:    
        pool.starmap(_score, arglist)
    
    # get score statistics
    score_files = glob.glob(os.path.join(output_path, "score_*.sc"))
    df = []
    for sf in score_files:
        with open(sf, "r") as of:
            lines = [line.strip() for line in of.readlines()[2:]]
        for line in lines:
            l = line.split()
            cst_score = float(l[2])
            iface_score = float(l[20])
            cluster_idx = int(l[-1].split("_")[1])
            pos_idx = int(l[-1].split("_")[2][1:])
            df.append((
                pos_idx,
                cluster_idx, 
                cst_score,
                iface_score
            ))
    
    df = pd.DataFrame(df)
    df.columns = [
        "position", 
        "cluster",
        "cst_score",
        "iface_score"
    ]
    score_stats_file = os.path.join(output_path, "scan_statistics.csv")
    df.to_csv(score_stats_file, index=False)
    return score_stats_file


def filter_hits(scan_csv_file, pdb_file, chain="A", 
                top_percentile=25, n_samples=10,
                out_prefix="matcher_selected_positions"):
    """
    Filter a requested top percentile of the matcher solutions according
    to protein-protractor interface Rosetta energy and the match constraint
    score. Good solutions are ones that satisfy the top percentile cutoff
    for both these scores. Several good solutions may correspond to the same
    residue, so further randomly subsample a given number of good solutions
    per residue. Cluster ids of matched pdb files corresponding to these 
    good solutions are recorded.  
    
    Args:
        scan_csv_file (str): CSV file containing the final scoring of all
        matcher solutions after Rosetta relaxation.
        
        pdb_file (str): (Relaxed) PDB file.
        
        chain (str, optional): Chain along which acylation positions were
        scanned. Defaults to "A".
        
        top_percentile (float, optional): Top percentile of matcher solutions to filter
        along both axes ie. interface score and matcher constraint score.
        
        n_samples (int, optional): Max. number of conformers to randomly sample for a single
        acylation location. If the actual number of conformers corresponding to good solutions
        for that residue is lower, pick the lower value.
        
        out_prefix (str, optional): Output prefix. This will be used to name:
        - csv file containing matched residue ids and indexes to pdbs produced from 
        Rosetta relaxation of matcher solutions. 
        
        - png file containing a plot of the matcher solution landscape.
        Defaults to "matcher_selected_positions".
    """
    
    # 1. remove raw data beyond baseline cutoffs
    df = pd.read_csv(scan_csv_file)
    df = df[
        (df["cst_score"] < MATCH_CST_SCORE_CUTOFF) & \
        (df["iface_score"] < INTERFACE_SCORE_CUTOFF)
    ]
    
    # 2. filter viable region according to top percentiles
    X_field, Y_field = "iface_score", "cst_score"
    X_bound = np.percentile(df[X_field].values, q=top_percentile)
    Y_bound = np.percentile(df[Y_field].values, q=top_percentile)
    
    select = []
    for (X, Y) in zip(df[X_field].values, df[Y_field].values):
        decision = int((X <= X_bound and Y <= Y_bound))
        select.append(decision)
    df = df.assign(select=select)
    
    # 3. plot the matcher solution landscape
    plt.style.use("ggplot")
    jg = sns.jointplot(
        data=df,
        x=X_field, y=Y_field, hue="select",
        kind="scatter", edgecolor=None, palette="crest"
    )
    
    for ax in [jg.ax_joint, jg.ax_marg_x]:
        ax.axvline(x=X_bound, color="black", lw=2)
    
    for ax in [jg.ax_joint, jg.ax_marg_y]:
        ax.axhline(y=Y_bound, color="black", lw=2)
    
    jg.figure.savefig(out_prefix + ".png", bbox_inches="tight", dpi=200)
    
    # 4. select a residue if it has at least one conformation in the viable region
    model = PDBParser(QUIET=True).get_structure("x", pdb_file)[0]
    residues = [r for r in model[chain].get_residues() if r.id[0] == " "]
    residues = sorted(residues, key=lambda r : r.id[1])
    res_map = {r.id[1]: r for r in residues}
    sel_map = {r.id[1]: 0 for r in residues}
    
    for i in range(len(df)):
        pos = df.iloc[i]["position"]
        sel_map[pos] += df.iloc[i]["select"]
    
    df_res = []
    for resid in sorted(sel_map):
        resname = seq3to1(res_map[resid].resname.upper())
        if int(sel_map[resid] >= NUM_CONFORMERS_CUTOFF):
            df_res.append((resid, resname, sel_map[resid]))
    df_res = pd.DataFrame(df_res, columns=["resid", "aa", "n_matches"])
    
    # 4. for each selected residue, extract <num_samples> random clusters
    cluster_ids = []
    for i in range(len(df_res)):
        resid = df_res.iloc[i]["resid"]
        this_df = df[(df["position"] == resid) & (df["select"] == 1)]
        this_df_samples = this_df.sample(n=min(n_samples, len(this_df)))
        this_cluster_ids = list(this_df_samples["cluster"].values)
        cluster_ids.append(this_cluster_ids)
    df_res = df_res.assign(cluster_ids=cluster_ids)
    
    df_res.to_csv(out_prefix + ".csv", index=False)
    