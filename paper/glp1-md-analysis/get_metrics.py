import os
import sys
import pickle
import numpy as np

import mdmetric

# constants
RELOAD = False
N_BLOCKS = 4
DISTANCE_CUTOFF = 0.45 # nm
GMX_OUTPUT_PATH = os.path.expanduser("path_to_glp1_md_results")

# user input
acylation_resid = int(sys.argv[1])
OUTPUT_PATH = os.path.abspath(sys.argv[2])

# create output path
os.system(f"mkdir -p {OUTPUT_PATH}")


# get traj info
def _get_traj_info(acylation_resid, max_n_conformers=10):
    this_gmx_path = os.path.join(GMX_OUTPUT_PATH, f"K{acylation_resid}")
    
    traj_paths = []
    traj_files = []
    pdb_files = []
    out_prefixes = []
    
    for i in range(max_n_conformers+1):
        this_traj_path = os.path.join(this_gmx_path, f"conf_{i+1}")
        this_traj_file = os.path.join(this_gmx_path, f"conf_{i+1}", "step3_production_prot_memb.xtc")
        this_pdb_file = os.path.join(this_gmx_path, f"conf_{i+1}", f"K{acylation_resid}_conf_{i+1}_protein_memb.pdb")
        this_out_prefix = os.path.join(OUTPUT_PATH, f"K{acylation_resid}", f"conf{i+1}")
        
        if os.path.isfile(this_traj_file) and os.path.isfile(this_pdb_file):
            traj_paths.append(this_traj_path)
            traj_files.append(this_traj_file)
            pdb_files.append(this_pdb_file)
            out_prefixes.append(this_out_prefix)


    return list(zip(traj_paths, traj_files, pdb_files, out_prefixes))


# get indices corresponding to different sections ie. peptide, receptor, protractor, membrane, etc
def _get_indices(topology):
    indices = dict()
    heavy_indices = topology.select("((symbol != H) or (symbol != NA) or (symbol != CL))")
    
    # protractor
    indices_ = topology.select("(resname == SEM) and ((name != N) and (name != CA) and (name != C) and (name != O))")
    indices["protractor"] = [i for i in indices_ if i in heavy_indices]
    
    # peptide
    indices["peptide"] = [i for i in topology.select("chainid == 0") if (i in heavy_indices) and (i not in indices["protractor"])]
    indices["peptide_bb"] = [i for i in topology.select("(chainid == 0) and (backbone)") if i in heavy_indices]
    
    # peptide indices mapped to residue
    peptide_residues = [r for r in topology.residues if r.chain.index == 0]
    indices["peptide_res_dict"] = dict()
    for r in peptide_residues:
        indices["peptide_res_dict"][r.resSeq] = [a.index for a in r.atoms if a.name in ["N", "CA", "C", "O"]]
    
    # receptor
    indices["receptor"] = [i for i in topology.select("chainid == 1") if i in heavy_indices]
    indices["receptor_bb"] = [i for i in topology.select("(chainid == 1) and (backbone)") if i in heavy_indices]
    
    # membrane
    indices["membrane"] = [i for i in topology.select("resname == DPPC") if i in heavy_indices]
    indices["membrane_head_group"] = [i for i in topology.select("resname == DPPC and name == P")]
    return indices


# calculate metrics for a single traj (corresponding to a given conformation)
def _get_single_traj_metrics(trj, ref_trj, indices):
    metrics = dict()
    
    # membrane thickness
    print("> membrane thickness")
    metrics["membrane_thickness"] = mdmetric.get_membrane_thickness(trj, indices["membrane_head_group"])
    
    # protractor rg
    print("> protractor Rg")
    metrics["protractor_rg"] = mdmetric.get_protractor_Rg(trj, indices["protractor"])
    
    # protractor ree
    print("> protractor Ree")
    metrics["protractor_ree"] = mdmetric.get_protractor_Ree(trj, indices["protractor"])
    
    # protractor-peptide contacts
    print("> protractor-peptide contacts")
    metrics["protractor_peptide_contacts"] = mdmetric.get_protractor_contacts(
        trj, indices["protractor"], indices["peptide"], cutoff=DISTANCE_CUTOFF
    )
    
    # protractor-receptor contacts
    print("> protractor-receptor contacts")
    metrics["protractor_receptor_contacts"] = mdmetric.get_protractor_contacts(
        trj, indices["protractor"], indices["receptor"], cutoff=DISTANCE_CUTOFF
    )
    
    # protractor-membrane contacts
    metrics["protractor_membrane_contacts"] = mdmetric.get_protractor_contacts(
        trj, indices["protractor"], indices["membrane"], cutoff=DISTANCE_CUTOFF
    )
    
    # peptide backbone rmsd
    print("> peptide RMSD")
    metrics["peptide_rmsd"] = mdmetric.get_binder_rmsd(
        trj, ref_trj, indices["peptide_bb"], indices["receptor_bb"]
    )
    
    # peptide per residue rmsd
    print("> peptide per residue RMSD")
    metrics["peptide_per_residue_rmsd"] = mdmetric.get_binder_per_residue_rmsd(
            trj, ref_trj, indices["peptide_res_dict"], indices["receptor_bb"]
    )
    
    return metrics


def _get_single_traj_energies(traj_path):
    print("> energies")
    metric_file_map = {
        "protractor_peptide_coul_energy": "Coul_SR_GLP1_SEM",
        "protractor_receptor_coul_energy": "Coul_SR_GLP1R_SEM",
        "protractor_membrane_coul_energy": "Coul_SR_MEMB_SEM",
        "protractor_peptide_lj_energy": "LJ_SR_GLP1_SEM",
        "protractor_receptor_lj_energy": "LJ_SR_GLP1R_SEM",
        "protractor_membrane_lj_energy": "LJ_SR_MEMB_SEM",
    }
    
    metrics = dict()
    for k, v in metric_file_map.items():
        filename = os.path.join(traj_path, v + ".npz")
        if os.path.isfile(filename):
            data = np.load(os.path.join(traj_path, v + ".npz"))
            metrics[k] = data[v]
        else:
            metrics[k] = []
    return metrics


def _average_metrics(metrics_list):
    metrics = dict()
    
    # 1D quantities
    keys_1D = [
        "membrane_thickness", 
        "protractor_rg", 
        "protractor_ree", 
        "protractor_peptide_contacts",
        "protractor_receptor_contacts",
        "protractor_membrane_contacts",
        "peptide_rmsd",
        "protractor_peptide_coul_energy",
        "protractor_membrane_coul_energy",
        "protractor_receptor_coul_energy",
        "protractor_peptide_lj_energy",
        "protractor_membrane_lj_energy",
        "protractor_receptor_lj_energy"
    ]
    
    for k in keys_1D:
        x = np.hstack([m[k] for m in metrics_list])
        metrics[k] = mdmetric.block_average(x, nblocks=N_BLOCKS)
    
    # 2D metric, i.e RMSF
    x = np.vstack([m["peptide_per_residue_rmsd"] for m in metrics_list])
    metrics["peptide_rmsf"] = np.std(x, axis=0, ddof=1)
    
    return metrics


#### MAIN #### 
print("getting trajectory information")
traj_info = _get_traj_info(acylation_resid)
    
print("reading trajectory and calculating metrics...")
metrics_list = []
resids = None
    
for i in range(len(traj_info)):
    print("conformation ", i+1)
    
    # get trajectory info
    traj_path, traj_file, pdb_file, out_prefix = traj_info[i]
    out_path = os.path.dirname(out_prefix)
    if not os.path.isdir(out_path):
        os.makedirs(out_path, exist_ok=True) 
    
    # read traj
    print("> reading trajectory...")
    trj, ref_trj = mdmetric.read_traj(traj_file, pdb_file, out_prefix, reload=RELOAD)
    indices = _get_indices(trj.topology)
        
    # calculating metrics and energies
    this_metrics = _get_single_traj_metrics(trj, ref_trj, indices)
    this_metrics.update(_get_single_traj_energies(traj_path))   
    metrics_list.append(this_metrics)
    if resids is None:
        resids = sorted(indices["peptide_res_dict"])

# avg metrics
print("> averaging...")
metrics = _average_metrics(metrics_list)
metrics["resids"] = resids

# save metrics to file
metrics_file = os.path.join(OUTPUT_PATH, f"K{acylation_resid}", "metrics.pkl")
with open(metrics_file, "wb") as of:
    pickle.dump(metrics, of)
    

#resids = [14, 15, 16, 18, 19, 20, 22, 23, 26, 27, 30, 31, 33, 34, 35, 36]