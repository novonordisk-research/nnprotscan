import os
import pickle
import itertools
import numpy as np

import mdtraj
from numba import njit, prange, types

def _superpose(trj, ref_trj, align_indices):
    ref_trj = ref_trj.center_coordinates()
    trj = trj.superpose(ref_trj, atom_indices=align_indices)
    return trj, ref_trj


def block_average(samples, nblocks=4):
    if len(samples) == 1:
        samples = samples.rehape(-1,1)
    nframes = samples.shape[0]
    blocksize = int(nframes / nblocks)
    
    # avg within blocks
    samples_block = np.zeros((nblocks, *samples.shape[1:]))
    for i in range(nblocks):
        start = i * blocksize
        stop = (i+1) * blocksize if i < nblocks - 1 else nframes
        samples_block[i] = np.mean(samples[start:stop], axis=0)
    
    # avg across blocks
    samples_mean = np.mean(samples_block, axis=0)
    samples_err = np.std(samples_block, axis=0)
    return samples_mean, samples_err


def read_traj(traj_file, pdb_file, out_prefix, nframes=None, 
              has_water=False, reload=False):
    
    traj_pkl_file = out_prefix + ".trj.pkl"
    ref_traj_pkl_file = out_prefix + ".reftrj.pkl"
    
    if not (os.path.isfile(traj_pkl_file) and os.path.isfile(ref_traj_pkl_file)) or reload:
        trj = mdtraj.load(traj_file, top=pdb_file)
        ref_trj = mdtraj.load(pdb_file)[0]
        
        # clean traj if necessary
        if nframes is not None:
            trj = trj[-nframes:1]
        if has_water:
            indices = trj.topology.select("not ((water) or (name == NA) or (name == CL))")
            trj = trj.atom_slice(indices)
            ref_trj = ref_trj.atom_slice(indices)
        
        with open(traj_pkl_file, "wb") as of:
            pickle.dump(trj, of)
        
        with open(ref_traj_pkl_file, "wb") as of:
            pickle.dump(ref_trj, of)
    
    else:
        with open(traj_pkl_file, "rb") as of:
            trj = pickle.load(of)
        
        with open(ref_traj_pkl_file, "rb") as of:
            ref_trj = pickle.load(of)
            
    return trj, ref_trj


def get_protractor_Rg(trj, indices):
    return mdtraj.compute_rg(trj.atom_slice(indices))

def get_protractor_Ree(trj, indices):
    pairs = np.array([[min(indices), max(indices)]])
    return mdtraj.compute_distances(trj, pairs).flatten()


def get_protractor_contacts(trj, indices, indices_other, cutoff=0.45):
    pairs = list(itertools.product(indices, indices_other))
    pair_distances = mdtraj.compute_distances(trj, pairs)
    contacts = np.sum(pair_distances <= cutoff, axis=1)
    return contacts


def get_binder_rmsd(trj, ref_trj, indices, align_indices):
    # align based on the receptor
    trj, ref_trj = _superpose(trj, ref_trj, align_indices)
    # compute the un-aligned RMSD (i.e. the dRMS)
    xyz = trj.atom_slice(indices).xyz
    xyz0 = ref_trj.atom_slice(indices).xyz
    rsq = np.sum( (xyz - xyz0) * (xyz - xyz0), axis=-1)
    rmsds = np.sqrt(np.mean(rsq, axis=-1))
    return rmsds


def get_binder_per_residue_rmsd(trj, ref_trj, res_indices_dict, align_indices):
    # align based on the receptor
    trj, ref_trj = _superpose(trj, ref_trj, align_indices)
    # get avg dRMS for each residue
    rmsds = []
    for r in sorted(res_indices_dict):
        indices = res_indices_dict[r]
        xyz = trj.atom_slice(indices).xyz
        xyz0 = ref_trj.atom_slice(indices).xyz
        this_msd = np.sum((xyz - xyz0) * (xyz - xyz0), axis=-1)
        this_rmsd = np.sqrt(np.mean(this_msd, axis=-1))
        rmsds.append(this_rmsd)
    return np.array(rmsds).transpose() # nframes X nres
    

@njit(types.Array(types.float32, 1, "C")(types.Array(types.float32, 2, "A"), types.Array(types.float32, 1, "A")),
    parallel=True,
    fastmath=True,
    nogil=True
)
def _get_membrane_thickness(
    z: types.Array(types.float32, 2, "A"),
    boxz: types.Array(types.float32, 1, "A")
) -> types.Array(types.float32, 1, "C"):
    
    # minimum image
    nframes, nlipids = z.shape
    zmin = np.zeros((nframes, nlipids), dtype=np.float32)
    for i in prange(nframes):
        this_boxz = boxz[i]
        this_inv_boxz = np.float32(1.0) / this_boxz
        for j in prange(nlipids):
            zmin[i,j] = z[i,j] - this_boxz * np.floor(z[i,j] * this_inv_boxz)
            
    # find center of membrane
    z_center = np.zeros(nframes, dtype=np.float32)
    for i in prange(nframes):
        for j in prange(nlipids):
            z_center[i] += zmin[i,j]
    z_center /= nlipids
    
    # calculate thickness
    thickness = np.zeros(nframes, dtype=np.float32)
    for i in prange(nframes):
        z_plus, z_minus = 0.0, 0.0
        n_plus, n_minus = 0, 0
        for j in prange(nlipids):
            if zmin[i,j] > z_center[i]:
                z_plus += zmin[i,j]
                n_plus += 1
            elif zmin[i,j] < z_center[i]:
                z_minus += zmin[i,j]
                n_minus += 1
        z_plus /= n_plus
        z_minus /= n_minus
        thickness[i] = abs(z_plus - z_minus)
    return thickness


def get_membrane_thickness(trj, indices):
    z_lipid = trj.atom_slice(indices).xyz[:, :, -1]
    boxz = trj.unitcell_lengths[:, -1]
    return _get_membrane_thickness(z_lipid, boxz)

