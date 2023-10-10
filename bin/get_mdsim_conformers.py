import os
import glob
import argparse
import pandas as pd

from nnprotscan import mdprep
from nnprotscan.linkage_model import LinkageModel

# supress pesky biopython warnings
import warnings
warnings.simplefilter ("ignore")

_PH = 7.0

# user args
parser = argparse.ArgumentParser(description="Get initial conformers for MD simulations from Rosetta Matcher outputs.")

parser.add_argument("matcher_path", help="Output path containing results of the matcher pipeline.")
parser.add_argument("-c", "--chain", help="Binder chain id in the matcher PDB files.")
parser.add_argument("-l", "--linkage", default="amide", help="Type of acylation.")
parser.add_argument("-o", "--output_path", default="md_conf", help="Output path for writing PDB files.")
parser.add_argument("-gmx", "--to_gmx", action="store_true", help="If true convert to a gromacs compatible PDB file.")
parser.add_argument("-ff", "--forcefield", help="If requesting gromacs conversion, supply forcefield name.")

args = parser.parse_args()

# check if ff has been provided
if args.to_gmx and (args.forcefield is None):
    raise TypeError("Forcefield required if conversion to gromacs compatible PDB is requested.")

# create output path
output_path = os.path.abspath(args.output_path)
if not os.path.isdir(output_path):
    os.makedirs(output_path, exist_ok=True)

# paths / filenames etc
matcher_path = os.path.abspath(args.matcher_path)
matcher_pdb_path = os.path.join(matcher_path, "scores", "pdbs")
match_result_csv_file = os.path.join(matcher_path, "match_result.csv")
df = pd.read_csv(match_result_csv_file)

# extract each conformer for a given residue and convert to hybrid residue type PDB file
print("Writing conformers for:")
for i in range(len(df)):
    resid = df.iloc[i]["resid"]
    cluster_ids = eval(df.iloc[i]["cluster_ids"])
    print(f"> residue {resid}...")
    for j, c in enumerate(cluster_ids):
        # copy appropriate matcher output to destination
        match_res_name = LinkageModel(args.linkage).match_res
        
        glob_str = os.path.join(matcher_pdb_path, f"UM_{c}_{match_res_name}{resid}_*_1_0001.pdb")
        src_file = glob.glob(glob_str)
        assert len(src_file) == 1
        src_file = src_file[0]
        
        tar_file = os.path.join(output_path, f"matcher_{match_res_name}{resid}_conf_{j+1}.pdb")
        os.system(f"cp {src_file} {tar_file}")
        
        # build topology
        match_site = args.chain + str(resid)
        t = mdprep.TopologyBuilder(tar_file, match_site=match_site, linkage=args.linkage)
        t.protonate(pH=_PH)
        t.acylate()
        
        t.write(tar_file)
        
        if args.to_gmx:
            mdprep.convert_to_gmx(pdb_file=tar_file, chain_id=args.chain, ff=args.forcefield)
            