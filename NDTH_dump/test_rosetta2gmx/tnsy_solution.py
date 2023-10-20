import sys
sys.path.insert(0, "/nfs_home/users/tnsy/Projects/method_dev/nnprotscan")

from nnprotscan.mdprep import convert_to_gmx

convert_to_gmx("input.pdb", "P", "pyRED_amber03", out_pdb_file="output_2.pdb")

