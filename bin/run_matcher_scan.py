import os
import time
import argparse
from nnprotscan import rosetta_matcher as matcher

_MATCHER_RMS_CUTOFF = 3.0 # A

# user args
parser = argparse.ArgumentParser(description="Rosetta matcher pipeline for protractor placement.")

parser.add_argument("pdb_file", help="PDB file containing binder-receptor co-complex structure.")
parser.add_argument("smiles", help="SMILES string of protractor molecule (usually fatty acid + linker).")
parser.add_argument("-c", "--chain", help="Binder chain in the PDB file.")
parser.add_argument("-l", "--linkage", default="amide", help="Type of acylation.")
parser.add_argument("-f", "--fixed_positions", default="", help="Comma separated list of positions that are fixed and hence not scanned.")
parser.add_argument("-n_conf", "--n_conformers", type=int, default=10, help="Number of initial protractor conformers to dock against protein.")
parser.add_argument("-p", "--top_percentile", type=float, default=25, help="Top percentile to filter along interface scores and match constraint scores.")
parser.add_argument("-n_samples", "--n_samples", type=int, default=10, help="Max. number of random viable protractor conformations generated per acylation residue.")
parser.add_argument("-o", "--output_path", default="matcher_output", help="Output path.")
parser.add_argument("-np", "--n_processors", type=int, default=os.cpu_count(), help="Number of parallel procs for matcher run.")

args = parser.parse_args()

# create output location
output_path = os.path.abspath(args.output_path)
if not os.path.isdir(output_path):
    os.makedirs(output_path, exist_ok=True)

t0 = time.time()

# generate sdf file for protractor
matcher.generate_small_mol_protractor_conformations(
    smiles=args.smiles,
    out_sdf_file=os.path.join(output_path, "protractor.sdf"),
    n=args.n_conformers
)

# generate small mol params for rosetta
matcher.prepare_linker_params(
    sdf_file=os.path.join(output_path, "protractor.sdf"),
    linkage=args.linkage,
    out_prefix=os.path.join(output_path, "protractor")
)

# relax the given protein structure
matcher.relax_protein(
    pdb_file=os.path.abspath(args.pdb_file),
    out_prefix=os.path.join(output_path, "protein")
)

# set fixed positions if provided
if not args.fixed_positions:
    fixed_positions = []
else:
    fixed_positions = [int(p) for p in args.fixed_positions.split(",")]

# run Rosetta Matcher protocol
matcher.run(
    pdb_file=os.path.join(output_path, "protein.pdb"),
    params_file=os.path.join(output_path, "protractor.params"),
    constraint_file=os.path.join(output_path, "protractor.cst"),
    chain=args.chain,
    fixed_positions=fixed_positions,
    rms_cutoff=_MATCHER_RMS_CUTOFF,
    n_cores=args.n_processors,
    output_path=os.path.join(output_path, "matches")
)

# score the matches
scan_csv_file = matcher.score(
    matched_pdb_path=os.path.join(output_path, "matches"),
    chain=args.chain,
    params_file=os.path.join(output_path, "protractor.params"),
    constraint_file=os.path.join(output_path, "protractor.cst"),
    output_path=os.path.join(output_path, "scores"),
    n_cores=args.n_processors
)

# cluster the hits
matcher.filter_hits(
    scan_csv_file=scan_csv_file,
    pdb_file=os.path.join(output_path, "protein.pdb"),
    chain=args.chain,
    top_percentile=args.top_percentile,
    n_samples=args.n_samples,
    out_prefix=os.path.join(output_path, "match_result")
)

dt = time.time() - t0
print("ELAPSED TIME=", dt)