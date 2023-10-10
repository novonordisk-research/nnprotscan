#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --output=glp1.log
#SBATCH --job-name=matcher
#SBATCH --mem=60G

#module load rosetta/2021
export ROSETTA_HOME=/nfs_home/projects/departments/cdd/programs/rosetta_src_2021.16.61629_bundle/main/source
export PYTHONPATH=$PYTHONPATH:$HOME/Projects/method_dev/nnprotscan

# conda environment
eval "$(conda shell.bash hook)"
conda activate nnprotscan

PDB_FILE=../data/glp1/6x18_glp1-7-36_glp1r-full_nowater_nocap.pdb
SMILES='O=CCOCCOCCNC(=O)COCCOCCNC(=O)CC[C@@H](NC(=O)CCCCCCCCCCCCCCCCC(=O)[O-])C(=O)[O-]'

python ../../bin/run_matcher_scan.py $PDB_FILE $SMILES -o outputs -c P -l amide -np 16

