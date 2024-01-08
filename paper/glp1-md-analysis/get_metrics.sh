#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --output=logs/mdmetrics_%a.log
#SBATCH --job-name=mdmetrics
#SBATCH --mem=60G
# SBATCH --array=14,15,16,18,19,20,22,23,26,27,30,31,33,34,35,36
#SBATCH --array=22,26

export PYTHONPATH=$PYTHONPATH:$HOME/path_to_nnprotscan/nnprotscan

# conda environment
eval "$(conda shell.bash hook)"
conda activate nnprotscan

python get_metrics.py ${SLURM_ARRAY_TASK_ID} .

