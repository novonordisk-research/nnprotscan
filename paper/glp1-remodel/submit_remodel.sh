#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -o remodel.log
#SBATCH --job-name=remodel

path_to_rosetta/bin/remodel.static.linuxgccrelease \
-s 6x18_chR.pdb \
-remodel:hb_lrbb 2.0 \
-remodel:ss_pair 2.0 \
-remodel:rsigma 2.0 \
-remodel::quick_and_dirty \
-remodel:blueprint 6x18_blueprint_chR \
-remodel::num_trajectory 5 \
-nstruct 20 \