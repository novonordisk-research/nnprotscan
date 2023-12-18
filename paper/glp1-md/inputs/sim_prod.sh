#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --job-name=mdrun
#SBATCH -p gpu
#SBATCH -n 32 -N 1
#SBATCH --gres=gpu:1
#SBATCH --mem 32G

module purge

module use --append /nfs_home/software/spack/environments/prod/shared/share/spack/modules/linux-ubuntu20.04-broadwell

module load gromacs/2021.4-gcc-9.4.0-pwff4g4

module load openmpi/4.1.4
module load cuda/11.7

input_top="topol.top"
input_mdp="./mdps"
input_ndx="index.ndx"
cores=32

nvidia-smi -L

echo "Starting minimization"

#gmx_mpi grompp -f ${input_mdp}/step2.0_minimization.mdp -o step2.0_minimization.tpr -c ${input_str} -r ${input_str} -p ${input_top} -maxwarn -1

gmx_mpi mdrun -deffnm step2.0_minimization -v

echo "Starting equilibration: stage 1"

gmx_mpi grompp -f ${input_mdp}/step2.1_equilibration.mdp -o step2.1_equilibration.tpr -c step2.0_minimization.gro -r step2.0_minimization.gro -p ${input_top} -n ${input_ndx} -maxwarn -1

gmx_mpi mdrun -pme gpu -pin on -pinoffset 0 -pinstride 1 -nb gpu -bonded gpu -v -deffnm step2.1_equilibration

echo "Starting equilibration: stage 2"

gmx_mpi grompp -f ${input_mdp}/step2.2_equilibration.mdp -o step2.2_equilibration.tpr -c step2.1_equilibration.gro -r step2.1_equilibration.gro -t step2.1_equilibration.cpt -p ${input_top} -n ${input_ndx} -maxwarn -1

gmx_mpi mdrun -pme gpu -pin on -pinoffset 0 -pinstride 1 -nb gpu -bonded gpu -v -deffnm step2.2_equilibration

echo "Starting equilibration: stage 3"

gmx_mpi grompp -f ${input_mdp}/step2.3_equilibration.mdp -o step2.3_equilibration.tpr -c step2.2_equilibration.gro -r step2.2_equilibration.gro -t step2.2_equilibration.cpt -p ${input_top} -n ${input_ndx} -maxwarn -1

gmx_mpi mdrun -pme gpu -pin on -pinoffset 0 -pinstride 1 -nb gpu -bonded gpu -v -deffnm step2.3_equilibration

echo "Starting equilibration: stage 4"

gmx_mpi grompp -f ${input_mdp}/step2.4_equilibration.mdp -o step2.4_equilibration.tpr -c step2.3_equilibration.gro -r step2.3_equilibration.gro -t step2.3_equilibration.cpt -p ${input_top} -n ${input_ndx} -maxwarn -1

gmx_mpi mdrun -pme gpu  -pin on -pinoffset 0 -pinstride 1 -nb gpu -bonded gpu -v -deffnm step2.4_equilibration

echo "Starting equilibration: stage 5"

gmx_mpi grompp -f ${input_mdp}/step2.5_equilibration.mdp -o step2.5_equilibration.tpr -c step2.4_equilibration.gro -r step2.4_equilibration.gro -t step2.4_equilibration.cpt -p ${input_top} -n ${input_ndx} -maxwarn -1

gmx_mpi mdrun -pme gpu  -pin on -pinoffset 0 -pinstride 1 -nb gpu -bonded gpu -v -deffnm step2.5_equilibration

echo "Starting equilibration: stage 6"

gmx_mpi grompp -f ${input_mdp}/step2.6_equilibration.mdp -o step2.6_equilibration.tpr -c step2.5_equilibration.gro -r step2.5_equilibration.gro -t step2.5_equilibration.cpt -p ${input_top} -n ${input_ndx} -maxwarn -1

gmx_mpi mdrun  -pme gpu  -pin on -pinoffset 0 -pinstride 1 -nb gpu -bonded gpu -v -deffnm step2.6_equilibration

echo "Finished equilibration: stage 6"

echo "Starting Production run"

gmx_mpi grompp -f ${input_mdp}/step3_production.mdp -o step3_production.tpr -c step2.6_equilibration.gro -r step2.6_equilibration.gro -t step2.6_equilibration.cpt -p ${input_top} -n ${input_ndx} -maxwarn -1

gmx_mpi mdrun -pme gpu  -pin on -pinoffset 0 -pinstride 1 -nb gpu -bonded gpu -v -deffnm step3_production
