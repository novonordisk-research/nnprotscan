#!/bin/bash

###############################################################################################################################################################
## Author: Nidhin Thomas (NDTH)
## This script organizes the PDB outputs from Rosetta Match Protocol and runs MD simulations on all selected structures.

## bash protein_preparation_c_term.sh path/to/target/directory path/to/pdb/directory starting_atom_ID/for/mutating/residue ending_atom_ID/for/mutating/residue
###############################################################################################################################################################

root_dir="Path to project directory"
target_dir=$1
PDB_dir=$2
atomID_start=$3
atomID_end=$4
input_mdps="./mdps/"

cd ${PDB_dir}

for input_PDB in *.pdb;
do
pdbname=${input_PDB%.*}
residue=$(echo $pdbname | awk -F'_' '{ print $2 }')
conformation=$(echo $pdbname | awk -F'_' '{ print $4 }')

mkdir -p ${target_dir}/${residue}/conf_${conformation}

cp ${input_PDB}  ${target_dir}/${residue}/conf_${conformation}/
cp ${root_dir}/inputs_c_term/* -rf ${target_dir}/${residue}/conf_${conformation}/

cd ${target_dir}/${residue}/conf_${conformation}/

## converting the Rosetta PDB file to GMX version
python ${root_dir}/Rosetta2GMX_c_term.py ${input_PDB} ${root_dir}/CSEM.pdb ./Rosetta2GMX_output.pdb ${atomID_start} ${atomID_end}

## GROMACS pdb2gmx to correctly represent the protonation state and hydrogen atom naming.
echo 7 | gmx_mpi pdb2gmx -f ./Rosetta2GMX_output.pdb -ff amber03 -ignh -o ${pdbname}_GMX.pdb -water tip3p -p topol.top

## Editing the posre.itp file replace restraint of 1000 into a variable that can be controlled during equilibration
sed -i 's/1000  1000  1000/POSRES_PROT  POSRES_PROT  POSRES_PROT/g' posre_Protein_chain_A.itp
sed -i 's/1000  1000  1000/POSRES_PROT  POSRES_PROT  POSRES_PROT/g' posre_Protein_chain_B.itp

## renumbered the atoms to ensure that the PDB is correct.
gmx_mpi editconf -f ${pdbname}_GMX.pdb -resnr -1 -o ${pdbname}_resnr.pdb

## Next step is to define the system box size. Typically, cubic box is used. It can be dodecahedron, triclinic or octahedron.
gmx_mpi editconf -f ${pdbname}_resnr.pdb -bt cubic -d 1.0 -c -o protein_box.gro

## Solvating the protein. The size of the simulation box should be defined prior to solvating the protein. 
gmx_mpi solvate -cp protein_box.gro -cs spc216.gro -p topol.top -o protein_solv.gro

## Adding ions to neutralize and replicate physiological buffer concentration. 
gmx_mpi grompp -f ${input_mdps}/genion.mdp -o protein_genion.tpr -c protein_solv.gro -r protein_solv.gro -maxwarn 1

## Creating index groups so that water/SOL molecules can be replaced by ions. 
echo "q" | gmx_mpi make_ndx -f protein_genion.tpr -o index_genion.ndx

## Adding ions to the system. Here, SOD and CLA ions are added. topol.top is automatically upated. 
echo "13" | gmx_mpi genion -s protein_genion.tpr -n index_genion.ndx -p topol.top -pname NA -nname CL -conc 0.15 -neutral -o protein_ion.gro

echo "Starting energy minimization"

gmx_mpi grompp -f ${input_mdps}/step2.0_minimization.mdp -o step2.0_minimization.tpr -c protein_ion.gro -r protein_ion.gro -p topol.top

## Creating index file
echo -e "1\nname 19 PROT\n11\nname 20 SOLV\nq" | gmx_mpi make_ndx -f step2.0_minimization.tpr -o index.ndx

## Production MD with minimization and equilibration
sbatch sim_prod.sh

cd ${PDB_dir}

done
