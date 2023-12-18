#!/bin/bash

###############################################################################################################################################################
## Author: Nidhin Thomas (NDTH)
## This script organizes the PDB outputs from Rosetta Match Protocol and runs MD simulations on all selected structures.

## bash protein_preparation.sh path/to/target/directory path/to/pdb/directory starting_atom_ID/for/mutating/residue ending_atom_ID/for/mutating/residue
###############################################################################################################################################################

root_dir="Path to the directory you want to run MD simulations"
target_dir=$1
PDB_dir=$2
atomID_start=$3
atomID_end=$4

cd ${PDB_dir}

for input_PDB in *.pdb;
do
pdbname=${input_PDB%.*}
residue=$(echo $pdbname | awk -F'_' '{ print $2 }')
conformation=$(echo $pdbname | awk -F'_' '{ print $4 }')

mkdir -p ${target_dir}/${residue}/conf_${conformation}

cp ${input_PDB}  ${target_dir}/${residue}/conf_${conformation}/
cp ${root_dir}/inputs/* -rf ${target_dir}/${residue}/conf_${conformation}/

cd ${target_dir}/${residue}/conf_${conformation}/

## converting the Rosetta PDB file to GMX version
python ${root_dir}/Rosetta2GMX.py ${input_PDB} ${root_dir}/SEM.pdb ./Rosetta2GMX_output.pdb ${atomID_start} ${atomID_end}

## GROMACS pdb2gmx to correctly represent the protonation state and hydrogen atom naming.
echo 7 | gmx_mpi pdb2gmx -f ./Rosetta2GMX_output.pdb -ff amber03 -ignh -o ${pdbname}_GMX.pdb

## Aligning conformations to protein embedded into the membrane. 
python ${root_dir}/protein_alignment.py ${root_dir}/GLP1_GLP1R_reference.pdb ${pdbname}_GMX.pdb ${pdbname}_aligned.pdb

cat ${pdbname}_aligned.pdb ${root_dir}/GLP1_GLP1R_GMX_membrane.pdb >> ${pdbname}_converted.pdb
temp_file="temp.pdb"

## Insert the new line at the beginning of the file
sed '1iCRYST1  115.734  135.217  150.000  90.00  90.00  90.00               1' ${pdbname}_converted.pdb > ${temp_file}

## Replace the original file with the temporary one
mv ${temp_file} ${pdbname}_converted.pdb

## Editing the posre.itp file replace restraint of 1000 into a variable that can be controlled during equilibration
sed -i 's/1000  1000  1000/POSRES_PROT  POSRES_PROT  POSRES_PROT/g' posre_Protein_chain_P.itp
sed -i 's/1000  1000  1000/POSRES_PROT  POSRES_PROT  POSRES_PROT/g' posre_Protein_chain_R.itp

cp ./itp/* ./

## renumbered the atoms to ensure that the PDB is correct.
gmx_mpi editconf -f ${pdbname}_converted.pdb -resnr -1 -o ${pdbname}_resnr.pdb

## Adding ions to neutralize and replicate physiological buffer concentration. 
gmx_mpi grompp -f ./mdps/genion.mdp -o ${pdbname}_genion.tpr -c ${pdbname}_resnr.pdb -r ${pdbname}_resnr.pdb -maxwarn -1

## Adding ions to the system. Here, NA and CL ions are added. topol.top is automatically upated. 
echo "17" | gmx_mpi genion -s ${pdbname}_genion.tpr -p topol.top -pname Na -nname Cl -neutral -o ${pdbname}_ion.pdb

## Check the files by preparing the input files for MD
gmx_mpi grompp -f ./mdps/step2.0_minimization.mdp -c ${pdbname}_ion.pdb -r ${pdbname}_ion.pdb -maxwarn -1 -o step2.0_minimization.tpr

## Run initial minimization
#gmx_mpi mdrun -deffnm step2.0_minimization -v

## Creating index file
echo -e "1\nname 24 PROT\n20\nname 25 MEMB\n23\nname 26 SOLV\nq" | gmx_mpi make_ndx -f step2.0_minimization.tpr -o index.ndx

## Production MD with minimization and equilibration
sbatch sim_prod.sh

cd ${PDB_dir}

done
