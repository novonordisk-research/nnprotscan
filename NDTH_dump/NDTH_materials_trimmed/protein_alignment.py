import pymol
from pymol import cmd
import argparse

def align_and_save(input1, input2, output):
    # Start PyMOL
    pymol.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
    pymol.finish_launching()

    # Load the first PDB
    cmd.load(input1, 'protein1')

    # Load the second PDB
    cmd.load(input2, 'protein2')

    # Align protein2 onto protein1
    cmd.align('protein2', 'protein1')
    
    cmd.set('retain_order', 1)

    # Save the aligned structure of protein2 only
    cmd.save(output, 'protein2')

    # Quit PyMOL
    cmd.quit()

if __name__ == "__main__":
    # Setup argparse
    parser = argparse.ArgumentParser(description='Align two PDB files and save the aligned structure of the second protein.')
    parser.add_argument('input1', help='Path to the first input PDB file (the structure to align to)')
    parser.add_argument('input2', help='Path to the second input PDB file (the structure to be aligned)')
    parser.add_argument('output', help='Path to the output PDB file (the aligned structure of the second protein)')

    args = parser.parse_args()

    # Call the function
    align_and_save(args.input1, args.input2, args.output)



