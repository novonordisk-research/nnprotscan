import argparse

#########################################################################################################################################################################
## This is a custom script used to convert the Rosetta output PDB file into GROMACS readable PDB file. The purpose of the script is to rename and reorder the NCA 
## residue into gromacs format SEM. SEM residue is custom built using QM method and the force field is generated for AMBER03 format. 

##Example: python Rosetta2GMX_full.py <input.pdb> SEM.pdb <output.pdb>
#########################################################################################################################################################################

def rename_atoms(pdb_file, ref_file, output,SEM_starting_atomID, SEM_ending_atomID):
    atom_mapping_dict = {
        "N": "N",
        "CA": "CA",
        "C": "C",
        "O": "O",
        "OXT": "OXT",
        "CB": "CB",
        "CG": "CG",
        "CD": "CD",
        "CE": "CE",
        "NZ": "NZ",
        "H": "HN",
        "HA": "HA",
        "HB2": "HB1",
        "HB3": "HB2",
        "HG2": "HG1",
        "HG3": "HG2",
        "HD2": "HD1",
        "HD3": "HD2",
        "HE2": "HE1",
        "HE3": "HE2",
        "HZ3": "HZ1",
        "C15x": "C24",
        "C14x": "C23",
        "C13x": "C21",
        "N2x": "N20",
        "C12x": "C19",
        "C11x": "C18",
        "O6x": "O17",
        "C10x": "C16",
        "C9x": "C15",
        "O5x": "O14",
        "C8x": "C13",
        "C7x": "C11",
        "N1x": "N10",
        "C6x": "C9",
        "C5x": "C8",
        "O3x": "O7",
        "C4x": "C6",
        "C3x": "C5",
        "O2x": "O4",
        "C2x": "C3",
        "C1x": "C1",
        "O1x": "O2",
        "O4x": "O12",
        "O7x": "O22",
        "C16x": "C25",
        "N3x": "N29",
        "C17x": "C30",
        "O8x": "O31",
        "C18x": "C32",
        "C19x": "C33",
        "C20x": "C34",
        "C21x": "C35",
        "C22x": "C36",
        "C23x": "C37",
        "C24x": "C38",
        "C25x": "C39",
        "C26x": "C40",
        "C27x": "C41",
        "C28x": "C42",
        "C29x": "C43",
        "C30x": "C44",
        "C31x": "C45",
        "C32x": "C46",
        "C33x": "C47",
        "C34x": "C48",
        "O9x": "O49",
        "O10x": "O50",
        "C35x": "C26",
        "O11x": "O27",
        "O12x": "O28",
        "H26x": "H241",
        "H27x": "H242",
        "H24x": "H231",
        "H25x": "H232",
        "H23x": "H201",
        "H21x": "H191",
        "H22x": "H192",
        "H19x": "H181",
        "H20x": "H182",
        "H17x": "H161",
        "H18x": "H162",
        "H15x": "H151",
        "H16x": "H152",
        "H13x": "H131",
        "H14x": "H132",
        "H12x": "H101",
        "H10x": "H91",
        "H11x": "H92",
        "H8x": "H81",
        "H9x": "H82",
        "H6x": "H61",
        "H7x": "H62",
        "H4x": "H51",
        "H5x": "H52",
        "H2x": "H31",
        "H3x": "H32",
        "H28x": "H251",
        "H29x": "H291",
        "H30x": "H321",
        "H31x": "H322",
        "H32x": "H331",
        "H33x": "H332",
        "H34x": "H341",
        "H35x": "H342",
        "H36x": "H351",
        "H37x": "H352",
        "H38x": "H361",
        "H39x": "H362",
        "H40x": "H371",
        "H41x": "H372",
        "H42x": "H381",
        "H43x": "H382",
        "H44x": "H391",
        "H45x": "H392",
        "H46x": "H401",
        "H47x": "H402",
        "H48x": "H411",
        "H49x": "H412",
        "H50x": "H421",
        "H51x": "H422",
        "H52x": "H431",
        "H53x": "H432",
        "H54x": "H441",
        "H55x": "H442",
        "H56x": "H451",
        "H57x": "H452",
        "H58x": "H461",
        "H59x": "H462",
        "H60x": "H471",
        "H61x": "H472",
    }
    
    prefix_lines = []
    output_lines = []
    suffix_lines = []
    
    SEM_atomID = []
    ref_lines = []
    output_atom_number = []

    # Read the reference file and store the atom order
    with open(ref_file, 'r') as ref:
        for line in ref:
            if line.startswith('ATOM'):
                ref_lines.append(line)

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('REMARK') or line.startswith('END') or line.startswith('CONECT'):
                continue
            elif ( int(line[6:11].strip()) < SEM_starting_atomID ):
                prefix_lines.append(line)
            elif line.startswith('HETATM'):
                atom_name = line[11:16].strip()
                output_atom_number.append(line[6:11].strip())
                if atom_name in atom_mapping_dict.keys():
                    new_atom_name = atom_mapping_dict[atom_name]
                    line = "ATOM" + "{:>7s}".format(output_atom_number[-1]) + "{:>5s}".format(new_atom_name) + "{:>4s}".format("SEM") + line[20:]
                    SEM_atomID.append(int((line[6:11]).strip()))
                    output_lines.append(line)
            else:
                suffix_lines.append(line)

    # Rearrange the atoms based on the reference atom order
    reordered_lines = []
    for ref_line in ref_lines:
        ref_atom_name = ref_line[11:16].strip()
        for line in output_lines:
            atom_name = line[11:16].strip()
            if atom_name == ref_atom_name:
                reordered_lines.append(line)
                print(line)
                break
                
    # Write to output file
    with open(output, 'w') as file:
        file.writelines(prefix_lines + reordered_lines + suffix_lines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rename atoms in a PDB file based on a reference file.')
    parser.add_argument('input_file', help='Path to the input PDB file')
    parser.add_argument('ref_file', help='Path to the reference PDB file')
    parser.add_argument('output_file', help='Path to the output PDB file')
    parser.add_argument('atomID_start', help='Atom number of the first atom in NCA residue', type=int)
    parser.add_argument('atomID_end', help='Atom number of the last atom in NCA residue', type=int)
    args = parser.parse_args()

    rename_atoms(args.input_file, args.ref_file, args.output_file, args.atomID_start, args.atomID_end)
