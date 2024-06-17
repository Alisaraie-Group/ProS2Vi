import argparse
from visual import VisualMap

def main():
    parser = argparse.ArgumentParser(description='A script that takes in the PDB code and path to the PDB file, and creates a visualization of the secondary structure assignments. Note that you need to have DSSP installed in your system, for this script to work.')

    # Arguments
    parser.add_argument('pdb_name', type=str, help='The PDB code of the file you want to parse.')
    parser.add_argument('pdb_file_path', type=str, help='The path to the PDB file.')
    parser.add_argument('pdb_file_type', type=str, help='The file type of the protein structure.')
    parser.add_argument('-r', dest='residues_per_line', type=int, default=50, help='The number of residues per each line. If no arguments is provided, it defaults to 50 residues per line.')
    parser.add_argument('-o', dest='output_image_name', type=str, default='', help='The name of the output image, you can use the png, jpg or gif extensions. If no argument is provided, it defaults to PDB_CODE.jpg')
    parser.add_argument('-d', dest='dpi', type=int, default=100, help='The DPI of the output image, as an integer. If no argument provided, it defaults to 100.')
    parser.add_argument('-pdf', action='store_true', default=False, help='A boolean flag indicating whether a PDF output is needed. If --pdf argument is sent, a PDF is generated.')

    args = parser.parse_args()

    vs = VisualMap(pdb_name=args.pdb_name, file_path=args.pdb_file_path, file_type=args.pdb_file_type)
    vs.generate_visual(residues_per_line=args.residues_per_line, output_image_name=args.output_image_name, dpi=args.dpi, pdf=args.pdf)



if __name__ == '__main__':
    main()