#   Copyright 2024 Muhammad Luckman Qasim,  Laleh Alisaraie
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import argparse
from src import visual

def main():
    parser = argparse.ArgumentParser(description='A script that takes in the PDB code and path to the PDB file, and creates a visualization of the secondary structure assignments. Note that you need to have DSSP installed in your system, for this script to work.')

    # Arguments
    parser.add_argument('pdb_file_path', type=str, help='The path to the PDB file.')
    parser.add_argument('-n', dest='pdb_name', type=str, default=None, help='The PDB code of the file you want to parse.')
    parser.add_argument('-s', dest='subtitle', type=str, default=None, help='A short informative description of the protein.')
    parser.add_argument('-sn', dest='scientific_name', type=str, default=None, help='The source scientific name.')
    parser.add_argument('-r', dest='residues_per_line', type=int, default=50, help='The number of residues per each line. If no arguments is provided, it defaults to 50 residues per line.')
    parser.add_argument('-o', dest='output_image_name', type=str, default='', help='The name of the output image, you can use the png or jpg. If no argument is provided, it defaults to PDB_CODE.png')
    parser.add_argument('-d', dest='dpi', type=int, default=100, help='The DPI of the output image, as an integer. If no argument provided, it defaults to 100.')
    parser.add_argument('-pdf', action='store_true', default=False, help='A boolean flag indicating whether a PDF output is needed. If --pdf argument is sent, a PDF is generated.')

    args = parser.parse_args()

    vs = visual.VisualMap(file_path=args.pdb_file_path, pdb_name=args.pdb_name, subtitle=args.subtitle, scientific_name=args.scientific_name)
    vs.generate_visual(residues_per_line=args.residues_per_line, output_image_name=args.output_image_name, dpi=args.dpi, pdf=args.pdf)


if __name__ == '__main__':
    main()
