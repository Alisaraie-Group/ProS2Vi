#   Copyright 2024-2025 Muhammad Luckman Qasim,  Laleh Alisaraie
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

from jinja2 import Environment, FileSystemLoader
from markupsafe import Markup
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.DSSP import DSSP
import math
import pdfkit
import imgkit
import pdf2image
import requests
import os
from src.assets import ICONS


class VisualMap:
    '''
    Computes a Visualization mapping, residues with their secondary structure types.

    '''
    COLORS = {
        'H_COLOR': '#0000ff',
        'E_COLOR': '#ff0000',
        'E_A_COLOR': '#ff0000',
        'B_COLOR': '#228B22',
        'B_A_COLOR': '#228B22',
        'I_COLOR': '#ff00ff',
        'T_COLOR': '#ffff00',
        'S_COLOR': '#ffa500',
        'G_COLOR': '#ff0000',
        '-_COLOR': '#b7b7b7',
        'P_COLOR': '#228B22'
    }
    
    def __init__(self, file_path: str, pdb_name: str = None, subtitle: str = None, scientific_name: str = None) -> None:
        '''
        
        Args:
            pdb_name (str): The PDB code of the file.
            file_path (str): The file path of the PDB file.

        '''
        self.structure_list = self._get_dssp_output(pdb_name, file_path)
        self.file_path = file_path
        self.pdb_name = pdb_name
        self.subtitle = subtitle
        self.scientific_name = scientific_name

    def _update_template(self, residues_per_line: int = 50) -> str:
        '''
        A private Method that updates the template and returns a string representing the updated template.

        '''
        environment = Environment(loader=FileSystemLoader("templates/"))
        template = environment.get_template("template.html.jinja")
        unitprot_data = self.get_uniprot_data(self.pdb_name)

        residue_tables = ''

        for chain_id, chain in self.structure_list.items():
            table = ''
            struc_annotation_row = ''
            structure_row = ''
            residue_row = ''
            last_count = 0
            curr_structure_start = 0
            index_count = {'H': 0, 'B': 0, 'E': 0, 'I': 0, 'G': 0, 'S': 0, 'T': 0}

            for i, res in enumerate(chain):
                # Get the structure icon
                if (i == (len(chain) - 1) or chain[i]['res_struc'] != chain[i+1]['res_struc']) and chain[i]['res_struc'] in ['B', 'E']:
                    svg_string = ICONS[f"{res['res_struc']}_A"]
                    color = VisualMap.COLORS[f"{res['res_struc']}_A_COLOR"]
                    structure_icon = self._inject_svg_color(svg_string, color)
                else:
                    svg_string = ICONS[res['res_struc']]
                    color = VisualMap.COLORS[f"{res['res_struc']}_COLOR"]
                    structure_icon = self._inject_svg_color(svg_string, color)

                # Create the table cells for the indexed secondary structure assignment
                if i == (len(chain) - 1) or chain[i]['res_struc'] != chain[i+1]['res_struc'] or (i + 1) % residues_per_line == 0:
                    col_span = i + 1 - curr_structure_start
                    curr_structure_start = i + 1
                    if chain[i]["res_struc"] in index_count.keys() and chain[i]['res_struc'] != chain[i+1]['res_struc']:
                        index_count[chain[i]['res_struc']] += 1
                        struc_annotation_row += f'<td id="annotation" colspan="{col_span}">{chain[i]["res_struc"]}{index_count[chain[i]["res_struc"]]}</td>' 
                    else:
                        struc_annotation_row += f'<td id="annotation" colspan="{col_span}"></td>'
                structure_row += f'<td id="icon_row">{structure_icon}</td>'
                residue_row += f'<td id="res_row">{res["res_name"]}</td>'

                if (i + 1) % residues_per_line == 0 or i == (len(chain) - 1):
                    count = i
                    while (count + 1) % residues_per_line != 0:
                        structure_row += f'<td></td>'
                        residue_row += f'<td></td>'
                        struc_annotation_row += f'<td id="annotation"></td>'
                        count += 1
                    table += f'<tr id="annotation_row"><td></td>{struc_annotation_row}<td></td></tr>'
                    table += f'<tr><td id="icon_row"></td>{structure_row}<td></td></tr>'
                    start_res_num = chain[last_count]['res_num']
                    end_res_num = chain[i]['res_num']
                    table += f'<tr><td id="count_element_start">{start_res_num}</td>{residue_row}<td id="count_element_end">{end_res_num}</td></tr>'
                    struc_annotation_row = ''
                    structure_row = ''
                    residue_row = ''
                    last_count = i + 1
            
            uniprot_id = self.get_uniprot_id_by_chain_id(unitprot_data, chain_id)
            residue_tables += f'<div id="chain">Chain {chain_id}:</div>'
            if uniprot_id:
                residue_tables += f'<div>Uniprot ID: {uniprot_id}</div>'
            residue_tables += f'<table id="res"><tbody id="res">{table}</tbody></table>'
            last_count = 0
            index_count = {'H': 0, 'B': 0}

        rcsb_data = self.get_rcsb_entry_data(self.pdb_name)
        if self.pdb_name:
            self.pdb_name = self.pdb_name.upper()        
        if self.subtitle == None and rcsb_data:
            self.subtitle = rcsb_data['struct']['title']
        if self.scientific_name == None:
            self.scientific_name = self.get_scientific_name(self.pdb_name, rcsb_data)

        content = template.render(
            pdb_name = self.pdb_name,
            pdb_title = self.subtitle,
            scientific_name = self.scientific_name,
            residue_tables = residue_tables,
            H = self._inject_svg_color(ICONS['H'], VisualMap.COLORS['H_COLOR']),
            B = self._inject_svg_color(ICONS['B_A'], VisualMap.COLORS['B_A_COLOR']),
            E = self._inject_svg_color(ICONS['E_A'], VisualMap.COLORS['E_A_COLOR']),
            G = self._inject_svg_color(ICONS['G'], VisualMap.COLORS['G_COLOR']),
            I = self._inject_svg_color(ICONS['I'], VisualMap.COLORS['I_COLOR']),
            T = self._inject_svg_color(ICONS['T'], VisualMap.COLORS['T_COLOR']),
            S = self._inject_svg_color(ICONS['S'], VisualMap.COLORS['S_COLOR']),
            P = self._inject_svg_color(ICONS['P'], VisualMap.COLORS['P_COLOR']),
            U = self._inject_svg_color(ICONS['-'], VisualMap.COLORS['-_COLOR'])
        )
            
        return content
    
    def _get_width_and_height(self, residues_per_line):
        '''
        Private method that returns the width and height of the resulting visualization in the form of a table. Eg. (output_width, output_height)

        '''
        output_width = (20 * residues_per_line + 48 + 100 + 50) / 2
        output_height = (len(self.structure_list) * 76 + sum(math.ceil(len(i) / residues_per_line) for i in self.structure_list.values()) * 96 + 100 + 77 + 79 + 400) / 2
        return (int(output_width), int(output_height))
    
    @staticmethod
    def _inject_svg_color(svg_string, color):
        '''
        Helper function that injects a color into an SVG string and wraps it for inline display.
        
        Args:
            svg_string (str): The SVG template string with {color} placeholder
            color (str): The color code to inject (e.g., '#0000ff')
        
        Returns:
            Markup: The SVG string with color injected, marked as safe HTML for Jinja2 rendering
        '''
        # Inject the color into the SVG
        svg_with_color = svg_string.format(color=color)
        
        # Return as Markup to prevent HTML escaping in Jinja2 template
        # Strip whitespace to avoid spacing issues
        return Markup(svg_with_color.strip())
    
    def _fix_pdb_header(self, file_path):
        with open(file_path, 'r') as f:
            lines = f.readlines()
        if not any(line.startswith('HEADER') for line in lines):
            lines.insert(0, "HEADER    AUTO-GENERATED PDB FILE\n")
        with open(file_path, 'w') as f:
            f.writelines(lines)
    
    def _get_dssp_output(self, pdb_name, file_path):
        '''
        Private method that returns DSSP secondary structure assignments from a PDB or mmCIF file.
        Handles PDB files without headers gracefully.

        Returns:
            A dictionary containing chain IDs as keys and lists of residues with secondary structures as values.
            Example: {"A": [{"chain": "A", "res_num": 1, "res_name": "A", "res_struc": "H"}, {}, {}]}
        '''

        # Parser selection based on file type
        if file_path.lower().endswith('.pdb'):
            parser = PDBParser(QUIET=True)
        elif file_path.lower().endswith('.cif'):
            parser = MMCIFParser(QUIET=True)
        else:
            raise Exception('File type must be either PDB or mmCIF')

        # Load the structure
        structure = parser.get_structure(pdb_name, file_path)
        model = structure[0]

        try:
            dssp = DSSP(model, file_path)
        except Exception as e:
            self._fix_pdb_header(file_path)
            try:
                dssp = DSSP(model, file_path)
            except Exception as e:
                raise Exception(f"DSSP failed after header fix: {e}")

        # Parse DSSP output into a structured dictionary
        structure_out = {}
        for key in dssp.keys():
            chain_id = key[0][0]
            res_dict = {
                'chain': chain_id,
                'res_num': key[1][1],
                'res_name': dssp[key][1],
                'res_struc': dssp[key][2]
            }
            structure_out.setdefault(chain_id, []).append(res_dict)

        return structure_out
    
    def get_uniprot_data(self, identifier):
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{identifier}"
        response = requests.get(url)
        
        if response.status_code == 200:
            json_data = response.json() 
            return json_data
        else:
            print(f"Error: Unable to fetch data (status code: {response.status_code})")
            return None
        
    def get_rcsb_entry_data(self, identifier):
        url = f"https://data.rcsb.org/rest/v1/core/entry/{identifier}"
        response = requests.get(url)
        
        if response.status_code == 200:
            json_data = response.json() 
            return json_data
        else:
            print(f"Error: Unable to fetch data (status code: {response.status_code})")
            return None
        
    def get_scientific_name(self, identifier, rcsb_data):
        if identifier and rcsb_data:
            entity_ids = rcsb_data['rcsb_entry_container_identifiers']['polymer_entity_ids']

            for entity_id in entity_ids:
                url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{identifier}/{entity_id}"
                response = requests.get(url)
                
                if response.status_code == 200:
                    json_data = response.json() 
                else:
                    print(f"Error: Unable to fetch data (status code: {response.status_code})")
                    continue
                
                if 'rcsb_entity_source_organism' in json_data:
                    return json_data['rcsb_entity_source_organism'][0].get('scientific_name', None)
                else:
                    continue
        return None
    
    def get_uniprot_id_by_chain_id(self, data, target_chain_id):
        if data and target_chain_id:
            for structure_id, structure_data in data.items():
                for uniprot_id, uniprot_data in structure_data['UniProt'].items():
                    for mapping in uniprot_data['mappings']:
                        if mapping['chain_id'] == target_chain_id:
                            return uniprot_id
        return None

    def generate_visual(self, residues_per_line: int = 50, output_image_name: str = '', dpi: int = 100, pdf: bool = False) -> None:
        '''
        Generates a visualization with the residues mapped to their secondar structure types.

        Args:
            residues_per_line (int): An integer respresnting the number of residues per line, in the visualization, the default values is 60.
            output_image_name (string): Representing the output image path, valid types: "jpg, png, and gif".
            pdf (bool): A boolean value indicating if a PDF should be generated too. The default value is False.
            color: TO BE ADDED

        Raises:
            TypeError: If the any of the parameters are not of the correct type
            ValueError: If any of the parameters are outside the range

        '''
        if not isinstance(residues_per_line, int):
            raise TypeError('The residues_per_line parameter must be an integer')
        if not isinstance(output_image_name, str):
            raise TypeError('The ouput_image parameter must be a string')
        if not isinstance(pdf, bool):
            raise TypeError('The pdf parameter must be a boolean value (True/False)')

        output_folder = 'output'
        # Create folder if it doesn't exist
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        if output_image_name == '':
            output_image_name = f'{os.path.splitext(os.path.basename(self.file_path))[0]}.png'
        if '.' not in output_image_name or output_image_name.split('.')[-1].lower() not in ['jpg', 'png']:
            raise ValueError('The output image name must end with one of the following extensions: "JPG", "PNG"')
        
        updated_template = self._update_template(residues_per_line=residues_per_line)
        output_width, output_height = self._get_width_and_height(residues_per_line=residues_per_line)
        output_width, output_height = self._get_width_and_height(residues_per_line=residues_per_line)

        options = {'page-height': f'{output_height}px',
                   'page-width': f'{output_width}px',
                   'margin-top': '0',
                   'margin-bottom': '0',
                   'margin-left': '0',
                   'margin-right': '0',
                   'enable-local-file-access': '',
                   'print-media-type': True,
                   'disable-smart-shrinking': True}
        if dpi == 100:
            imgkit.from_string(updated_template, f'{output_folder}/{output_image_name}', css='templates/output_styles.css', options={'quality': 100})
        else:
            pdf_bytes = pdfkit.from_string(updated_template, css='templates/output_styles.css', options=options)
            pages = pdf2image.convert_from_bytes(pdf_bytes, dpi=dpi)
            if len(pages) == 1:
                pages[0].save(output_image_name)
            else:
                for count, page in enumerate(pages):
                    page.save(f'{output_folder}/{output_image_name}')

        if pdf:
            pdfkit.from_string(updated_template, f'{self.pdb_name}.pdf', css='templates/output_styles.css', options=options)
