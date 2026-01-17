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
import logging
import tempfile
import shutil
import re
import uuid
from typing import Optional, Dict, Any, List, Tuple
from src.assets import ICONS

# Rendering dimension constants
CELL_WIDTH = 20  # Width of each residue cell in pixels
ROW_HEIGHT = 96  # Height of each data row in pixels
CHAIN_HEADER_HEIGHT = 76  # Height of chain header section in pixels
COUNT_COLUMN_WIDTH = 48  # Width of the residue count column in pixels
UNIPROT_COLUMN_WIDTH = 100  # Width of UniProt ID column in pixels
PADDING_WIDTH = 50  # Additional padding width in pixels
TITLE_HEIGHT = 100  # Height of title section in pixels
LEGEND_HEIGHT = 77  # Height of legend section in pixels
FOOTER_HEIGHT = 79  # Height of footer section in pixels
EXTRA_PADDING_HEIGHT = 400  # Extra vertical padding in pixels
SCALE_FACTOR = 2  # Division factor for final output dimensions

# API timeout for network requests (seconds)
API_TIMEOUT = 10


def _check_poppler_available() -> bool:
    """
    Check if poppler-utils is installed (required for pdf2image).

    Returns:
        True if poppler is available, False otherwise.
    """
    return shutil.which('pdftoppm') is not None or shutil.which('pdftocairo') is not None


class VisualMap:
    '''
    Computes a Visualization mapping, residues with their secondary structure types.

    '''
    # Package-relative paths for templates
    _TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), '..', 'templates')
    _CSS_FILE = os.path.join(os.path.dirname(__file__), '..', 'templates', 'output_styles.css')

    # Cached Jinja2 environment (class-level for performance)
    _jinja_env = None

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

    @classmethod
    def _get_jinja_env(cls) -> Environment:
        """Get cached Jinja2 environment, creating it if necessary."""
        if cls._jinja_env is None:
            cls._jinja_env = Environment(loader=FileSystemLoader(cls._TEMPLATE_DIR))
        return cls._jinja_env
    
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

    def _prepare_chain_data(self, residues_per_line: int = 50) -> List[Dict[str, Any]]:
        '''
        Prepare structured data for template rendering.

        Returns:
            List of chain data dictionaries, each containing:
            - chain_id: Chain identifier
            - uniprot_id: UniProt ID if available
            - rows: List of row data, each with annotation_cells, structure_cells, residue_cells, start_res_num, end_res_num
        '''
        unitprot_data = self.get_uniprot_data(self.pdb_name)
        chains_data = []

        for chain_id, chain in self.structure_list.items():
            rows = []
            
            for row_start in range(0, len(chain), residues_per_line):
                row_end = min(row_start + residues_per_line, len(chain))
                row_residues = chain[row_start:row_end]
                
                # Build cells for this row
                annotation_cells = []
                structure_cells = []
                residue_cells = []
                
                curr_structure_start = 0
                index_count = {'H': 0, 'B': 0, 'E': 0, 'I': 0, 'G': 0, 'S': 0, 'T': 0}
                
                for i, res in enumerate(row_residues):
                    global_index = row_start + i
                    
                    # Get the structure icon
                    is_last_in_chain = global_index == len(chain) - 1
                    is_structure_change = is_last_in_chain or chain[global_index]['res_struc'] != chain[global_index + 1]['res_struc']
                    
                    if is_structure_change and chain[global_index]['res_struc'] in ['B', 'E']:
                        svg_string = ICONS[f"{res['res_struc']}_A"]
                        color = VisualMap.COLORS[f"{res['res_struc']}_A_COLOR"]
                    else:
                        svg_string = ICONS[res['res_struc']]
                        color = VisualMap.COLORS[f"{res['res_struc']}_COLOR"]
                    
                    structure_icon = self._inject_svg_color(svg_string, color)
                    
                    # Handle annotation cells (with colspan)
                    is_row_end = i == len(row_residues) - 1
                    if is_structure_change or is_row_end:
                        col_span = i + 1 - curr_structure_start
                        curr_structure_start = i + 1
                        
                        annotation_text = ''
                        if chain[global_index]["res_struc"] in index_count.keys() and is_structure_change and not is_last_in_chain:
                            if chain[global_index]['res_struc'] != chain[global_index + 1]['res_struc']:
                                index_count[chain[global_index]['res_struc']] += 1
                                annotation_text = f"{chain[global_index]['res_struc']}{index_count[chain[global_index]['res_struc']]}"
                        
                        annotation_cells.append({'colspan': col_span, 'text': annotation_text})
                    
                    # Add structure and residue cells
                    structure_cells.append({'icon': structure_icon})
                    residue_cells.append({'name': res["res_name"]})
                
                # Pad to residues_per_line
                padding_needed = residues_per_line - len(row_residues)
                for _ in range(padding_needed):
                    structure_cells.append({'icon': ''})
                    residue_cells.append({'name': ''})
                    annotation_cells.append({'colspan': 1, 'text': ''})
                
                # Get start and end residue numbers
                start_res_num = chain[row_start]['res_num']
                end_res_num = chain[row_end - 1]['res_num']
                
                rows.append({
                    'annotation_cells': annotation_cells,
                    'structure_cells': structure_cells,
                    'residue_cells': residue_cells,
                    'start_res_num': start_res_num,
                    'end_res_num': end_res_num
                })
            
            uniprot_id = self.get_uniprot_id_by_chain_id(unitprot_data, chain_id)
            chains_data.append({
                'chain_id': chain_id,
                'uniprot_id': uniprot_id,
                'rows': rows
            })
        
        return chains_data

    def _update_template(self, residues_per_line: int = 50) -> str:
        '''
        A private Method that updates the template and returns a string representing the updated template.

        '''
        template = self._get_jinja_env().get_template("template.html.jinja")
        
        # Prepare structured data for chains
        chains_data = self._prepare_chain_data(residues_per_line)
        
        # Get metadata
        rcsb_data = self.get_rcsb_entry_data(self.pdb_name)
        if self.pdb_name:
            self.pdb_name = self.pdb_name.upper()
        if self.subtitle is None and rcsb_data:
            struct_data = rcsb_data.get('struct')
            if struct_data and 'title' in struct_data:
                self.subtitle = struct_data['title']
        if self.scientific_name is None:
            self.scientific_name = self.get_scientific_name(self.pdb_name, rcsb_data)

        content = template.render(
            pdb_name = self.pdb_name,
            pdb_title = self.subtitle,
            scientific_name = self.scientific_name,
            chains_data = chains_data,
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
    
    def _get_width_and_height(self, residues_per_line: int) -> Tuple[int, int]:
        '''
        Private method that returns the width and height of the resulting visualization in the form of a table. Eg. (output_width, output_height)

        '''
        output_width = (CELL_WIDTH * residues_per_line + COUNT_COLUMN_WIDTH + UNIPROT_COLUMN_WIDTH + PADDING_WIDTH) / SCALE_FACTOR
        num_rows = sum(math.ceil(len(chain) / residues_per_line) for chain in self.structure_list.values())
        output_height = (
            len(self.structure_list) * CHAIN_HEADER_HEIGHT +
            num_rows * ROW_HEIGHT +
            TITLE_HEIGHT + LEGEND_HEIGHT + FOOTER_HEIGHT + EXTRA_PADDING_HEIGHT
        ) / SCALE_FACTOR
        return (int(output_width), int(output_height))
    
    @staticmethod
    def _inject_svg_color(svg_string, color):
        '''
        Helper function that injects a color into an SVG string and wraps it for inline display.
        Also makes gradient IDs unique to prevent conflicts when multiple SVGs are on the same page.

        Args:
            svg_string (str): The SVG template string with {color} placeholder
            color (str): The color code to inject (e.g., '#0000ff')

        Returns:
            Markup: The SVG string with color injected, marked as safe HTML for Jinja2 rendering
        '''
        # Inject the color into the SVG
        svg_with_color = svg_string.format(color=color)

        # Generate a unique suffix for gradient IDs to prevent conflicts
        # when multiple SVGs of the same type are rendered on one page
        unique_suffix = uuid.uuid4().hex[:8]

        # Find all gradient IDs and make them unique
        # Pattern matches id="grad-..." in definitions
        def replace_gradient_id(match):
            original_id = match.group(1)
            return f'id="{original_id}-{unique_suffix}"'

        # Pattern matches url(#grad-...) in fill references
        def replace_gradient_ref(match):
            original_id = match.group(1)
            return f'url(#{original_id}-{unique_suffix})'

        # Pattern matches xlink:href="#grad-..." in gradient references
        def replace_xlink_ref(match):
            original_id = match.group(1)
            return f'xlink:href="#{original_id}-{unique_suffix}"'

        # Apply replacements for gradient definitions, fill references, and xlink references
        svg_with_color = re.sub(r'id="(grad-[^"]+)"', replace_gradient_id, svg_with_color)
        svg_with_color = re.sub(r'url\(#(grad-[^)]+)\)', replace_gradient_ref, svg_with_color)
        svg_with_color = re.sub(r'xlink:href="#(grad-[^"]+)"', replace_xlink_ref, svg_with_color)

        # Return as Markup to prevent HTML escaping in Jinja2 template
        # Strip whitespace to avoid spacing issues
        return Markup(svg_with_color.strip())
    
    def _fix_pdb_header(self, file_path: str) -> str:
        """
        Create a temporary copy of the PDB file with a HEADER line if missing.

        Returns:
            Path to the temporary file with the fixed header.
        """
        with open(file_path, 'r') as f:
            lines = f.readlines()

        if any(line.startswith('HEADER') for line in lines):
            # No fix needed, return original path
            return file_path

        # Create temp file with header added
        suffix = os.path.splitext(file_path)[1]
        temp_fd, temp_path = tempfile.mkstemp(suffix=suffix)
        try:
            with os.fdopen(temp_fd, 'w') as f:
                f.write("HEADER    AUTO-GENERATED PDB FILE\n")
                f.writelines(lines)
        except Exception:
            os.unlink(temp_path)
            raise
        return temp_path
    
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

        temp_file_path = None
        try:
            dssp = DSSP(model, file_path)
        except Exception:
            # Try with a header-fixed temporary copy
            temp_file_path = self._fix_pdb_header(file_path)
            if temp_file_path != file_path:
                # Need to reload structure from temp file
                structure = parser.get_structure(pdb_name, temp_file_path)
                model = structure[0]
            try:
                dssp = DSSP(model, temp_file_path)
            except Exception as e:
                raise Exception(f"DSSP failed after header fix: {e}")

        # Parse DSSP output into a structured dictionary
        structure_out = {}
        try:
            for key in dssp.keys():
                chain_id = key[0][0]
                res_dict = {
                    'chain': chain_id,
                    'res_num': key[1][1],
                    'res_name': dssp[key][1],
                    'res_struc': dssp[key][2]
                }
                structure_out.setdefault(chain_id, []).append(res_dict)
        finally:
            # Clean up temporary file if created
            if temp_file_path and temp_file_path != file_path:
                try:
                    os.unlink(temp_file_path)
                except OSError:
                    pass

        return structure_out
    
    def get_uniprot_data(self, identifier: str) -> Optional[Dict[str, Any]]:
        """Fetch UniProt mapping data from PDBe API."""
        if not identifier:
            return None
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{identifier}"
        try:
            response = requests.get(url, timeout=API_TIMEOUT)
            if response.status_code == 200:
                return response.json()
            else:
                logging.warning(f"Unable to fetch UniProt data (status code: {response.status_code})")
                return None
        except requests.RequestException as e:
            logging.warning(f"Request failed for UniProt data: {e}")
            return None
        
    def get_rcsb_entry_data(self, identifier: str) -> Optional[Dict[str, Any]]:
        """Fetch entry data from RCSB PDB API."""
        if not identifier:
            return None
        url = f"https://data.rcsb.org/rest/v1/core/entry/{identifier}"
        try:
            response = requests.get(url, timeout=API_TIMEOUT)
            if response.status_code == 200:
                return response.json()
            else:
                logging.warning(f"Unable to fetch RCSB entry data (status code: {response.status_code})")
                return None
        except requests.RequestException as e:
            logging.warning(f"Request failed for RCSB entry data: {e}")
            return None
        
    def get_scientific_name(self, identifier: str, rcsb_data: Optional[Dict[str, Any]]) -> Optional[str]:
        """Fetch scientific name from RCSB polymer entity API."""
        if not identifier or not rcsb_data:
            return None

        # Safely access nested entity IDs
        container_ids = rcsb_data.get('rcsb_entry_container_identifiers')
        if not container_ids:
            return None
        entity_ids = container_ids.get('polymer_entity_ids', [])

        for entity_id in entity_ids:
            url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{identifier}/{entity_id}"
            try:
                response = requests.get(url, timeout=API_TIMEOUT)
                if response.status_code != 200:
                    logging.warning(f"Unable to fetch polymer entity data (status code: {response.status_code})")
                    continue
                json_data = response.json()
            except requests.RequestException as e:
                logging.warning(f"Request failed for polymer entity data: {e}")
                continue

            source_organisms = json_data.get('rcsb_entity_source_organism')
            if source_organisms and len(source_organisms) > 0:
                return source_organisms[0].get('scientific_name')

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
            imgkit.from_string(updated_template, f'{output_folder}/{output_image_name}', css=self._CSS_FILE, options={'quality': 100})
        else:
            # High-DPI rendering requires poppler for pdf2image
            if not _check_poppler_available():
                raise RuntimeError(
                    "High-DPI rendering (dpi != 100) requires poppler-utils to be installed. "
                    "Install it via: apt-get install poppler-utils (Debian/Ubuntu), "
                    "brew install poppler (macOS), or choco install poppler (Windows)."
                )
            try:
                pdf_bytes = pdfkit.from_string(updated_template, css=self._CSS_FILE, options=options)
                pages = pdf2image.convert_from_bytes(pdf_bytes, dpi=dpi)
            except Exception as e:
                raise RuntimeError(f"Failed to render high-DPI image: {e}")
            if len(pages) == 1:
                pages[0].save(f'{output_folder}/{output_image_name}')
            else:
                for count, page in enumerate(pages):
                    page.save(f'{output_folder}/{output_image_name}')

        if pdf:
            pdfkit.from_string(updated_template, f'{self.pdb_name}.pdf', css=self._CSS_FILE, options=options)
