from src import visual
from flask import Flask, request, render_template, redirect, url_for, send_from_directory
import webbrowser
from threading import Timer
import os
import json
from PIL.Image import DecompressionBombError
from Bio.PDB import PDBList

app = Flask(__name__)

UPLOAD_FOLDER = os.path.join(os.getcwd(), 'uploads')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    visual.VisualMap.COLORS = {
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
        'P_COLOR': '#b7b7b7'
    }
    return render_template('index.html')

@app.route('/submit', methods=['POST'])
def submit():
    pdb_name = request.form['pdb_name'] if request.form['pdb_name'] else None
    file = request.files['file']
    title = request.form['title'] if request.form['title'] else None
    subtitle = request.form['subtitle'] if request.form['subtitle'] else None
    scientific_name = request.form['scientific_name'] if request.form['scientific_name'] else None
    residues_per_line = int(request.form['residues_per_line']) if request.form['residues_per_line'] else 50
    dpi = int(request.form['dpi']) if request.form['dpi'] else 100
    global output_image
    output_image = request.form['output_image']
    generate_pdf = 'checkbox' in request.form
    
    if request.headers.get('Content-Type') == 'application/json':
        request_data = request.get_json()
    else:
        request_data = request.form.to_dict()
    elements_colors_json = request_data.get('element_colors')
    colors_parsed = json.loads(elements_colors_json)

    filename = file.filename

    if filename == '':
        if not os.path.isfile(f'uploads/{pdb_name}.cif'):
            try:
                pdbl = PDBList()
                pdbl.retrieve_pdb_file(pdb_name, pdir='uploads', file_format='mmCif')
            except Exception as e:
                return redirect(url_for('error_page'))

        file_path = f'uploads/{pdb_name.lower()}.cif'
    else:
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(file_path)

    if output_image == '':
        output_image = f'{pdb_name}.jpg'

    pdb_name = title if title else pdb_name

    vs = visual.VisualMap(pdb_name=pdb_name, file_path=file_path, subtitle=subtitle, scientific_name=scientific_name)

    if 'helix' in colors_parsed.keys():
        visual.VisualMap.COLORS['H_COLOR'] = colors_parsed['helix']
    if 'beta' in colors_parsed.keys():
        visual.VisualMap.COLORS['B_COLOR'] = colors_parsed['beta']
        visual.VisualMap.COLORS['B_A_COLOR'] = colors_parsed['beta']
    if 'strand' in colors_parsed.keys():
        visual.VisualMap.COLORS['E_COLOR'] = colors_parsed['strand']
        visual.VisualMap.COLORS['E_A_COLOR'] = colors_parsed['strand']
    if 'pihelix' in colors_parsed.keys(): 
        visual.VisualMap.COLORS['I_COLOR'] = colors_parsed['pihelix']
    if '310helix' in colors_parsed.keys():
        visual.VisualMap.COLORS['G_COLOR'] = colors_parsed['310helix']
    if 'turn' in colors_parsed.keys():
        visual.VisualMap.COLORS['T_COLOR'] = colors_parsed['turn']
    if 'bend' in colors_parsed.keys():
        visual.VisualMap.COLORS['S_COLOR'] = colors_parsed['bend']
    if 'unsolved' in colors_parsed.keys():
        visual.VisualMap.COLORS['-_COLOR'] = colors_parsed['unsolved']
        visual.VisualMap.COLORS['P_COLOR'] = colors_parsed['unsolved']

    try:
        vs.generate_visual(residues_per_line=residues_per_line, output_image_name=output_image, dpi=dpi, pdf=generate_pdf)
    except DecompressionBombError:
        return redirect(url_for('error_page_size'))
    except Exception as e:
        return redirect(url_for('error_page'))

    return redirect(url_for('result'))

@app.route('/result')
def result():
    return render_template("display_image.html", filename=output_image)

@app.route('/result_image/<filename>')
def result_image(filename):
    return send_from_directory(os.getcwd(), filename)

@app.route('/error_page_size')
def error_page_size():
    return render_template('error_template_size.html')

@app.route('/error_page')
def error_page():
    return render_template('error_template.html')

def open_browser():
      webbrowser.open_new("http://127.0.0.1:3000")

if __name__ == "__main__":
      Timer(1, open_browser).start()
      app.run(port=3000)

