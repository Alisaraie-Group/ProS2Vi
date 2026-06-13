[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12554830.svg)](https://doi.org/10.5281/zenodo.12554830)
[![DOI](http://img.shields.io/badge/DOI-10.48550/arXiv.2408.03436-B31B1B.svg)](https://doi.org/10.48550/arXiv.2408.03436)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# ProS<sup>2</sup>Vi - a Python Tool for Visualizing Proteins Secondary Structure

## Quick Start with Docker (Recommended)

The easiest way to run ProS<sup>2</sup>Vi is using Docker. This works on Windows, macOS, and Linux without needing to install Python, DSSP, or other dependencies.

### Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop/)

### Run ProS<sup>2</sup>Vi

**Quick Start:**
```bash
docker run -p 3000:3000 alisaraiegroup/pros2vi
```

Then open http://localhost:3000 in your browser.

**Using Docker Compose:**
```bash
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi
docker-compose up
```

Then open http://localhost:3000 in your browser.

### Save Output Files Locally
To persist generated images to your local machine:
```bash
docker run -p 3000:3000 -v $(pwd)/output:/app/output alisaraiegroup/pros2vi
```
Generated images will be saved to the `output/` folder in your current directory.

### Build from Source (Optional)
If you prefer to build the Docker image yourself:

**Using Docker:**
```bash
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi
docker build -t pros2vi .
docker run -p 3000:3000 pros2vi
```

**Using Docker Compose:**
```bash
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi
docker-compose up --build
```

---

## Alternative Setup Methods

If you prefer not to use Docker, you can set up ProS<sup>2</sup>Vi using either Conda or Python's venv.

### Option 1: Conda Setup

Works on Linux, macOS, and Windows (via WSL).

#### Prerequisites
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download)

#### Installation

**Method 1: Using mamba (Recommended - faster and more reliable)**

```bash
# Clone the repository
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi

# Install mamba (if not already installed)
conda install -n base -c conda-forge mamba -y

# Create the environment with mamba
mamba env create -f environment.yml

# Activate the environment
conda activate pros2vi

# Install system wkhtmltopdf (required for image generation)
# Ubuntu/Debian:
sudo apt-get install -y wkhtmltopdf
# macOS:
# brew install wkhtmltopdf
```

**Method 2: Using conda directly**

```bash
# Clone the repository
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi

# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate pros2vi

# Install system wkhtmltopdf (required for image generation)
# Ubuntu/Debian:
sudo apt-get install -y wkhtmltopdf
# macOS:
# brew install wkhtmltopdf
```

**Method 3: If you encounter network issues**

Create a minimal conda environment and install packages via pip:

```bash
# Create environment with Python and system dependencies
conda create -n pros2vi python=3.11 -y
conda activate pros2vi

# Install DSSP and poppler via conda
conda install -c conda-forge dssp poppler -y

# Install system wkhtmltopdf (for proper SVG/icon rendering)
# Ubuntu/Debian:
sudo apt-get install -y wkhtmltopdf
# macOS:
# brew install wkhtmltopdf

# Install Python packages via pip
pip install -r requirements.txt
```

> **Important**: The conda-forge version of `wkhtmltopdf` (0.12.4) has limited SVG rendering support which causes icons to appear blank. Install the system version (0.12.6+) for proper icon display.

### Option 2: Python venv Setup

Works on Linux, macOS, and Windows with Python 3.11 installed.

#### Prerequisites
- Python 3.11 or higher
- System package manager (apt, brew, etc.)

#### Installation

**Ubuntu/Debian:**

```bash
# Clone the repository
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi

# Install system dependencies
sudo apt-get update
sudo apt-get install -y dssp poppler-utils wkhtmltopdf

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install Python packages
pip install -r requirements.txt
```

**macOS:**

```bash
# Clone the repository
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi

# Install system dependencies via Homebrew
brew install dssp poppler wkhtmltopdf

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install Python packages
pip install -r requirements.txt
```

**Windows (WSL):**

Follow the Ubuntu/Debian instructions above in your WSL environment.

### Running ProS<sup>2</sup>Vi

> **Note**: Before running ProS<sup>2</sup>Vi, ensure your environment is activated:
> - For Conda: `conda activate pros2vi`
> - For venv: `source venv/bin/activate` (Linux/macOS) or `venv\Scripts\activate` (Windows)

#### Using the Flask-based GUI

1. **Start the GUI**:
   ```bash
   python pros2vi_gui.py
   ```
   This will launch a browser window at http://localhost:3000.

2. **Add the Protein Structure**: Upload your PDB/mmCIF file or enter the PDB code.

3. **Configure Options**: Enter values in the optional input fields if needed.

4. **Submit and Visualize**: Press `Submit` to create the visualization. Output files are saved in the `output/` directory.

#### Using the Command-line Interface

1. **Basic Usage**: Provide the path to the PDB/mmCIF file:
   ```bash
   python pros2vi_cli.py pdb_folder/1fat.pdb
   ```

2. **With Options**:
   - **30 Residues per Line and PDF Output**:
     ```bash
     python pros2vi_cli.py pdb_folder/1fat.pdb -r 30 -pdf
     ```

   - **200 DPI Resolution and PNG Output**:
     ```bash
     python pros2vi_cli.py pdb_folder/1fat.pdb -o output/test.png -d 200
     ```

## DSSP Version Compatibility

ProS<sup>2</sup>Vi fully supports **DSSP version 4**, which includes important updates:

- **Poly-Proline II (PPII) Helix Detection**: DSSP v4 extends the secondary structure alphabet with code **'P'** to mark κ-helices (also known as Poly-Proline II helices), which are left-handed helical structures commonly found in protein loops and flexible regions.
- **mmCIF Format Support**: Primary input/output format is now mmCIF (macromolecular Crystallographic Information File) for FAIR compliance.
- **Legacy PDB Compatibility**: Maintains backward compatibility with traditional PDB file format.

### Secondary Structure Codes

ProS<sup>2</sup>Vi visualizes the following DSSP secondary structure assignments:

| Code | Structure Type | Description |
|------|----------------|-------------|
| **H** | α-helix | Alpha helix (4-12 residues per turn) |
| **G** | 3₁₀-helix | 3-10 helix (3 residues per turn) |
| **I** | π-helix | Pi helix (5 residues per turn) |
| **E** | β-strand | Extended strand in parallel/antiparallel β-sheet |
| **B** | β-bridge | Isolated β-bridge residue |
| **T** | Turn | Hydrogen bonded turn |
| **S** | Bend | Bend (high curvature) |
| **P** | PPII helix | Poly-Proline II helix / κ-helix (DSSP v4+) |
| **-** | Coil/Loop | Unstructured coil or loop |

**Note**: The 'P' code for PPII helices is a new feature in DSSP v4. These structures comprise approximately 1.9% of all residues in the PDB and are typically very short (half are ≤3 residues).

### Checking Your DSSP Version

To verify your DSSP installation supports PPII helices:

```bash
mkdssp --version  # Should show version 4.0.0 or higher
```

**Note**: The Docker image, Conda environment, and system packages (via apt/brew) all include DSSP v4+ by default, so no manual installation is needed.

      
## Citing ProS<sup>2</sup>Vi

If you use ProS<sup>2</sup>Vi in academic work, please cite the software as

```
@software{qasim_2024_12554831,
  author       = {Qasim, Muhammad Luckman and
                  Alisaraie, Laleh},
  title        = {Alisaraie-Group/ProS2Vi: ProS2Vi 1.1.0},
  month        = jun,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v1.1.0},
  doi          = {10.5281/zenodo.12554831},
  url          = {https://doi.org/10.5281/zenodo.12554831}
}
```

and the article as

```
@misc{qasim2024pros2vipythontoolvisualizing,
      title={ProS2Vi: a Python Tool for Visualizing Proteins Secondary Structure}, 
      author={Luckman Qasim and Laleh Alisaraie},
      year={2024},
      eprint={2408.03436},
      archivePrefix={arXiv},
      primaryClass={q-bio.BM},
      url={https://arxiv.org/abs/2408.03436}, 
}
```
