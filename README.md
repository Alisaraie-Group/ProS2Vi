[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12554830.svg)](https://doi.org/10.5281/zenodo.12554830)
[![DOI](http://img.shields.io/badge/DOI-10.48550/arXiv.2408.03436-B31B1B.svg)](https://doi.org/10.48550/arXiv.2408.03436)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# ProS<sup>2</sup>Vi - a Python Tool for Visualizing Proteins Secondary Structure

## Quick Start with Docker (Recommended)

The easiest way to run ProS<sup>2</sup>Vi is using Docker. This works on Windows, macOS, and Linux without needing to install Python, DSSP, or other dependencies.

### Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop/)

### Run ProS<sup>2</sup>Vi
```bash
docker run -p 3000:3000 alisaraiegroup/pros2vi
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
```bash
git clone https://github.com/Alisaraie-Group/ProS2Vi.git
cd ProS2Vi
docker build -t pros2vi .
docker run -p 3000:3000 pros2vi
```

---

## Manual Setup Instructions

Follow these steps to set up ProS<sup>2</sup>Vi on your system.

### Clone the Repository

1. **Clone the GitHub Repository**: Open your terminal and clone the repository using the following command:
   ```bash
   git clone https://github.com/Alisaraie-Group/ProS2Vi.git
   ```

2. **Navigate to the Project Directory**: Change to the directory containing the cloned repository:
   ```bash
   cd ProS2Vi
   ```

### Prerequisites

1. **Preferred Operating System**: ProS<sup>2</sup>Vi needs to be run on a Linux system. For Windows users, you need to install Windows Subsystem for Linux (WSL) from the [Microsoft website](https://docs.microsoft.com/en-us/windows/wsl/install).
2. **Installation Methods**: You can install all the required components manually or automatically using the provided bash script.
2. **Installation Methods**: You can install all the required components manually or automatically using the provided shell script.

- **Manual Installation**:
  1. **Install Python**: Ensure Python is installed on your system.
     ```bash
     sudo apt-get install python3
     ```
  2. **Install DSSP** (version 4.x recommended for full PPII helix support):
     ```bash
     sudo apt-get install dssp
     ```
  3. **Install wkhtmltopdf**:
     ```bash
     sudo apt-get install wkhtmltopdf
     ```
  4. **Set Up a Virtual Environment**:
     ```bash
     python -m venv venv
     source venv/bin/activate
     ```
  5. **Install Python Requirements**:
     ```bash
     pip install -r requirements.txt
     ```

- **Automatic Installation**: Run the `install_dependencies.sh` script that check ands install all dependencies including Python, DSSP, wkhtmltopdf, and then creates a virtual environment and installs the required Python libraries:
  ```bash
  chmod +x install_dependencies.sh
  ./install_dependencies.sh
  ```

### Running ProS<sup>2</sup>Vi

You can use ProS<sup>2</sup>Vi either through a command-line parser or via a Flask-based GUI. Follow the detailed steps below for each method.

### Using the Flask-based GUI

1. **Run Main Script**: Execute `python3 pros2vi_gui.py` from the terminal. This will launch a browser window.
2. **Add the Protein Structure**: In the browser, upload your PDB file or enter the PDB code.
3. **Optional Arguments**: Enter values in the optional input fields if needed.
4. **Submit and Visualize**: Press `Submit` to create the visualization. The visualization will be saved in the main directory.

### Using the Command-line Parser

1. **Run the Script**: Provide the path to the PDF/mmCIF file and optional arguments. Example:
    ```bash
    python3 pros2vi_cli.py pdb_folder/1fat.pdb
    ```

2. **Specify Additional Arguments**: Use positional arguments for additional options. Example:
    - **30 Residues per Line and PDF Output**:
      ```bash
      python visual.py pdb_folder/1fat.pdb -r 30 -pdf
      ```
      Here, `-r 30` sets 30 residues per line, and `-pdf` generates a PDF output.

    - **200 DPI Resolution and PNG Output**:
      ```bash
      python visual.py pdb_folder/1fat.pdb -o output/test.png -d 200
      ```
      Here, `-o output/test.png` sets the output path and filename to `test.png`, and `-d 200` sets the DPI to 200.

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

If you have an older version, update DSSP to benefit from PPII helix detection:

```bash
sudo apt-get update
sudo apt-get install --upgrade dssp
```

Alternatively, install via conda-forge:

```bash
conda install -c conda-forge dssp
```

      
## Citing ProS<sup>2</sup>Vi

If you use ProS<sup>2</sup>Vi in academic work, please cite the software as

```
@software{qasim_2024_12554831,
  author       = {Qasim, Muhammad Luckman and
                  Alisaraie, Laleh},
  title        = {Alisaraie-Group/ProS2Vi: ProS2Vi 1.0.3},
  month        = jun,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v1.0.3},
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
