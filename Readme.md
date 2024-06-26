[![DOI](https://zenodo.org/badge/816348234.svg)](https://zenodo.org/doi/10.5281/zenodo.12554830)

# ProS<sup>2</sup>Vi - a Python Tool for Visualizing Proteins Secondary Structure

## Setup Instructions

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
  2. **Install DSSP**:
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
      
