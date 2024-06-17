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

1. **Preferred Operating System**: It is recommended to run this tool on a Linux system or Windows Subsystem for Linux (WSL) for optimal performance and compatibility.
2. **Installation Methods**: You can install all the required components manually or automatically using the provided bash script.

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
  4. **Install Python Requirements**:
     ```bash
     pip install -r requirements.txt
     ```

- **Automatic Installation**: Run the `requirements.sh` script to check and install all dependencies including Python, DSSP, wkhtmltopdf, and the required Python libraries:
  ```bash
  chmod +x requirements.sh
  ./requirements.sh
  ```

### Running ProS<sup>2</sup>Vi

You can use ProS<sup>2</sup>Vi either through a command-line parser or via a Flask-based GUI. Follow the detailed steps below for each method.

### Using the Flask-based GUI

1. **Run Main Script**: Execute `python app.py` from the terminal. This will launch a browser window.
2. **Upload PDB File**: In the browser, upload your PDB file and enter the PDB code.
3. **Optional Arguments**: Enter values in the optional input fields if needed.
4. **Submit and Visualize**: Press `Submit` to create the visualization. The visualization will be saved in the main directory.

### Using the Command-line Parser

1. **Run the Script**: Provide the required arguments such as the PDB name/code and the path to the PDB file. Example:
    ```bash
    python visual.py 1FAT pdb_folder/1fat.pdb
    ```

2. **Specify Additional Arguments**: Use positional arguments for additional options. Example:
    - **30 Residues per Line and PDF Output**:
      ```bash
      python visual.py 1FAT pdb_folder/1fat.pdb -r 30 -pdf
      ```
      Here, `-r 30` sets 30 residues per line, and `-pdf` generates a PDF output.

    - **200 DPI Resolution and PNG Output**:
      ```bash
      python visual.py 1FAT pdb_folder/1fat.pdb -o output/test.png -d 200
      ```
      Here, `-o output/test.png` sets the output path and filename to `test.png`, and `-d 200` sets the DPI to 200.
      
