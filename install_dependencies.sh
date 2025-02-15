#!/bin/bash

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Install DSSP
if ! command_exists dssp; then
    echo "dssp not found. Installing dssp..."
    sudo apt-get update
    sudo apt-get install dssp -y
else
    echo "dssp is already installed."
fi

# Install wkhtmltopdf
if ! command_exists wkhtmltopdf; then
    echo "wkhtmltopdf not found. Installing wkhtmltopdf..."
    sudo apt-get update
    sudo apt-get install wkhtmltopdf -y
else
    echo "wkhtmltopdf is already installed."
fi

# Install Python
if ! command_exists python3; then
    echo "Python not found. Installing Python..."
    sudo apt-get update
    sudo apt-get install python3 -y
else
    echo "Python is already installed."
fi

# Install pip
if ! command_exists pip3; then
    echo "pip not found. Installing pip..."
    sudo apt-get update
    sudo apt-get install python3-pip -y
else
    echo "pip is already installed."
fi

# Install python3-venv
if ! dpkg -s python3-venv >/dev/null 2>&1; then
    echo "python3-venv not found. Installing python3-venv..."
    sudo apt-get update
    sudo apt-get install python3-venv -y
else
    echo "python3-venv is already installed."
fi

# Create virtual environment
VENV_DIR="venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv "$VENV_DIR" || { echo "Failed to create virtual environment"; exit 1; }
else
    echo "Python virtual environment already exists."
fi

# Activate virtual environment
echo "Activating the Python virtual environment..."
source "$VENV_DIR/bin/activate" || { echo "Failed to activate virtual environment"; exit 1; }

# Upgrade pip inside the virtual environment
pip install --upgrade pip

# Install dependencies if requirements.txt exists
if [ -f "requirements.txt" ]; then
    echo "Installing Python libraries from requirements.txt..."
    pip install --no-cache-dir -r requirements.txt || { echo "Failed to install requirements"; exit 1; }
else
    echo "requirements.txt not found."
fi

echo "All installations and checks are complete."
