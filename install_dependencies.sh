#!/bin/bash

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

if ! command_exists dssp; then
    echo "dssp not found. Installing dssp..."
    sudo apt-get update
    sudo apt-get install dssp -y
else
    echo "dssp is already installed."
fi

if ! command_exists wkhtmltopdf; then
    echo "wkhtmltopdf not found. Installing wkhtmltopdf..."
    sudo apt-get update
    sudo apt-get install wkhtmltopdf -y
else
    echo "wkhtmltopdf is already installed."
fi

if ! command_exists python3; then
    echo "Python not found. Installing Python..."
    sudo apt-get update
    sudo apt-get install python3 -y
else
    echo "Python is already installed."
fi

if ! command_exists pip3; then
    echo "pip not found. Installing pip..."
    sudo apt-get update
    sudo apt-get install python3-pip -y
else
    echo "pip is already installed."
fi

if [ -f "requirements.txt" ]; then
    echo "Installing Python libraries from requirements.txt..."
    pip3 install -r requirements.txt
else
    echo "requirements.txt not found."
fi

echo "All installations and checks are complete."
