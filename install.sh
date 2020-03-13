#!/bin/bash

echo "---------------------"
echo "----   INSTALL   ----"
echo "---------------------"
echo ""

echo "Checking python version..."
# This line is to force everyone to use a similar version of Python
VERSION=$(python3.6 --version)
if [[ ${VERSION} != *"Python 3.6"* ]]; then
    echo "ERROR! Python 3.6.x is required"
    exit 1
fi

if [[ -d "venv" ]]; then
    echo "Reinstalling environment..."
    rm -rf venv
else
    echo "Installing environment..."
fi

python3 -m pip install --user virtualenv
python3 -m venv venv

echo "Activating environment..."
source venv/bin/activate

echo "Installing packages..."
python3 -m pip install --upgrade pip
pip install -r requirements.txt
if [[ $? != 0 ]]; then
    echo "ERROR! Something went wrong"
    exit $?
fi

echo "Succesfully installed!"
echo ""
