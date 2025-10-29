#!/bin/bash

echo "y" | conda create -n tidyscreen python=3.12

echo "Installing TidyScreen"

conda run -n tidyscreen pip install git+https://github.com/alfredoq/TidyScreen_v2@develop

echo "Installing Dependencies"

echo "y" | conda install -n tidyscreen -c conda-forge ambertools==23.6 espaloma espaloma_charge chemicalite visidata vmd-python

conda run -n tidyscreen pip install git+https://github.com/forlilab/Meeko@develop

echo "y" | conda install -n tidyscreen -c bioconda autodock autogrid

echo "y" | conda install -n tidyscreen -c conda-forge redis-server

echo "Installing Ersilia Hub"

# Try to delete if already exists
if [ -d "/tmp/ersilia" ]; then
    rm -rf /tmp/ersilia
    echo "folder deleted"
fi

cd /tmp && git clone https://github.com/ersilia-os/ersilia.git

cd ersilia && conda run -n tidyscreen pip install -e .


