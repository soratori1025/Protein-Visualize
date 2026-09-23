#!/usr/bin/env bash
# render-build.sh

# Exit on error
set -o errexit

echo "Installing Python dependencies..."
pip install -r requirements.txt

# DSSP/libcifpp requires compound dictionaries to parse mmCIF/PDB files.
# main.py already looks for tools/share/libcifpp, so we will download them there.
echo "Setting up libcifpp dictionaries for DSSP..."
mkdir -p tools/share/libcifpp

# Download the required dictionary files
echo "Downloading components.cif..."
curl -s -o tools/share/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif

echo "Downloading mmcif_pdbx_v50.dic..."
curl -s -o tools/share/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic

echo "Downloading mmcif_ma.dic..."
curl -s -o tools/share/libcifpp/mmcif_ma.dic https://github.com/ihmwg/ModelCIF/raw/master/dist/mmcif_ma.dic

echo "Build complete! Dictionaries installed to tools/share/libcifpp"
