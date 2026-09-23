#!/usr/bin/env bash
# render-build.sh

# Exit on error
set -o errexit

echo "Installing Python dependencies..."
pip install -r requirements.txt

echo "Setting up backend tools (DSSP & STRIDE)..."
mkdir -p tools/bin
mkdir -p tools/share/libcifpp

# 1. Use Micromamba (fast conda alternative) to install dssp and stride locally without root
if [ ! -f "tools/bin/micromamba" ]; then
    echo "Downloading Micromamba..."
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C tools bin/micromamba
fi

echo "Installing DSSP and STRIDE via bioconda..."
./tools/bin/micromamba create -y -p ./tools/conda-env -c conda-forge -c bioconda dssp stride

# Symlink the executables so main.py can find them in tools/bin
ln -sf $(pwd)/tools/conda-env/bin/mkdssp tools/bin/mkdssp
ln -sf $(pwd)/tools/conda-env/bin/stride tools/bin/stride

# 2. DSSP/libcifpp requires compound dictionaries to parse mmCIF/PDB files.
# main.py already looks for tools/share/libcifpp, so we will download them there.
echo "Downloading components.cif..."
curl -s -o tools/share/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif

echo "Downloading mmcif_pdbx_v50.dic..."
curl -s -o tools/share/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic

echo "Downloading mmcif_ma.dic..."
curl -s -o tools/share/libcifpp/mmcif_ma.dic https://github.com/ihmwg/ModelCIF/raw/master/dist/mmcif_ma.dic

echo "Build complete! Tools installed to tools/bin and dictionaries to tools/share/libcifpp"
