#!/usr/bin/env bash
# render-build.sh

# Exit on error
set -o errexit

echo "Installing Python dependencies..."
pip install -r requirements.txt

echo "Setting up backend tools (DSSP & STRIDE)..."
mkdir -p ../tools/bin

# 1. Use Micromamba (fast conda alternative) to install dssp and stride locally without root
if [ ! -f "../tools/bin/micromamba" ]; then
    echo "Downloading Micromamba..."
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ../tools bin/micromamba
fi

# Install dssp=3.0.0 (which avoids the heavy libcifpp dictionaries and RAM exhaustion on Render)
echo "Installing DSSP 3.0.0 and STRIDE via bioconda..."
../tools/bin/micromamba create -y -p ../tools/conda-env -c conda-forge -c bioconda dssp=3.0.0 stride

# Create wrapper scripts to ensure the binaries can find Conda's shared libraries
echo "Creating wrapper scripts..."

cat << 'EOF' > ../tools/bin/mkdssp
#!/usr/bin/env bash
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
export LD_LIBRARY_PATH="$DIR/../conda-env/lib:$LD_LIBRARY_PATH"
# Fallback to dssp since bioconda's dssp=3.0.0 binary is named dssp
exec "$DIR/../conda-env/bin/dssp" "$@"
EOF
chmod +x ../tools/bin/mkdssp

cat << 'EOF' > ../tools/bin/dssp
#!/usr/bin/env bash
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
export LD_LIBRARY_PATH="$DIR/../conda-env/lib:$LD_LIBRARY_PATH"
exec "$DIR/../conda-env/bin/dssp" "$@"
EOF
chmod +x ../tools/bin/dssp

cat << 'EOF' > ../tools/bin/stride
#!/usr/bin/env bash
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
export LD_LIBRARY_PATH="$DIR/../conda-env/lib:$LD_LIBRARY_PATH"
exec "$DIR/../conda-env/bin/stride" "$@"
EOF
chmod +x ../tools/bin/stride

echo "Build complete! Tools installed to ../tools/bin"
