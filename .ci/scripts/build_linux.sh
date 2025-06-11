#!/bin/bash
set -e

echo "Installing Linux dependencies..."

# Check if we're in a container (no sudo needed) or on host system
if [ -f /.dockerenv ] || [ -n "$container" ]; then
    # Inside container - run without sudo
    if command -v yum &> /dev/null; then
        yum install -y cmake3 make gcc-c++ openblas-devel glpk-devel gsl-devel fftw-devel git wget
    elif command -v apt-get &> /dev/null; then
        apt-get update
        apt-get install -y cmake build-essential libopenblas-dev libglpk-dev libgsl-dev libfftw3-dev git wget
    fi
else
    # On host system - use sudo
    if command -v yum &> /dev/null; then
        sudo yum install -y cmake3 make gcc-c++ openblas-devel glpk-devel gsl-devel fftw-devel git wget
    elif command -v apt-get &> /dev/null; then
        sudo apt-get update
        sudo apt-get install -y cmake build-essential libopenblas-dev libglpk-dev libgsl-dev libfftw3-dev git wget
    fi
fi


# Build OSQP using bash command
bash .ci/library_builders/build_osqp.sh

echo "Linux dependencies installed successfully"
