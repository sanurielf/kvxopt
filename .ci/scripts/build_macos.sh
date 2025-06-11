#!/bin/bash
set -e

echo "Installing macOS dependencies..."

# Install Homebrew packages
brew install cmake openblas glpk gsl fftw suite-sparse

# Build OSQP using bash command instead of direct execution[1][6]
bash .ci/library_builders/build_osqp.sh

echo "macOS dependencies installed successfully"