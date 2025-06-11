#!/bin/bash
set -e

echo "Installing Windows dependencies for architecture: ${ARCH:-x64}"

# Source Windows-specific configuration
source .ci/config/windows.env

# Build OpenBLAS with architecture awareness
.ci/library_builders/build_openblas.sh

# Install CMake
choco install cmake --installargs 'ADD_CMAKE_TO_PATH=System'

# Build OSQP
.ci/library_builders/build_osqp.sh

echo "Windows dependencies installed successfully for ${ARCH:-x64}"
