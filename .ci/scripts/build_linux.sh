#!/bin/bash
set -e

echo "Installing Linux dependencies..."


# Detect package manager and install packages[2][3]
if command -v yum &> /dev/null; then
    # CentOS/RHEL with yum
    if [ "$ARCH" == "aarch64" ]; then
        yum install -y epel-release
    fi
    yum install -y cmake3 make gcc-c++ openblas-devel glpk-devel gsl-devel fftw-devel git wget suitesparse-devel

elif command -v dnf &> /dev/null; then
    # Fedora/RHEL 8+ with dnf
    dnf install -y cmake make gcc-c++ openblas-devel glpk-devel gsl-devel fftw-devel git wget suitesparse-devel

elif command -v apt-get &> /dev/null; then
    # Debian/Ubuntu with apt[2]
    apt-get update
    apt-get install -y cmake build-essential libopenblas-dev libglpk-dev libgsl-dev libfftw3-dev git wget libsuitesparse-dev

elif command -v apk &> /dev/null; then
    # Alpine with apk
    apk add cmake build-base openblas-dev glpk-dev gsl-dev fftw-dev git wget suitesparse-dev
    mv /usr/include/suitesparse/suitesparse/* /usr/include/suitesparse/
else
    echo "Error: No supported package manager found (yum, dnf, apk, or apt-get)"
    exit 1
fi

# Build OSQP using bash command
bash .ci/library_builders/build_osqp.sh

echo "Linux dependencies installed successfully"
