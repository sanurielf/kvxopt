#!/bin/bash
set -e

# Source global configuration (variables already loaded in GitHub Actions)
# But source again for local testing
if [ -f ".ci/config/versions.env" ]; then
    source .ci/config/versions.env
fi

platform=$1
arch=${2:-x64}

# Export ARCH for use by library builders
export ARCH="$arch"

echo "Installing dependencies for platform: $platform, architecture: $arch"
echo "OpenBLAS version: $OPENBLAS_VERSION"
echo "SuiteSparse version: $SUITESPARSE_VERSION"
echo "OSQP version: $OSQP_VERSION"

case $platform in
  linux)
    bash .ci/scripts/build_linux.sh
    ;;
  macos)
    source .ci/config/macos.env
    bash .ci/scripts/build_macos.sh
    ;;
  *)
    echo "Unknown platform: $platform"
    exit 1
    ;;
esac
