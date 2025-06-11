#!/bin/bash
set -e

echo "Building OSQP ${OSQP_VERSION} from source..."

# Clone and checkout specific version
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
    OSQP_DIR="C:/osqp"
    INSTALL_PREFIX="C:/osqp-install"
else
    OSQP_DIR="/tmp/osqp"
    INSTALL_PREFIX="/usr/local"
fi

git clone --recursive https://github.com/oxfordcontrol/osqp.git "$OSQP_DIR"
cd "$OSQP_DIR"
git checkout "v${OSQP_VERSION}"
git submodule sync --recursive
git -c protocol.version=2 submodule update --init --force --depth=1 --recursive

# Build
mkdir build && cd build
cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" -DCMAKE_BUILD_TYPE=Release ..

if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
    cmake --build . --config Release --target install
else
    make -j$(nproc) && make install
    ldconfig || true
fi

echo "OSQP ${OSQP_VERSION} built and installed"
