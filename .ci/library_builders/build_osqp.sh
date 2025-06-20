#!/bin/bash
set -e

echo "Building OSQP ${OSQP_VERSION} from source.... OS: ${OSTYPE}, ARCH: ${ARCH:-x64} CMAKE_TARGET: ${WINDOWS_CMAKE_TARGET:-x64}"

# Clone and checkout specific version
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OSTYPE" == "win32" ]]; then
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

if [[ "$OSTYPE" == "win32" ]]; then
    echo "Generating build files for OSQP with cmake for Windows..."
    cmake -G "Visual Studio 17 2022" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -A $WINDOWS_CMAKE_TARGET -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" -DCMAKE_BUILD_TYPE=Release ..
else
    echo "Generating build files for OSQP with cmake..."
    cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" -DCMAKE_BUILD_TYPE=Release ..
fi


if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin"  || "$OSTYPE" == "win32" ]]; then
    echo "Compiling OSQP with cmake for Windows..."
    cmake --build . --config Release --target install
    cp $INSTALL_PREFIX/bin/*.dll ../../src/python/.libs
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Compiling OSQP with make for macOS..."
    make -j$(sysctl -n hw.ncpu) && sudo make install
else
    echo "Compiling OSQP with make for Linux..."
    make -j$(nproc) && make install
    ldconfig || true
fi



echo "OSQP ${OSQP_VERSION} built and installed"
