#!/bin/bash
set -e

# Source configuration
source .ci/config/versions.env

if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OSTYPE" == "win32" ]]; then
    echo "Building OpenBLAS ${OPENBLAS_VERSION} for Windows..."

    # Determine architecture and corresponding checksum
    case "${ARCH:-x64}" in
        AMD64|x86_64)
            OPENBLAS_ARCH="x64"
            OPENBLAS_SHA256="${OPENBLAS_SHA256_X64}"
            INSTALL_DIR="C:/openblas-install"
            ;;
        x86|Win32|i386)
            OPENBLAS_ARCH="x86"
            OPENBLAS_SHA256="${OPENBLAS_SHA256_X86}"
            INSTALL_DIR="C:/openblas-install"
            ;;
        *)
            echo "Error: Unknown architecture '${ARCH}'"
            exit 1
            ;;
    esac

    echo "Downloading OpenBLAS-${OPENBLAS_VERSION}-${OPENBLAS_ARCH}.zip..."

    # Download OpenBLAS binary
    # https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.30/OpenBLAS-0.3.30-x64.zip
    curl -L -o "openblas-${OPENBLAS_ARCH}.zip" \
        "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OPENBLAS_VERSION}/OpenBLAS-${OPENBLAS_VERSION}-${OPENBLAS_ARCH}.zip"

    # Verify download was successful
    if [ ! -f "openblas-${OPENBLAS_ARCH}.zip" ] || [ $(stat -c%s "openblas-${OPENBLAS_ARCH}.zip" 2>/dev/null || echo 0) -lt 1000 ]; then
        echo "Error: OpenBLAS download failed or file is too small"
        exit 1
    fi

    echo "Verifying SHA256 checksum using certutil..."
    ACTUAL_HASH=$(certutil -hashfile "openblas-${OPENBLAS_ARCH}.zip" SHA256 | grep -v "SHA256" | grep -v "CertUtil" | tr -d ' \r\n' | tr '[:upper:]' '[:lower:]')
    EXPECTED_HASH=$(echo "${OPENBLAS_SHA256}" | tr '[:upper:]' '[:lower:]')

    echo "Expected: ${EXPECTED_HASH}"
    echo "Actual:   ${ACTUAL_HASH}"

    if [ "${ACTUAL_HASH}" != "${EXPECTED_HASH}" ]; then
        echo "Error: SHA256 mismatch!"
        echo "Expected: ${EXPECTED_HASH}"
        echo "Actual:   ${ACTUAL_HASH}"
        exit 1
    else
        echo "SHA256 verification passed"
    fi

    # Extract using PowerShell
    echo "Extracting OpenBLAS to ${INSTALL_DIR}..."
    powershell -Command "Expand-Archive -Path 'openblas-${OPENBLAS_ARCH}.zip' -DestinationPath '${INSTALL_DIR}' -Force"

    dir "${INSTALL_DIR}"
    dir "${INSTALL_DIR}/bin"
    cp "${INSTALL_DIR}/bin/libopenblas.dll" "src/python/.libs/libopenblas.dll"
    dir "${INSTALL_DIR}/lib"
    dir "${INSTALL_DIR}/include"
    if [ ! -d "${INSTALL_DIR}/bin" ] || [ ! -d "${INSTALL_DIR}/lib" ] || [ ! -d "${INSTALL_DIR}/include" ]; then
        echo "Error: Expected directories not found after extraction"
        exit 1
    fi

    # Set environment variables for GitHub Actions
    echo "KVXOPT_BLAS_LIB_DIR=${INSTALL_DIR}/lib" >> $GITHUB_ENV
    echo "KVXOPT_BLAS_INC_DIR=${INSTALL_DIR}/include" >> $GITHUB_ENV
    echo "KVXOPT_LAPACK_LIB_DIR=${INSTALL_DIR}/lib" >> $GITHUB_ENV
    echo "${INSTALL_DIR}/bin" >> $GITHUB_PATH

    echo "OpenBLAS ${OPENBLAS_VERSION} (${OPENBLAS_ARCH}) installed to ${INSTALL_DIR}"

else
    echo "OpenBLAS installation handled by system package manager on non-Windows platforms"
fi
