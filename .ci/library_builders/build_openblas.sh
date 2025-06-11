#!/bin/bash
set -e

# Source configuration
source .ci/config/versions.env

if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
    echo "Building OpenBLAS ${OPENBLAS_VERSION} for Windows..."
    
    # Determine architecture and corresponding checksum[1]
    case "${ARCH:-x64}" in
        x64|AMD64|x86_64)
            OPENBLAS_ARCH="x64"
            OPENBLAS_SHA256="${OPENBLAS_SHA256_X64}"
            INSTALL_DIR="C:/openblas-x64"
            ;;
        x86|Win32|i386)
            OPENBLAS_ARCH="x86"
            OPENBLAS_SHA256="${OPENBLAS_SHA256_X86}"
            INSTALL_DIR="C:/openblas-x86"
            ;;
        *)
            echo "Error: Unknown architecture '${ARCH}'"
            exit 1
            ;;
    esac
    
    echo "Downloading OpenBLAS-${OPENBLAS_VERSION}-${OPENBLAS_ARCH}.zip..."
    
    # Download OpenBLAS binary for specific architecture[1]
    curl -L -o "openblas_${OPENBLAS_ARCH}.zip" \
        "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OPENBLAS_VERSION}/OpenBLAS-${OPENBLAS_VERSION}_${OPENBLAS_ARCH}.zip"
    
    # Verify checksum[1]
    echo "${OPENBLAS_SHA256}  openblas_${OPENBLAS_ARCH}.zip" > "openblas_${OPENBLAS_ARCH}.sha256"
    certutil -hashfile "openblas_${OPENBLAS_ARCH}.zip" SHA256 | findstr /v "hash" | findstr /v "CertUtil" > computed.sha256
    fc "openblas_${OPENBLAS_ARCH}.sha256" computed.sha256
    
    # Extract to architecture-specific directory
    powershell -Command "Expand-Archive -Path openblas_${OPENBLAS_ARCH}.zip -DestinationPath ${INSTALL_DIR}"
    
    # Set environment variables for this architecture
    echo "KVXOPT_BLAS_LIB_DIR=${INSTALL_DIR}/lib" >> $GITHUB_ENV
    echo "KVXOPT_BLAS_INC_DIR=${INSTALL_DIR}/include" >> $GITHUB_ENV  
    echo "KVXOPT_LAPACK_LIB_DIR=${INSTALL_DIR}/lib" >> $GITHUB_ENV
    echo "PATH=${INSTALL_DIR}/bin;$PATH" >> $GITHUB_ENV
    
    echo "OpenBLAS ${OPENBLAS_VERSION} (${OPENBLAS_ARCH}) installed to ${INSTALL_DIR}"
    
else
    echo "OpenBLAS installation handled by system package manager on non-Windows platforms"
fi
