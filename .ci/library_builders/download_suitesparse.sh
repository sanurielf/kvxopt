#!/bin/bash
set -e

echo "Downloading SuiteSparse ${SUITESPARSE_VERSION} for Windows..."

# Download using curl (available in Git Bash on Windows)
echo "Downloading from GitHub..."
curl -L -o "v${SUITESPARSE_VERSION}.tar.gz" \
    "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/v${SUITESPARSE_VERSION}.tar.gz"

ls


# Verify download was successful
if [ ! -f "v${SUITESPARSE_VERSION}.tar.gz" ]; then
    echo "Error: SuiteSparse download failed"
    exit 1
fi

# Get file size to ensure it's not empty
FILE_SIZE=$(stat -c%s "v${SUITESPARSE_VERSION}.tar.gz" 2>/dev/null || echo 0)
if [ "$FILE_SIZE" -lt 1000 ]; then
    echo "Error: Downloaded file is too small (${FILE_SIZE} bytes)"
    exit 1
fi

# Verify SHA256 checksum using certutil
echo "Verifying SHA256 checksum using certutil..."
ACTUAL_HASH=$(certutil -hashfile "v${SUITESPARSE_VERSION}.tar.gz" SHA256 | grep -v "SHA256" | grep -v "CertUtil" | tr -d ' \r\n' | tr '[:upper:]' '[:lower:]')
EXPECTED_HASH=$(echo "${SUITESPARSE_SHA256}" | tr '[:upper:]' '[:lower:]')

echo "Expected: ${EXPECTED_HASH}"
echo "Actual:   ${ACTUAL_HASH}"

if [ "${ACTUAL_HASH}" != "${EXPECTED_HASH}" ]; then
    echo "Error: SHA256 checksum mismatch!"
    echo "Expected: ${EXPECTED_HASH}"
    echo "Actual:   ${ACTUAL_HASH}"
    exit 1
else
    echo "SHA256 verification passed"
fi

# Extract using tar (available in Git Bash)
echo "Extracting SuiteSparse..."
tar -xf "v${SUITESPARSE_VERSION}.tar.gz"


# Move to expected location
if [ -d "SuiteSparse-${SUITESPARSE_VERSION}" ]; then
    mv "SuiteSparse-${SUITESPARSE_VERSION}" suitesparse
    echo "SuiteSparse ${SUITESPARSE_VERSION} downloaded and extracted successfully"
    pwd
    ls
    echo "Contents of the suitesparse directory:"
    ls suitesparse
else
    echo "Error: Expected directory SuiteSparse-${SUITESPARSE_VERSION} not found after extraction"
    exit 1
fi

# Clean up
rm "v${SUITESPARSE_VERSION}.tar.gz"


echo "SuiteSparse setup completed"
