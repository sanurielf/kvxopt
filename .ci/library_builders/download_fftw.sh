#!/bin/bash
set -e

echo "Downloading FFFTW ${FFTW_VERSION} for Windows..."

# Download using curl (available in Git Bash on Windows)
curl -L -o "fftw-${FFTW_VERSION}-dll${WINDOWS_FFTW_TARGET}.zip" \
    "ftp://ftp.fftw.org/pub/fftw/fftw-${FFTW_VERSION}-dll${WINDOWS_FFTW_TARGET}.zip"

ls


# Verify download was successful
if [ ! -f "fftw-${FFTW_VERSION}-dll${WINDOWS_FFTW_TARGET}.zip" ]; then
    echo "Error: FFTW download failed"
    exit 1
fi

# Get file size to ensure it's not empty
FILE_SIZE=$(stat -c%s "fftw-${FFTW_VERSION}-dll${WINDOWS_FFTW_TARGET}.zip" 2>/dev/null || echo 0)
if [ "$FILE_SIZE" -lt 1000 ]; then
    echo "Error: Downloaded file is too small (${FILE_SIZE} bytes)"
    exit 1
fi

# Verify SHA256 checksum using certutil
echo "Verifying SHA256 checksum using certutil..."
ACTUAL_HASH=$(certutil -hashfile "fftw-${FFTW_VERSION}-dll${WINDOWS_FFTW_TARGET}.zip" SHA256 | grep -v "SHA256" | grep -v "CertUtil" | tr -d ' \r\n' | tr '[:upper:]' '[:lower:]')

if [ ${WINDOWS_FFTW_TARGET} == "64" ]; then
    FFTW_SHA256=${FFTW_SHA256_X64}
else
    FFTW_SHA256=${FFTW_SHA256_X86}
fi

EXPECTED_HASH=$(echo "${FFTW_SHA256}" | tr '[:upper:]' '[:lower:]')

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
echo "Extracting FFTW..."
mkdir C:/fftw-install

unzip "fftw-${FFTW_VERSION}-dll${WINDOWS_FFTW_TARGET}.zip" -d C:/fftw-install

dir C:/fftw-install

# Clean up
rm "fftw-${FFTW_VERSION}-dll${WINDOWS_FFTW_TARGET}.zip"

cp C:/fftw-install/*.dll ./src/python/.libs


echo "FFTW setup completed"
