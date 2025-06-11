@echo off
REM Set GLPK version (default to 5.0 if not set)
set "GLPK_VERSION=%1"
if "%GLPK_VERSION%"=="" set "GLPK_VERSION=5.0"

REM Set archive name and download URL
set "GLPK_TAR_GZ=glpk-%GLPK_VERSION%.tar.gz"
set "GLPK_TAR=glpk-%GLPK_VERSION%.tar"
set "GLPK_URL=http://ftp.gnu.org/gnu/glpk/glpk-%GLPK_VERSION%.tar.gz"

REM Download GLPK
echo Downloading GLPK version %GLPK_VERSION%...
curl -v -L -o "%GLPK_TAR_GZ%" "%GLPK_URL%"

echo Curl exited with errorlevel %errorlevel%

if %errorlevel% neq 0 (
    echo Error: GLPK download failed
    exit /b 1
)

REM Check if file exists and has reasonable size (>100KB)
for %%I in ("%GLPK_TAR_GZ%") do set "FILESIZE=%%~zI"
if %FILESIZE% LSS 100000 (
    echo Error: Downloaded file is too small, likely failed.
    del "%GLPK_TAR_GZ%"
    exit /b 1
)

REM Extract .tar from .tar.gz using 7-Zip
echo Extracting %GLPK_TAR_GZ% to %GLPK_TAR%...
7z x "%GLPK_TAR_GZ%" -aoa
if errorlevel 1 (
    echo Error: Failed to extract .tar from .tar.gz
    exit /b 1
)

REM Extract contents from .tar using 7-Zip
echo Extracting %GLPK_TAR% to glpk-%GLPK_VERSION%...
mkdir "glpk-%GLPK_VERSION%"
7z x "%GLPK_TAR%" -o. -aoa
if errorlevel 1 (
    echo Error: Failed to extract files from .tar
    exit /b 1
)

REM Rename to only glpk
move /Y "glpk-%GLPK_VERSION%" "glpk"

REM Cleanup
del "%GLPK_TAR_GZ%"
del "%GLPK_TAR%"

echo GLPK %GLPK_VERSION% downloaded and extracted successfully.

dir "glpk"
exit /b 0
