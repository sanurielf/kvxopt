
echo Initiallizing Visual Studio environment for architecture: %WINDOWS_VC_TARGET%

call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" %WINDOWS_VC_TARGET%


echo Creating .lib import files for FFTW for architecture: %WINDOWS_CMAKE_TARGET%

C:

cd C:\fftw-install\


lib /machine:%ARCH% /def:libfftw3-3.def

if errorlevel 1 (
    echo Error: lib command failed
    exit /b 1
)

echo Creating import files for FFTW completed successfully.

dir C:\fftw-install

D:
