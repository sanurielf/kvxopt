@echo off
setlocal enabledelayedexpansion

REM Load environment variables from versions.env
if exist ..\config\versions.env (
    for /f "tokens=*" %%i in (..\config\versions.env) do set %%i
)

REM Get platform and architecture arguments
set platform=%1
set arch=%2

if "%arch%"=="" (
    set arch=x64
)

REM Export ARCH for use by library builders
set ARCH=%arch%

echo Installing dependencies for platform: %platform%, architecture: %arch%
echo OpenBLAS version: %OPENBLAS_VERSION%
echo SuiteSparse version: %SUITESPARSE_VERSION%
echo OSQP version: %OSQP_VERSION%

@echo on
call .ci\scripts\build_windows.bat

@echo off

endlocal
