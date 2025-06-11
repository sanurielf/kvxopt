@echo off
setlocal

echo Initiallizing Visual Studio environment for architecture: %WINDOWS_VC_TARGET%

call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" %WINDOWS_VC_TARGET%


echo Building GLPK for architecture: %WINDOWS_VC_TARGET% and GLP target %WINDOWS_GLPK_TARGET%


pushd glpk\%WINDOWS_GLPK_TARGET%


copy config_VC config.h


nmake /f Makefile_VC glpk.lib
if errorlevel 1 (
    echo Error: nmake failed
    popd
    exit /b 1
)


if not exist C:\glpk-install\lib md C:\glpk-install\lib
if not exist C:\glpk-install\include md C:\glpk-install\include

copy glpk.lib C:\glpk-install\lib\glpk.lib
copy ..\src\glpk.h C:\glpk-install\include\


dir C:\glpk-install\include
dir C:\glpk-install\lib

popd

cp C:\glpk-install\bin\*.dll src/python/.libs
cp C:\glpk-install\lib\*.dll src/python/.libs

endlocal