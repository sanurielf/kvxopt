
echo Initiallizing Visual Studio environment for architecture: %WINDOWS_VC_TARGET%

call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" %WINDOWS_VC_TARGET%


echo Building GSL for architecture: %WINDOWS_CMAKE_TARGET%

git clone --recursive https://github.com/ampl/gsl.git

pushd gsl

git checkout %GSL_COMMIT_HASH%

md build

pushd build

cmake -G "Visual Studio 17 2022" -A %WINDOWS_CMAKE_TARGET% -DGSL_INSTALL_MULTI_CONFIG=ON -DBUILD_SHARED_LIBS=ON -DMSVC_RUNTIME_DYNAMIC=ON -DCMAKE_INSTALL_PREFIX="C:/gsl-install" ..

cmake --build . --target install


popd
popd

cp C:/gsl-install/bin/Debug/*.dll src/python/.libs
