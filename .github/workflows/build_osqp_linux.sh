# @Author: Uriel Sandoval
# @Date:   2021-03-26 18:10:10
# @Last Modified by:   Uriel Sandoval
# @Last Modified time: 2021-03-26 18:48:08


# Get OSQP
git clone --recursive https://github.com/oxfordcontrol/osqp.git
cd osqp
git checkout 0.6.2


# Get OSQP
cd osqp
mkdir build
cd build
cmake ..
sudo cmake --build . --target install
sudo ldconfig -v


