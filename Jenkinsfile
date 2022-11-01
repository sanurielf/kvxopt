pipeline {
  agent none
  environment {
    KVXOPT_BUILD_GSL = '1'
    KVXOPT_BUILD_FFTW = '1'
    KVXOPT_BUILD_GLPK = '1'
    KVXOPT_BUILD_DSDP = '1'
    KVXOPT_BUILD_OSQP = '1'
    SUITESPARSE_VERSION = '5.10.1'
    SUITESPARSE_SHA256 = 'acb4d1045f48a237e70294b950153e48dce5b5f9ca8190e86c2b8c54ce00a7ee'
    OSQP_VERSION = '0.6.2'
    OSQP_SHA256 = '0a7ade2fa19f13e13bc12f6ea0046ef764049023efb4997a4e72a76534f623ec'
    PREFIX_LINUX = '/usr/local'
  }
  stages {
    stage('Linux build') {
      parallel {
        stage ('Test in docker') {
          agent { docker { image 'python:3' } }
          stages {

            stage ('Set environment paths') {
              steps {
                script {
                  env.PREFIX = "${PREFIX_LINUX}"
                  env.KVXOPT_OSQP_LIB_DIR = "${PREFIX}/lib"
                  env.KVXOPT_OSQP_INC_DIR = "${PREFIX}/include/osqp"
                  env.LD_LIBRARY_PATH = "${PREFIX}/lib"
                }
              }
            }

            stage('Install python dependencies') {
              steps {
                sh '''python -m pip install --upgrade pip
                      pip install --upgrade pytest pytest-cov coveralls numpy'''
              }
            }

            stage('Install libraries') {
              steps {
                sh '''apt-get update'''
                sh '''DEBIAN_FRONTEND=noninteractive apt-get -yq install libopenblas-dev libfftw3-dev libglpk-dev libdsdp-dev libgsl0-dev libsuitesparse-dev'''
              }
            }

            stage('Install cmake') {
              steps {
                sh '''apt-get update
                      DEBIAN_FRONTEND=noninteractive apt-get -yq install cmake'''
              }
            }

            stage('Get OSQP source and compile library'){
              steps {
                sh '''wget https://github.com/osqp/osqp/releases/download/v${OSQP_VERSION}/complete_sources.tar.gz -O osqp-${OSQP_VERSION}.tar.gz
                      echo "${OSQP_SHA256}  osqp-${OSQP_VERSION}.tar.gz" > OSQP.sha256
                      shasum -a 256 -c OSQP.sha256
                      tar -xf osqp-${OSQP_VERSION}.tar.gz
                      cd osqp
                      mkdir build
                      cd build
                      cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ..
                      cmake --build . --target install'''
              }
            }

            stage('Install') {
              steps {
                sh '''python setup.py build
                      pip install .'''
              }
            }

            stage('Test') {
              steps {
                sh '''python -c 'from kvxopt import blas,dsdp,lapack,glpk,osqp,fftw,gsl,cholmod,umfpack,klu'
                      pytest --cov=kvxopt tests/'''
              }
            }

          }
        }
        stage('Test in Linux server') {
          agent { label 'node_local_debian' }

          environment {
            PATH = "${HOME}/.local/bin:$PATH"
            KVXOPT_BUILD_GRB = '1'
            KVXOPT_GRB_LIB_DIR = '/opt/gurobi912/linux64/lib'
            LD_LIBRARY_PATH = '/opt/gurobi912/linux64/lib'
            KVXOPT_GRB_INC_DIR = '/opt/gurobi912/linux64/include'
            KVXOPT_GRB_LIB = 'gurobi91'
          }

          stages {

            stage ('Create virtualenv') {
              steps {
                echo "PATH is: $PATH"
                echo "SHELL is: $SHELL"
                sh '''python3 -m pip install --upgrade virtualenv
                      virtualenv --python=python3 venv
                      . venv/bin/activate
                    '''
              }
            }

            stage ('Set environment paths') {
              steps {

                sh "printenv"
                script {
                  env.PREFIX = "${env.WORKSPACE}"
                  env.KVXOPT_OSQP_LIB_DIR = "${PREFIX}/lib"
                  env.KVXOPT_OSQP_INC_DIR = "${PREFIX}/include/osqp"
                  env.LD_LIBRARY_PATH = "${PREFIX}/lib"
                  sh 'mkdir -p ${PREFIX}/lib'
                  sh 'mkdir -p ${PREFIX}/include'
                }
              }
            }

            stage('Install python dependencies') {
              steps {
                sh '''python3 -m pip install --upgrade pip
                      pip3 install --upgrade pytest pytest-cov coveralls numpy'''
              }
            }

            stage('Get OSQP source and compile library'){
              steps {
                sh '''wget https://github.com/osqp/osqp/releases/download/v${OSQP_VERSION}/complete_sources.tar.gz -O osqp-${OSQP_VERSION}.tar.gz
                      echo "${OSQP_SHA256}  osqp-${OSQP_VERSION}.tar.gz" > OSQP.sha256
                      shasum -a 256 -c OSQP.sha256
                      tar -xf osqp-${OSQP_VERSION}.tar.gz
                      cd osqp
                      mkdir build
                      cd build
                      cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ..
                      cmake --build . --target install'''
              }
            }

            stage('Install') {
              steps {
                sh '''python3 setup.py build
                      pip3 install .'''
              }
            }

            stage('Test') {
              steps {
                sh '''python3 -c 'from kvxopt import blas,dsdp,lapack,glpk,osqp,fftw,gsl,cholmod,umfpack,klu'
                      pytest --cov=kvxopt tests/'''
              }
            }

          }
          post {
            // Clean after build
            always {
                cleanWs(cleanWhenNotBuilt: false,
                        deleteDirs: true,
                        disableDeferredWipeout: true,
                        notFailBuild: true)
            }
          }
        }
      }
    }
  }

}