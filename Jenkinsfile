pipeline {
  agent any
  stages {
    stage('Install Dependencies') {
      steps {
        sh '''python -m pip install --upgrade pip
pip install pytest pytest-cov coveralls numpy'''
      }
    }

    stage('Install libraries') {
      steps {
        sh '''sudo apt-get install libopenblas-dev libfftw3-dev libglpk-dev libdsdp-dev libgsl0-dev
'''
      }
    }

  }
  environment {
    KVXOPT_BUILD_GSL = '1'
  }
}