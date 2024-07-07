#!/bin/bash

gcc_env() {
    export PREREQ_DIR=${PREREQ_DIR:-"${PWD}/opt-gcc"}
    export INSTALL_DIR=${INSTALL_DIR:-"${PWD}/gcc-prereq-install"}
    export GCC_VERSION=$(gcc -dumpfullversion | awk -F. '{print $3+100*($2+100*$1)}')
    export BRAMS_INSTALL_DIR=${BRAMS_INSTALL_DIR:-"${HOME}/brams-6.0-gcc"}
}

intel_env() {
    export INTEL_COMPILER_VERSION=${INTEL_COMPILER_VERSION:-"2022.2.0"}
    export INTEL_KIT_PATH=${INTEL_KIT_PATH:-"${PWD}/"}
    export INTEL_COMPILER_DIR=${INTEL_COMPILER_DIR:-"${PWD}/inteloneapi-${INTEL_COMPILER_VERSION}/"}
}

ifort_env() {
    intel_env
    export BRAMS_INSTALL_DIR=${BRAMS_INSTALL_DIR:-"${HOME}/brams-6.0-ifort"}
    export PREREQ_DIR=${PREREQ_DIR:-"${PWD}/opt-intel-${INTEL_COMPILER_VERSION}"}
    export INSTALL_DIR=${INSTALL_DIR:-"${PWD}/intel-${INTEL_COMPILER_VERSION}-prereq-install"}
}

ifort_impi_env() {
    intel_env
    export BRAMS_INSTALL_DIR=${BRAMS_INSTALL_DIR:-"${HOME}/brams-6.0-ifort-impi"}
    export PREREQ_DIR=${PREREQ_DIR:-"${PWD}/opt-intel-impi-${INTEL_COMPILER_VERSION}"}
    export INSTALL_DIR=${INSTALL_DIR:-"${PWD}/intel-impi-${INTEL_COMPILER_VERSION}-prereq-install"}
}

export PREREQ_DL_DIR=${PREREQ_DL_DIR:-"${PWD}/prereq-download-dir/"}
export NUM_MAKE_JOBS=${NUM_MAKE_JOBS:-8}
export MPICH_VERSION=${MPICH_VERSION:-3.3.1}
export BRAMS_DIR=${BRAMS_DIR:-"${PWD}/../"}

case $1 in
  "ifort")
    echo Installing with ifort
    ifort_env

    . install-inteloneapi-offline.bash
    . install-prereq-ifort.bash
  ;;
  "ifort-impi")
    echo Installing with ifort and intel mpi
    ifort_impi_env

    . install-inteloneapi-offline.bash
    . install-prereq-ifort-impi.bash
  ;;
  "gcc")
    echo Installing with gcc
    gcc_env

    if [[ ${GCC_VERSION} -ge 100000 ]]
    then
        . install-prereq-gcc10+.bash
    else
        . install-prereq-gcc.bash
    fi
   ;;
  "download")
     echo Downloading prerequisites source files
    . download_prereq.bash
  ;;
  "brams-ifort")
     echo Installing BRAMS with ifort
     ifort_env

     . install-brams-ifort.bash
  ;;
  "brams-ifort-impi")
     echo Installing BRAMS with ifort and intel mpi
     ifort_impi_env

     . install-brams-ifort-impi.bash
  ;;
  "brams-gcc")
     echo Installing BRAMS with gcc
     gcc_env

     if [[ ${GCC_VERSION} -ge 100000 ]]
     then
        . install-brams-gcc10+.bash
     else
        . install-brams-gcc.bash
     fi
  ;;

esac