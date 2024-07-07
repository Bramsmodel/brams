#!/bin/bash
# Use in BRAMS-6.0 folder
PREREQ_DIR=${PREREQ_DIR:-"${PWD}/opt-intel-impi-${INTEL_COMPILER_VERSION}"}
INTEL_COMPILER_DIR=${INTEL_COMPILER_DIR:-"${PWD}/inteloneapi-${INTEL_COMPILER_VERSION}/"}
BRAMS_DIR=${BRAMS_DIR:-"${PWD}/../"}
BRAMS_INSTALL_DIR=${BRAMS_INSTALL_DIR:-"${HOME}/brams-6.0-ifort-impi"}

source ${INTEL_COMPILER_DIR}/setvars.sh --force
export PATH=${PREREQ_DIR}/bin:$PATH
export LD_LIBRARY_PATH=${PREREQ_DIR}/lib:$LD_LIBRARY_PATH

cd ${BRAMS_DIR}/build

./configure --program-prefix=BRAMS_6.0 --prefix=${BRAMS_INSTALL_DIR} --enable-jules \
  --with-chem=RELACS_TUV --with-aer=SIMPLE --with-fpcomp=mpiifort \
   --with-cpcomp=mpiicc --with-fcomp=ifort --with-ccomp=icc \
   --with-netcdff=${PREREQ_DIR} --with-netcdfc=${PREREQ_DIR} --with-wgrib2=${PREREQ_DIR}
make clean
make FPCOMP=mpiifort CPCOMP=mpiicc CCOMP=icc FCOMP=ifort -f Make_utils
make
make install
[[ $? -eq 0 ]] || { echo "Error while installing BRAMS" ; exit 1 ;}