#!/bin/bash
# Use in BRAMS-6.0 folder

FROM_SCRATCH=0

PREREQ_DIR=/home/mateusff/opt-gcc

BRAMS_DIR=${BRAMS_DIR:-"${PWD}/../"}
BRAMS_INSTALL_DIR=`pwd`

source /home/mateusff/modules.sh

export LD_LIBRARY_PATH=${PREREQ_DIR}/lib:/opt/spack/opt/spack/linux-rocky8-zen/gcc-8.5.0/libpciaccess-0.17-u577d6vkumv2ku5u3jljfmpl54ejgkxi/lib/:$LD_LIBRARY_PATH
export PATH=${PREREQ_DIR}/bin:$PATH
#export LD_LIBRARY_PATH=${PREREQ_DIR}/lib:$LD_LIBRARY_PATH


cd $BRAMS_DIR/build
if [ $FROM_SCRATCH == 1 ]; then
    echo "FROM_SCRATCH"
    ./configure --program-prefix=BRAMS_6.0 --prefix=${BRAMS_INSTALL_DIR} --enable-jules \
		--with-chem=RELACS_TUV --with-aer=SIMPLE --with-fpcomp=mpif90  \
		--with-cpcomp=mpicc --with-fcomp=gfortran --with-ccomp=gcc \
		--with-netcdff=${PREREQ_DIR} --with-netcdfc=${PREREQ_DIR} --with-wgrib2=${PREREQ_DIR}
    make clean
fi

make

