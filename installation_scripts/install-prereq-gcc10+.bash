#!/bin/bash

PREREQ_DIR=${PREREQ_DIR:-"${PWD}/opt-gcc"}
INSTALL_DIR=${INSTALL_DIR:-"${PWD}/gcc-prereq-install"}
PREREQ_DL_DIR=${PREREQ_DL_DIR:-"${PWD}/prereq-download-dir/"}
NUM_MAKE_JOBS=${NUM_MAKE_JOBS:-8}
MPICH_VERSION=${MPICH_VERSION:-3.3.1}

export PATH=${PREREQ_DIR}/bin:$PATH
export LD_LIBRARY_PATH=${PREREQ_DIR}/lib:$LD_LIBRARY_PATH

mkdir -p ${PREREQ_DIR}
mkdir -p ${PREREQ_DIR}/bin
mkdir -p ${PREREQ_DIR}/lib
mkdir -p ${PREREQ_DIR}/include
mkdir -p ${INSTALL_DIR}
cd ${INSTALL_DIR}

unset F90 F90FLAGS

if [[ ! -f .zlib.done ]]
then
  cp ${PREREQ_DL_DIR}/zlib-1.2.11.tar.gz .
  tar -xzvf zlib-1.2.11.tar.gz
  cd zlib-1.2.11/
  CC=gcc ./configure --prefix=${PREREQ_DIR}
  make clean && make -j ${NUM_MAKE_JOBS}
  make install
  [[ $? -eq 0 ]] || { echo "Error while installing zlib" ; exit 1 ; }
  touch ../.zlib.done
  cd ..
fi

if [[ ! -f .szip.done ]]
then
  cp ${PREREQ_DL_DIR}/szip-2.1.1.tar.gz .
  tar -xzvf szip-2.1.1.tar.gz
  cd szip-2.1.1/
  CC=gcc ./configure --prefix=${PREREQ_DIR}
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -eq 0 ]] || { echo "Error while installing szip" ; exit 1 ; }
  touch ../.szip.done
  cd ..
fi

if [[ ! -f .curl.done ]]
then
  cp ${PREREQ_DL_DIR}/curl-7.26.0.tar.gz .
  tar -xzvf curl-7.26.0.tar.gz
  cd curl-7.26.0/
  CC=gcc ./configure --prefix=${PREREQ_DIR} --without-libssh2
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -eq 0 ]] || { echo "Error while installing curl" ; exit 1 ; }
  touch ../.curl.done
  cd ..
fi

if [[ ! -f .mpich-${MPICH_VERSION}.done ]]
then
  cp ${PREREQ_DL_DIR}/mpich-${MPICH_VERSION}.tar.gz .
  tar -xzvf mpich-${MPICH_VERSION}.tar.gz
  cd mpich-${MPICH_VERSION}/
  ./configure -disable-fast FFLAGS="-O2 -fallow-argument-mismatch" \
    FCFLAGS="-O2 -fallow-argument-mismatch" CC=gcc FC=gfortran F77=gfortran CFLAGS=-O2 \
    CXXFLAGS=-O2  --prefix=${PREREQ_DIR} --with-device=ch3
  make clean &&  make -j $NUM_MAKE_JOBS
  make install
  [[ $? -eq 0 ]] || { echo "Error while installing mpich-${MPICH_VERSION}" ; exit 1 ; }
  touch ../.mpich-${MPICH_VERSION}.done
  cd ..
fi

if [[ ! -f .hdf5-mpich-${MPICH_VERSION}.done ]]
then
  cp ${PREREQ_DL_DIR}/hdf5-1_12_1.tar.gz .
  tar -xzvf hdf5-1_12_1.tar.gz
  cd hdf5-hdf5-1_12_1/
  ./configure --prefix=${PREREQ_DIR} FFLAGS="-fallow-argument-mismatch"  CC=${PREREQ_DIR}/bin/mpicc \
    FC=${PREREQ_DIR}/bin/mpif90 --with-zlib=${PREREQ_DIR} --with-szlib=${PREREQ_DIR} \
    --enable-parallel --enable-fortran
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -eq 0 ]] || { echo "Error while installing hdf5" ; exit 1 ; }
  touch ../.hdf5-mpich-${MPICH_VERSION}.done
  cd ..
fi

if [[ ! -f .netcdf-c-mpich-${MPICH_VERSION}.done ]]
then
  cp ${PREREQ_DL_DIR}/netcdf-c-4.8.1.tar.gz .
  tar -xzvf netcdf-c-4.8.1.tar.gz
  cd netcdf-c-4.8.1/
  CPPFLAGS=-I${PREREQ_DIR}/include LDFLAGS=-L${PREREQ_DIR}/lib CFLAGS='-O3'  CC=${PREREQ_DIR}/bin/mpicc \
   ./configure --prefix=${PREREQ_DIR} --enable-netcdf4 --enable-shared --enable-dap
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -eq 0 ]] || { echo "Error while installing netcdf-c" ; exit 1 ; }
  touch ../.netcdf-c-mpich-${MPICH_VERSION}.done
  cd ..
fi

if [[ ! -f .netcdf-fortran-mpich-${MPICH_VERSION}.done ]]
then
  cp ${PREREQ_DL_DIR}/netcdf-fortran-4.5.3.tar.gz .
  tar -xzvf netcdf-fortran-4.5.3.tar.gz
  cd netcdf-fortran-4.5.3/
  CPPFLAGS=-I${PREREQ_DIR}/include LDFLAGS=-L${PREREQ_DIR}/lib CFLAGS='-O3' FC=${PREREQ_DIR}/bin/mpif90 \
    CC=${PREREQ_DIR}/bin/mpicc FFLAGS="-fallow-argument-mismatch"  ./configure --prefix=${PREREQ_DIR}
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -eq 0 ]] || { echo "Error while installing netcdf-fortran" ; exit 1 ; }
  touch ../.netcdf-fortran-mpich-${MPICH_VERSION}.done
  cd ..
fi
