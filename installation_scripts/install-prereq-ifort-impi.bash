#!/bin/bash

INTEL_COMPILER_VERSION=${INTEL_COMPILER_VERSION:-"2022.2.0"}
INTEL_COMPILER_DIR=${INTEL_COMPILER_DIR:-"${PWD}/inteloneapi-${INTEL_COMPILER_VERSION}/"}
PREREQ_DL_DIR=${PREREQ_DL_DIR:-"${PWD}/prereq-download-dir/"}
PREREQ_DIR=${PREREQ_DIR:-"${PWD}/opt-intel-impi-${INTEL_COMPILER_VERSION}"}
INSTALL_DIR=${INSTALL_DIR:-"${PWD}/intel-impi-${INTEL_COMPILER_VERSION}-prereq-install"}
NUM_MAKE_JOBS=${NUM_MAKE_JOBS:-8}

source ${INTEL_COMPILER_DIR}/setvars.sh
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
  CC=icc ./configure --prefix=${PREREQ_DIR}
  make clean && make -j ${NUM_MAKE_JOBS}
  make install
  [[ $? -ne 0 ]] && { echo "Error while installing zlib" ; exit 1 ; }
  touch ../.zlib.done
  cd ..
fi

if [[ ! -f .szip.done ]]
then
  cp ${PREREQ_DL_DIR}/szip-2.1.tar.gz .
  tar -xzvf szip-2.1.tar.gz
  cd szip-2.1/
  CC=icc ./configure --prefix=${PREREQ_DIR}
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -ne 0 ]] && { echo "Error while installing szip" ; exit 1 ; }
  touch ../.szip.done
  cd ..
fi

if [[ ! -f .curl.done ]]
then
  cp ${PREREQ_DL_DIR}/curl-7.26.0.tar.gz .
  tar -xzvf curl-7.26.0.tar.gz
  cd curl-7.26.0/
  CC=icc ./configure --prefix=${PREREQ_DIR} --without-libssh2
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -ne 0 ]] && { echo "Error while installing curl" ; exit 1 ; }
  touch ../.curl.done
  cd ..
fi

if [[ ! -f .hdf5-impi.done ]]
then
  cp ${PREREQ_DL_DIR}/hdf5-1_12_1.tar.gz .
  tar -xzvf hdf5-1_12_1.tar.gz
  cd hdf5-hdf5-1_12_1/
  ./configure --prefix=${PREREQ_DIR} CC=mpiicc FC=mpiifort \
   --with-zlib=${PREREQ_DIR} --with-szlib=${PREREQ_DIR} --enable-parallel --enable-fortran
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -ne 0 ]] && { echo "Error while installing hdf5" ; exit 1 ; }
  touch ../.hdf5-impi.done
  cd ..
fi

if [[ ! -f .netcdf-c-impi.done ]]
then
  cp ${PREREQ_DL_DIR}/netcdf-c-4.8.1.tar.gz .
  tar -xzvf netcdf-c-4.8.1.tar.gz
  cd netcdf-c-4.8.1/
  CPPFLAGS=-I${PREREQ_DIR}/include LDFLAGS=-L${PREREQ_DIR}/lib CFLAGS='-O3'  CC=mpiicc \
   ./configure --prefix=${PREREQ_DIR} --enable-netcdf4 --enable-shared --enable-dap
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -ne 0 ]] && { echo "Error while installing netcdf-c" ; exit 1 ; }
  touch ../.netcdf-c-impi.done
  cd ..
fi

if [[ ! -f .netcdf-fortran-impi.done ]]
then
  cp ${PREREQ_DL_DIR}/netcdf-fortran-4.5.3.tar.gz .
  tar -xzvf netcdf-fortran-4.5.3.tar.gz
  cd netcdf-fortran-4.5.3/
  CPPFLAGS=-I${PREREQ_DIR}/include LDFLAGS=-L${PREREQ_DIR}/lib CFLAGS='-O3' FC=mpiifort \
    CC=mpiicc ./configure --prefix=${PREREQ_DIR}
  make clean && make -j $NUM_MAKE_JOBS
  make install
  [[ $? -ne 0 ]] && { echo "Error while installing netcdf-fortran" ; exit 1 ; }
  touch ../.netcdf-fortran-impi.done
  cd ..
fi

if [[ ! -f .wgrib2-impi.done ]]
then
  cp ${PREREQ_DL_DIR}/wgrib2.tgz .
  tar -xzvf wgrib2.tgz
  cd grib2

  sed -i 's/^USE_NETCDF3=.*/USE_NETCDF3=0/' makefile && \
  sed -i 's/^USE_NETCDF4=.*/USE_NETCDF4=0/' makefile && \
  sed -i 's/^USE_REGEX=.*/USE_REGEX=1/' makefile && \
  sed -i 's/^USE_TIGGE=.*/USE_TIGGE=1/' makefile && \
  sed -i 's/^USE_MYSQL=.*/USE_MYSQL=0/' makefile && \
  sed -i 's/^USE_IPOLATES=.*/USE_IPOLATES=3/' makefile && \
  sed -i 's/^USE_SPECTRAL=.*/USE_SPECTRAL=0/' makefile && \
  sed -i 's/^USE_UDF=.*/USE_UDF=0/' makefile && \
  sed -i 's/^USE_OPENMP=.*/USE_OPENMP=0/' makefile && \
  sed -i 's/^USE_PROJ4=.*/USE_PROJ4=0/' makefile && \
  sed -i 's/^USE_WMO_VALIDATION=.*/USE_WMO_VALIDATION=0/' makefile && \
  sed -i 's/^DISABLE_TIMEZONE=.*/DISABLE_TIMEZONE=0/' makefile && \
  sed -i 's/^USE_NAMES=NCE.*/USE_NAMES=NCEP/' makefile && \
  sed -i 's/^MAKE_FTN_API=.*/MAKE_FTN_API=1/' makefile && \
  sed -i 's/^DISABLE_ALARM=.*/DISABLE_ALARM=0/' makefile && \
  sed -i 's/^MAKE_SHARED_LIB=.*/MAKE_SHARED_LIB=0/' makefile && \
  sed -i 's/^USE_G2CLIB=.*/USE_G2CLIB=0/' makefile && \
  sed -i 's/^USE_PNG=.*/USE_PNG=0/' makefile && \
  sed -i 's/^USE_JASPER=.*/USE_JASPER=0/' makefile && \
  sed -i 's/^USE_OPENJPEG=.*/USE_OPENJPEG=0/' makefile && \
  sed -i 's/^USE_AEC=.*/USE_AEC=0/' makefile

  make CC=icc FC=ifort clean
  make CC=icc FC=ifort
  make CC=icc FC=ifort lib

  cp wgrib2/wgrib2 ${PREREQ_DIR}/bin/
  cp wgrib2/libwgrib2.a ${PREREQ_DIR}/lib/
  cp ./lib/*.a ${PREREQ_DIR}/lib/
  cp ./lib/*.mod ${PREREQ_DIR}/include/
  [[ $? -ne 0 ]] && { echo "Error while installing wgrib2" ; exit 1 ; }

  touch ../.wgrib2-impi.done
fi

exit 0