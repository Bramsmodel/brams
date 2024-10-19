#!/bin/bash
# 
# Load "spack load openmpi%gcc@12.2.0 gcc@12.2.0" via "source LoadBramsEnv.bash"
#
BASE_DIR=/scr2-exa/panetta/Furnas
YOUR_DIR=${YOUR_DIR:-"${BASE_DIR}/UtilsBRAMS_SPACK_12.1"}

BRAMS_DIR=${BRAMS_DIR:-"${BASE_DIR}/BuildParalelismo/B7LimpoSrc"}
RUN_DIR=${BASE_DIR}/BuildParalelismo/runB7Para

TAR_DIR=${TAR_DIR:-"${YOUR_DIR}/tgz"}
EXEC_NAME=brams-6.0

export PATH=${YOUR_DIR}/bin:$PATH
export LD_LIBRARY_PATH=${YOUR_DIR}/lib:$LD_LIBRARY_PATH

echo "*** which mpif90 ***"
which mpif90
echo "*** mpif90 --version ***"
mpif90 --version
echo "*** which gfortran ***"
which gfortran
echo "*** gfortran --version ***"
gfortran --version
echo "*** which gcc ***"
which gcc
echo "*** gcc --version ***"
gcc --version
echo "*** Brams Utils from $YOUR_DIR ***"
echo ""

mpich=0
zlib=0
szip=0
curl=0
hdf5=0
netcdfc=0
netcdff=0
grib2=0
bramsFromScratch=1
brams=1

if [ $mpich = 1 ]; then
    echo "building mpich"
    TAR_FILE=mpich-4.0a2.tar.gz 
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/mpich-4.0a2
    cd ${BUILD_DIR}
    echo "configure"
    ./configure -disable-fast CC=gcc FC=gfortran CFLAGS=-O2 FFLAGS=-O2 CXXFLAGS=-O2 FCFLAGS=-O2 --prefix=${YOUR_DIR} --with-device=ch3 &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
	exit 1
    fi
    echo "make"
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make install"
    make install &> make_install.out
    if [ $? != 0 ]; then
	echo "make install failed; check "${BUILD_DIR}"/make_install.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi



if [ $zlib = 1 ]; then
    echo "building zlib"
    TAR_FILE=zlib-1.2.8.tar.gz
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/zlib-1.2.8/
    cd ${BUILD_DIR}
    echo "configure"
    CC=gcc ./configure --prefix=${YOUR_DIR} &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
	exit 1
    fi
    echo "make"
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make install"
    make install &> make_install.out
    if [ $? != 0 ]; then
	echo "make_install failed; check "${BUILD_DIR}"/make_install.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi



if [ $szip = 1 ]; then
    echo "building szip"
    TAR_FILE=szip-2.1.tar.gz
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/szip-2.1
    cd ${BUILD_DIR}
    echo "configure"
    CC=gcc ./configure --prefix=${YOUR_DIR} &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
	exit 1
    fi
    echo "make"
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make install"
    make install &> make_install.out
    if [ $? != 0 ]; then
	echo "make_install failed; check "${BUILD_DIR}"/make_install.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi



if [ $curl = 1 ]; then
    echo "building curl"
    TAR_FILE=curl-7.26.0.tar.gz
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/curl-7.26.0
    cd ${BUILD_DIR}
    echo "configure"
    CC=gcc ./configure --prefix=${YOUR_DIR} --without-libssh2 &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
	exit 1
    fi
    echo "make"
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make install"
    make install &> make_install.out
    if [ $? != 0 ]; then
	echo "make_install failed; check "${BUILD_DIR}"/make_install.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi

if [ $hdf5 = 1 ]; then
    echo "building hdf5"
    TAR_FILE=hdf5-1.12.1.tar.gz
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/hdf5-1.12.1
    cd ${BUILD_DIR}
    echo "configure"
    ./configure --prefix=${YOUR_DIR} CC=mpicc FC=mpif90 --with-zlib=${YOUR_DIR}/lib/ --with-szlib=${YOUR_DIR}/lib --enable-parallel --enable-fortran &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
	exit 1
    fi
    echo "make"
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make install"
    make install &> make_install.out
    if [ $? != 0 ]; then
	echo "make_install failed; check "${BUILD_DIR}"/make_install.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi

if [ $netcdfc = 1 ]; then
    echo "building netcdf-c"
    TAR_FILE=netcdf-c-4.8.1.tar.gz
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/netcdf-c-4.8.1
    cd ${BUILD_DIR}
    echo "configure"
    CPPFLAGS=-I${YOUR_DIR}/include LDFLAGS=-L${YOUR_DIR}/lib CFLAGS='-O3'  CC=mpicc ./configure --prefix=${YOUR_DIR} --enable-netcdf4 --enable-shared --enable-dap &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
	exit 1
    fi
    echo "make"
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make install"
    make install &> make_install.out
    if [ $? != 0 ]; then
	echo "make_install failed; check "${BUILD_DIR}"/make_install.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi
                                                                                                                           



if [ $netcdff = 1 ]; then
    echo "building netcdf-f"
    TAR_FILE=netcdf-fortran-4.5.3.tar.gz
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/netcdf-fortran-4.5.3
    cd ${BUILD_DIR}
    echo "configure"
    CPPFLAGS=-I${YOUR_DIR}/include LDFLAGS=-L${YOUR_DIR}/lib CFLAGS='-O3' FC=mpif90  CC=mpicc ./configure --prefix=${YOUR_DIR} &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
	exit 1
    fi
    echo "make"
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make install"
    make install &> make_install.out
    if [ $? != 0 ]; then
	echo "make_install failed; check "${BUILD_DIR}"/make_install.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi


if [ $grib2 = 1 ]; then
    echo "building grib2"
    TAR_FILE=grib2.tar.gz
    cd ${YOUR_DIR}
    cp ${TAR_DIR}/${TAR_FILE} .
    echo "tar -xzvf ${TAR_FILE} &> tar.out"
    tar -xzvf ${TAR_FILE} &> tar.out
    if [ $? != 0 ]; then
	echo "tar failed; check "${YOUR_DIR}"/tar.out" 
	exit 1
    fi
    rm -f tar.out
    rm -f ${TAR_FILE}
    BUILD_DIR=${YOUR_DIR}/grib2
    cd ${BUILD_DIR}
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
    echo "make"
    make CC=gcc FC=gfortran &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "make lib"
    make CC=gcc FC=gfortran lib &> make_lib.out
    if [ $? != 0 ]; then
	echo "make_lib failed; check "${BUILD_DIR}"/make_lib.out" 
	exit 1
    fi
    cp wgrib2/wgrib2 ${YOUR_DIR}/bin/
    cp wgrib2/libwgrib2.a ${YOUR_DIR}/lib/
    cp ./lib/*.a ${YOUR_DIR}/lib/
    cp ./lib/*.mod ${YOUR_DIR}/include/
    cd ${YOUR_DIR}
    rm -rf ${BUILD_DIR}
fi




if [ $bramsFromScratch = 1 ]; then
    echo "building brams from scratch"
    echo "source at $BRAMS_DIR"
    echo "Brams utils from $YOUR_DIR"
echo "*** Results simbolic link at $RUN_DIR ***"
    BUILD_DIR=$BRAMS_DIR/build
    cd ${BUILD_DIR}
    echo "configure > configure.out"
    rm -f configure.out
#    ./configure --program-prefix=BRAMS_6.0 --prefix=${EXEC_NAME} --enable-jules  --with-chem=RELACS_TUV --with-aer=SIMPLE --with-fpcomp=mpif90  --with-cpcomp=mpicc --with-fcomp=gfortran --with-ccomp=gcc --with-netcdff=${YOUR_DIR} --with-netcdfc=${YOUR_DIR} &> configure.out
    ./configure --program-prefix=${EXEC_NAME} --enable-jules  --with-chem=RELACS_TUV --with-aer=SIMPLE --with-fpcomp=mpif90  --with-cpcomp=mpicc --with-fcomp=gfortran --with-ccomp=gcc --with-netcdff=${YOUR_DIR} --with-netcdfc=${YOUR_DIR} &> configure.out
    if [ $? != 0 ]; then
	echo "configure failed; check "${BUILD_DIR}"/configure.out" 
        exit 1
    fi
    echo "make clean > make_clean.out"
    rm -f make_clean.out
    make clean &> make_clean.out
    if [ $? != 0 ]; then
	echo "make clean failed; check "${BUILD_DIR}"/make_clean.out" 
	exit 1
    fi
    cd ${YOUR_DIR}
    brams=1
fi

if [ $brams = 1 ]; then
    BUILD_DIR=$BRAMS_DIR/build
    cd ${BUILD_DIR}
    echo "compiling brams from source at $BRAMS_DIR"
    echo "using utils from $YOUR_DIR"
    echo "results simbolic link at $RUN_DIR"
    echo "make > make.out"
    rm -f make.out
    make &> make.out
    if [ $? != 0 ]; then
	echo "make failed; check "${BUILD_DIR}"/make.out" 
	exit 1
    fi
    echo "rm -f ${RUN_DIR}/${EXEC_NAME}"
    rm -f ${RUN_DIR}/${EXEC_NAME}
    if [ $? != 0 ]; then
	echo "rm -f ${RUN_DIR}/${EXEC_NAME} falhou" 
	exit 1
    fi
    echo "ln -s ${BUILD_DIR}/${EXEC_NAME} ${RUN_DIR}" 
    ln -s ${BUILD_DIR}/${EXEC_NAME} ${RUN_DIR}
    if [ $? != 0 ]; then
	echo "ln -s ${BUILD_DIR}/${EXEC_NAME} ${RUN_DIR} falhou" 
	exit 1
    fi
    echo "*** TERMINOU CORRETAMENTE ***"
    cd ${YOUR_DIR}
fi
