# BRAMS 6.0 - INSTALL PREREQUISITES

## Read before you begin

This prerequisite installation manual is to informing you about installing packages required to build the BRAMS model. Many systems have the pre-requirements installed and you must point to the correct directories or use the references provided in their user manuals. The BRAMS template uses configure to mount makefiles. Thus, if your system has the libraries installed, just configure it properly according to the user's manuals. This manual guides you through the installation of the following libraries:

* **MPICH**

* **HDF5** (Required by NetCDF)

* **NetCDF** (C and Fortran)

* **GRIB2**

These libraries have some requirements such as zlib, szip, and curl. Most Linux and Unix systems have these libraries installed and you will not need to proceed with the installation defined in this manual. But it is necessary to adapt them according to the address of each.

> BRAMS can work without the GRIB2 (Wgrib2) library, however some template startup features may be unavailable. For initial tests with pre-generated data this library is not required.

The installation below addresses the pre requirements with the GCC/GFortran (GNU) compiler, but you can use other compilers such as INTEL, NVDIA, or any other. Just make sure that your libraries are built with the same compilers to avoid build problems or when you run the template.

In general, the machine's native build systems are more computationally efficient by running faster and allowing substantial gains. There is the possibility of stressing the compilation options by improving vectorization, unroling, use of mathematical or memory resources. All this may be able to give an extra gain in speed. But always keep in mind that the greater the depth of the options, you can obtained unespected results.

Before you install the prerequisites you must have at least one version of gfortran and gcc compiler. The instalation depends of the version of Linux and Kernel you have.

> We recomend to install the **gfortran revision 8.4.0** and **gcc revision 8.4.0**. You can install or use newest version of compilers but is by your own risk. You may try to use INTEL, NVIDIA or other compiler too. The strucuture of setup is almost the same.

> In body of text below we will use for your C compiler choice and for your Fortran compiler choice.

In most system there are a lot of libraries pre-installed. Yes, if possible use them. But, unfortunatelly, some times a library require other dependencies and the configure may be broken. The instructions below in this documents will guide you to install in a particular way. If you got others problems please, contact your system manager.

If you wanto more information about gfortran, please, see the manual [https://gcc.gnu.org/onlinedocs/gcc-8.4.0/gfortran.pdf](https://gcc.gnu.org/onlinedocs/gcc-8.4.0/gfortran.pdf)

> Notice:
> If you are installing in a fresh system maybe you need install basic packages from system repository: M4 and  make.

## 1.  Create the prerequisites folder's structure
---------------------

   You will need a folder (**{YOUR\_DIR}**) to put the libraries, includes and executables. You can choose a system folder or a local home folder. If you prefer a system folder remember that commands must be preceeded by a sudo command. Another folder (**{INSTALL\_DIR}**) will be used only for create the sctrucuture and can be erased after all instalation:

```bash
        mkdir {YOUR_DIR}
        mkdir {YOUR_DIR}/bin
        mkdir {YOUR_DIR}/lib
        mkdir {YOUR_DIR}/include
        mkdir {INSTALL_DIR}
        cd {INSTALL_DIR}
```

   In the example below we are using system folder **/opt/apps** as {YOUR\_DIR} and user home folder **/home/oscar** as {INSTALL\_DIR}:

```bash
        sudo mkdir /opt/apps
        sudo mkdir /opt/apps/bin
        sudo mkdir /opt/apps/lib
        sudo mkdir /opt/apps/include
        mkdir /home/oscar/install
        cd /home/oscar/install
```

   > If you will use differents compilers environment we recommend do create the folder with respective names like /opt/nvidia, /opt/intel or /opt/gnu instead /opt/apps as example above.

## 2.  Download all necessary prerequisites  
---------------
   In order to build full prerequisites a bunch of packages must be downloaded. You can do it in just one command by getting the packages from FTP BRAMS' area or getting every package individually.

   To get just one file (easyest way) with all libraries please, get it as show below:

```bash
        wget http://ftp.cptec.inpe.br/pesquisa/bramsrd/BRAMS/pre-requisites/prerequisites.tar
        tar -xvf prerequisites.tar
```

 | Package        | Site                                                   | File                        |
   |:--------------:|:------------------------------------------------------:|:---------------------------:|
   | mpich3         | <http://www.mpich.org/static/downloads/4.0a2/>           | mpich-4.0a2.tar.gz          |
   | zlib           | <ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/>        | zlib-1.2.8.tar.gz           |
   | szip           | <ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/>        | szip-2.1.tar.gz             |
   | curl           | <ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/>        | curl-7.26.0.tar.gz          |
   | netCDF C       | <https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/src/> | netcdf-c-4.8.1.tar.gz       |
   | netCDF Fortran | <https://www.unidata.ucar.edu/downloads/netcdf/ftp/>     | netcdf-fortran-4.5.3.tar.gz |
   | grib2          | <https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/>       | wgrib2.tgz                  |
   | hdf5           | <https://www.hdfgroup.org/downloads/hdf5/source-code/>   | hdf5-1.12.1.tar.gz

> Notices:  
> 1. The grib2 file has diferent names if you get it from NCEP or if you get it from BRAMS ftp.  
> 2. The packages may have new version. You can try to install the new versions. Check the web page of each one.

## 3.  Enviroment variables
---------------------
   You need to be sure about where are the libraries, includes and binaries you will build the system. Is necessary to make a symbolic link between the gfortran compiler (version 8) of your system using just "gfortran" in your bin area. The same for gcc. Let's setup the paths:

```bash
        export PATH={YOUR_DIR}/bin:$PATH
        export LD_LIBRARY_PATH={YOUR_DIR}/lib:$LD_LIBRARY_PATH
        ln -s /usr/bin/gfortran {YOUR_DIR}/bin/gfortran
        ln -s /usr/bin/gcc {YOUR_DIR}/bin/gcc
 ```

  > If you will use other compiler instead gfortran/gcc make the export for correct compiler.
  > Please, check if your path have in first part {YOUR\_DIR} and if the correct library is in first part of LD\_LIBRARY\_PATH.

```bash
        echo $PATH
        echo $LD_LIBRARY_PATH
```

  Please check if fortran and gcc version is the correct or make the same for the other compiler ou will use.

```bash
        gfortran --version
        gcc --version
```

The result must be something like. See that in this case we use 8.4.0.

        $ gfortran --version
        GNU Fortran (Ubuntu 8.4.0-4ubuntu1) 8.4.0
        Copyright (C) 2018 Free Software Foundation, Inc.
        This is free software; see the source for copying conditions.  There is NO
        warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
        $ gcc --version
        gcc (Ubuntu 8.4.0-4ubuntu1) 8.4.0
        Copyright (C) 2018 Free Software Foundation, Inc.
        This is free software; see the source for copying conditions.  There is NO
        warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

   > Notice: If you create scripts for install prerequisites maybe you need use full address for compilers and libraries.

## 4.  Building the mpich libraries and binaries
---------------------

   Now you must build all libraries one by one. Let's start with mpich:  

```bash
        tar -xzvf mpich-4.0a2.tar.gz 
        cd mpich-4.0a2/
        ./configure -disable-fast CC=<Your_C_compiler> FC=<Your_fortran_compiler> CFLAGS=-O2 FFLAGS=-O2 CXXFLAGS=-O2 FCFLAGS=-O2 --prefix={YOUR_DIR} --with-device=ch3
        make
        make install
        cd ..
```

   After the install check if the mpich is installed and all the library is present.

```bash
    ls {YOUR_DIR}/bin
```

    .rwxr-xr-x 593k root 19 Apr 15:44  hydra_nameserver
    .rwxr-xr-x 587k root 19 Apr 15:44  hydra_persist
    .rwxr-xr-x 819k root 19 Apr 15:44  hydra_pmi_proxy
    lrwxrwxrwx    6 root 19 Apr 15:44  mpic++ -> mpicxx
    .rwxr-xr-x  11k root 19 Apr 15:44  mpicc
    .rwxr-xr-x  17k root 19 Apr 15:44  mpichversion
    .rwxr-xr-x  11k root 19 Apr 15:44  mpicxx
    lrwxrwxrwx   13 root 19 Apr 15:44  mpiexec -> mpiexec.hydra
    .rwxr-xr-x 1,1M root 19 Apr 15:44  mpiexec.hydra
    lrwxrwxrwx    7 root 19 Apr 15:44  mpif77 -> mpifort
    lrwxrwxrwx    7 root 19 Apr 15:44  mpif90 -> mpifort
    .rwxr-xr-x  14k root 19 Apr 15:44  mpifort
    lrwxrwxrwx   13 root 19 Apr 15:44  mpirun -> mpiexec.hydra
    .rwxr-xr-x  35k root 19 Apr 15:44  mpivars
    .rwxr-xr-x 3,5k root 19 Apr 15:44  parkill

To check the version (compiler and options) of mpich:

```bash
    {YOUR_DIR}/bin/mpirun --version
```

    HYDRA build details:
        Version:                                 4.0b1
        Release Date:                            Mon Nov 15 10:22:52 CST 2021
        CC:                              gcc-8   -O2   
        Configure options:                       
        '--disable-option-checking' '--prefix=/home/lufla/apps' '-disable-fast' 'CC=gcc-8' 'FC=gfortran-8' 'CFLAGS=-O2 -O0' 'FFLAGS=-O2 -O0' 'CXXFLAGS=-O2 -O0' 'FCFLAGS=-O2 -O0' '--with-device=ch3' '--cache-file=/dev/null' '--srcdir=.' 'LDFLAGS=' 'LIBS=' 'CPPFLAGS=-D__HIP_PLATFORM_AMD__ -I/home/lufla/install/mpich-4.0b1/src/mpl/include -I/home/lufla/install/mpich-4.0b1/modules/json-c -D_REENTRANT -I/home/lufla/install/mpich-4.0b1/src/mpi/romio/include'
        Process Manager:                         pmi
        Launchers available:                     ssh rsh fork slurm ll lsf sge manual persist
        Topology libraries available:            hwloc
        Resource management kernels available:   user slurm ll lsf sge pbs cobalt
        Demux engines available:                 poll select

## 5.  Building zlib, szip and curl
---------------------

### zlib

```bash
        tar -xzvf zlib-1.2.8.tar.gz 
        cd zlib-1.2.8/
        CC=<Your_C_compiler> ./configure --prefix={YOUR_DIR}
        make
        make install
        cd ..
```

### szip

```bash
        tar -xzvf szip-2.1.tar.gz 
        cd szip-2.1/
        CC=<Your_C_compiler> ./configure --prefix={YOUR_DIR}
        make
        make install
        cd ..
```

### curl

```bash
        tar -xzvf curl-7.26.0.tar.gz 
        cd curl-7.26.0/
        CC=<Your_C_compiler> ./configure --prefix={YOUR_DIR} --without-libssh2  --with-openssl
        make
        make install
        cd ..
```

## 6.  Building HDF5
---------------------
```bash
        tar -xzvf hdf5-1.12.1.tar.gz 
        cd hdf5-1.12.1/
        ./configure --prefix={YOUR_DIR} CC={YOUR_DIR}/bin/mpicc FC={YOUR_DIR}/bin/mpif90 --with-zlib={YOUR_DIR} --with-szlib={YOUR_DIR} --enable-parallel --enable-fortran
        make
        make install
        cd ..
```

    \--with-openssl If you used **NVIDIA** (pgf90/pgcc) to build mpich, instead the configure above, make the configure using the command lines below

```bash
        ./configure --prefix={YOUR_DIR} FFLAGS=-fPIC FCFLAGS=-fPIC CC={YOUR_DIR}/bin/mpicc FC={YOUR_DIR}/bin/mpif90 --with-zlib={YOUR_DIR} --with-szlib={YOUR_DIR} --enable-parallel --enable-fortran
        make
        make install
        cd ..--with-openssl
```

## 7.  Building NetCDF libraries
---------------------

```bash
        tar -xzvf netcdf-c-4.8.1.tar.gz 
        cd netcdf-c-4.8.1/
        CPPFLAGS=-I{YOUR_DIR}/include LDFLAGS=-L{YOUR_DIR}/lib CFLAGS='-O3'  CC={YOUR_DIR}/bin/mpicc ./configure --prefix={YOUR_DIR} --enable-netcdf4 --enable-shared --enable-dap
        make
        make install
        cd ..
```

```bash
        tar -xzvf netcdf-fortran-4.5.3.tar.gz
        cd netcdf-fortran-4.5.3/
        CPPFLAGS=-I{YOUR_DIR}/include LDFLAGS=-L{YOUR_DIR}/lib CFLAGS='-O3' FC={YOUR_DIR}/bin/mpif90  CC={YOUR_DIR}/bin/mpicc ./configure --prefix={YOUR_DIR}
        make
        make install
        cd ..
```

   > Notice: Some system may require install M4 package.

   If you are using NVIDIA (pgf90/pgcc), to build the NetCDF-fortran, instead the configure above, make the configure using the command lines below:

```bash
        FFLAGS=-fPIC FCFLAGS=-fPIC CPPFLAGS=-I{YOUR_DIR}/include LDFLAGS=-L{YOUR_DIR}/lib CFLAGS='-O3' FC={YOUR_DIR}/bin/mpif90  CC={YOUR_DIR}/bin/mpicc ./configure --prefix={YOUR_DIR}
        make
        make install
        cd ..
```

## 8.  Building grib2 libraries and API
---------------------

```bash
       tar -xzvf grib2.tgz
       cd grib2 
 ```

   Before compile you must modify the makefile. Edit the makefile, find and change the following variables (use the values as show below):

    USE_NETCDF3=0
    USE_NETCDF4=0
    USE_REGEX=1
    USE_TIGGE=1
    USE_MYSQL=0
    USE_IPOLATES=3
    USE_SPECTRAL=0
    USE_UDF=0
    USE_OPENMP=0
    USE_PROJ4=0
    USE_WMO_VALIDATION=0
    DISABLE_TIMEZONE=0
    USE_NAMES=NCEP
    MAKE_FTN_API=1
    DISABLE_ALARM=0
    MAKE_SHARED_LIB=0
    
    USE_G2CLIB=0
    USE_PNG=0
    USE_JASPER=0
    USE_OPENJPEG=0
    USE_AEC=0

   If you are using **NVIDIA compilers** (pgf90/pgcc) you must change other parts of makefile. Find the commands sequence:

```bash
       system=$(shell uname -s)
       ifeq (${system},Linux)
          ifeq ($(findstring gcc,$(notdir $(CC))),gcc)
             COMP_SYS=gnu_linux
          endif
          ifeq ($(findstring g95,$(notdir $(FC))),g95)
             COMP_SYS=gnu_linux_g95
          endif
          ifeq ($(findstring clang,$(notdir $(CC))),clang)
             COMP_SYS=clang_linux
          endif
          ifeq ($(findstring icc,$(notdir $(CC))),icc)
             COMP_SYS=intel_linux
          endif
          ifeq ($(findstring nvc,$(notdir $(CC))),nvc)
             COMP_SYS=nvidia_linux
          endif
       endif
```

   Now add the lines below at end of the sequence you found

```bash
       ifeq ($(findstring pgcc,$(notdir $(CC))),pgcc)
         COMP_SYS=nvidia_linux
       endif
```

   The result must be:

```bash
       system=$(shell uname -s)
       ifeq (${system},Linux)
          ifeq ($(findstring gcc,$(notdir $(CC))),gcc)
             COMP_SYS=gnu_linux
          endif
          ifeq ($(findstring g95,$(notdir $(FC))),g95)
             COMP_SYS=gnu_linux_g95
          endif
          ifeq ($(findstring clang,$(notdir $(CC))),clang)
             COMP_SYS=clang_linux
          endif
          ifeq ($(findstring icc,$(notdir $(CC))),icc)
             COMP_SYS=intel_linux
          endif
          ifeq ($(findstring nvc,$(notdir $(CC))),nvc)
             COMP_SYS=nvidia_linux
          endif
          ifeq ($(findstring pgcc,$(notdir $(CC))),pgcc)
             COMP_SYS=nvidia_linux
            endif
       endif
```

   Find the following sequence:

```bash
    ifeq ($(need_ftn),1)
      ifeq ($(findstring gfortran,$(notdir $(FC))),gfortran)
        a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 1`\" >> ${CONFIG_H})
      else ifeq ($(findstring ifort,$(notdir $(FC))),ifort)
        a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 1`\" >> ${CONFIG_H})
      else ifeq ($(findstring nvfortran,$(notdir $(FC))),nvfortran)
        a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 2 | tail -n 1`\" >> ${CONFIG_H})
      else
        a:=$(shell echo '$Hdefine FORTRAN '\"${FC}\" >> ${CONFIG_H})
      endif
    endif
```

   add the sequence lines after nvfortran section:

```bash
    else ifeq ($(findstring pgf90,$(notdir $(FC))),pgf90)
       a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 2 | tail -n 1`\" >> ${CONFIG_H})
```

   The results must be:

```bash
    ifeq ($(need_ftn),1)
      ifeq ($(findstring gfortran,$(notdir $(FC))),gfortran)
        a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 1`\" >> ${CONFIG_H})
      else ifeq ($(findstring ifort,$(notdir $(FC))),ifort)
        a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 1`\" >> ${CONFIG_H})
      else ifeq ($(findstring nvfortran,$(notdir $(FC))),nvfortran)
        a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 2 | tail -n 1`\" >> ${CONFIG_H})
      else ifeq ($(findstring pgf90,$(notdir $(FC))),pgf90)
        a:=$(shell echo '$Hdefine FORTRAN '\"`${FC} --version | head -n 2 | tail -n 1`\" >> ${CONFIG_H})
      else
        a:=$(shell echo '$Hdefine FORTRAN '\"${FC}\" >> ${CONFIG_H})
      endif
    endif
```

   With the modifications you can now compile and install:

```bash
    make CC=<Your_C_compiler> FC=<Your_Fortran_compiler>
    make CC=<Your_C_compiler> FC=<Your_Fortran_compiler> lib  
    cp wgrib2/wgrib2 {YOUR_DIR}/bin
    cp ./lib/*.a {YOUR_DIR}/lib
    cp ./lib/*.mod {YOUR_DIR}/include
```

## 9.  Creating alias for use
---------------------

   If you compile for a specific compiler, gnu, nvidia or intel (or others) all the libraries and binaries will be build in the compiler. For correct use keep in mind that the PATHs must be pointed to correct binary. For example: you will need the mpirun binary for run the model. We recommend you create some alias to call before the use of each run model. Edit the .bashrc in your home area and put the folowwing lines:

```bash
    alias gnu8='export PATH=/opt/gnu8/bin:$PATH;export LD_LIBRARY_PATH=/opt/gnu8/lib:$LD_LIBRARY_PATH'
    alias nvidia='export PATH=/opt/nvidia/bin:$PATH;export LD_LIBRARY_PATH=/opt/nvidia/lib:$LD_LIBRARY_PATH'
```

   In the example we create two alias, the first "gnu8" adjust the path to search in /opt/gnu8/bin first and use the libs from /opt/gnu8/lib. At same way the nvidia alias but change the pointer to nvidia area.

   After edit and save just type on terminal:  

```bash
    source ~/.bashrc 
```

   and the two alias will be available. Now you can choose the correct before run the model. If you will run the gnu install execute the alias gnu8 on terminal and you will be ready to go.

* * *

  All the required prerequisites are installed. Have fun!
