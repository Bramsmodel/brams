#!/bin/bash

# Basekit and HPCKit Linux's offline installer urls
#https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=linux&distributions=offline
#https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&distributions=offline

INTEL_COMPILER_VERSION=${INTEL_COMPILER_VERSION:-"2022.2.0"}
INTEL_KIT_PATH=${INTEL_KIT_PATH:-"${PWD}/"}
INTEL_COMPILER_DIR=${INTEL_COMPILER_DIR:-"${PWD}/inteloneapi-${INTEL_COMPILER_VERSION}/"}

if [[ ! -f  ${INTEL_COMPILER_DIR}/.intelbasekit.done ]]
then
  sh ${INTEL_KIT_PATH}/l_BaseKit_p_${INTEL_COMPILER_VERSION}*_offline.sh \
   -r yes -a --silent --eula accept --install-dir ${INTEL_COMPILER_DIR} --instance "${INTEL_COMPILER_VERSION}" \
   --components intel.oneapi.lin.mkl.devel:intel.oneapi.lin.dpcpp_dbg

  [[ $? -eq 0 ]] || { echo "Error while installing Intel BaseKit" ; exit 1 ; }
  touch ${INTEL_COMPILER_DIR}/.intelbasekit.done
else
  echo Intel BaseKit is already installed
fi

if [[ ! -f  ${INTEL_COMPILER_DIR}/.intelhpckit.done ]]
then
  sh ${INTEL_KIT_PATH}/l_HPCKit_p_${INTEL_COMPILER_VERSION}*_offline.sh \
   -r yes -a --silent --eula accept --install-dir ${INTEL_COMPILER_DIR} --instance "${INTEL_COMPILER_VERSION}" \
   --components intel.oneapi.lin.mpi.devel:intel.oneapi.lin.dpcpp-cpp-compiler-pro:intel.oneapi.lin.ifort-compiler
  [[ $? -eq 0 ]] || { echo "Error while installing Intel HPCKit" ; exit 1 ; }
  touch ${INTEL_COMPILER_DIR}/.intelhpckit.done
else
  echo Intel HPCKit is already installed
fi