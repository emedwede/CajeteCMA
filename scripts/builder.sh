#!/bin/bash

arch=$(uname -m)
#set your compiler 
compiler=g++
compiler_version=$(g++ --version)
#this is the path to where the github repo is on your computer 
env_prefix=/home/erock/Documents/Code/Release
project_name=CajeteCMA #don't change me
build_name=build_test #change me to whatever file name you want
source_prefix=${env_prefix}/${project_name} #code location
build_base_dir=${env_prefix}/${build_name} #build base
build_prefix=${build_base_dir}/build/${arch}-${compiler}
#set up so that if you type `make install` it goes to a sane place
install_prefix=${build_base_dir}/install/${arch}-${compiler}
build_dir=${build_prefix}/${project_name}
install_dir=${install_prefix}/${project_name}

#sundials_dir=/home/erock/Extras/sundials-build

source_dir=${source_prefix}

#only type ./scripts/builder.sh -c if you want to clean up
#to be super safe, clean up manually
if [[ "$1" == "-c" ]]; then
    echo Removing ${build_base_dir} ...
    rm -rf ${build_base_dir}
fi

echo 'Build in ' ${build_dir}
echo 'Install in ' ${install_dir}

mkdir -p ${build_dir} && cd ${build_dir}

echo `pwd`
cmake   -D CMAKE_CXX_EXTENSIONS=Off \
        -D CMAKE_INSTALL_PREFIX=${install_path} \
        ${source_dir}
        #-D SUNDIALS_DIR=${sundials_dir} \
        #${source_dir}
