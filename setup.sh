#!/bin/bash

# Author: Andrew Abi-Mansour
# DOC: 08/29/2019
# Email: andrew.gaam [at] gmail [dot] com
# Description: entry script for compiling/testing tinker.gpu
# Usage example: 
#	export FFTW_DIR=/usr/lib/x86_64-linux-gnu
#	export TINKER_DIR=../../tinker/source
#	export ACC=pgc++
#	./setup.sh --install --fftw_dir ${FFTW_DIR} --tinker_dir ${TINKER_DIR} --acc ${ACC}

function create_build {

	if [ ! -d "build" ]; then mkdir "build"; fi
}

function copy_files {

	cp make/Makefile make/*.inc build
}

function pre_install {

	create_build
	copy_files
	echo "Installation files copied to build dir. Now run 'cd build' then 'make' to compile tinker.gpu."
}

function install {

	if [ ! -d "build" ]; then pre_install; fi

	cd build
	make fftw_include=$1 fftw_lib=$2 libtinker_dir=$3 acc=$4 host=$5
}

function clean {

	if [ -d "build" ]; then rm -r build && echo "Build dir recursively removed."; fi
}

function test {

	if [ ! -d "build" ]; then install $1 $2 $3 $4 $5; fi

	cd build
	make test fftw_include=$1 fftw_lib=$2 libtinker_dir=$3 acc=$4 host=$5
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--build)
    BUILD=1
    shift
    ;;
    -i|--install)
    INSTALL=1
    shift
    ;;
    -t|--test)
    TEST=1
    shift
    ;;
    -c|--clean)
    CLEAN=1
    shift
    ;;
    --tinker_dir)
    TINKER_DIR="$2"
    shift
    shift
    ;;
    --fftw_dir)
    FFTW_DIR="$2"
    shift
    shift
    ;;
    --fftw_include)
    FFTW_INCLUDE="$2"
    shift
    shift
    ;;
    --fftw_lib)
    FFTW_LIB="$2"
    shift
    shift
    ;;
    --acc)
    ACC="$2"
    shift
    shift
    ;;
    --host)
    HOST=1
    ;;
    --default)
    DEFAULT=YES
    shift
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# default args
# ============
if [ -z ${HOST} ]; then
	HOST=1
fi

if [ -z ${ACC} ]; then
        ACC=g++
fi

if [ -z ${FFTW_DIR} ]; then
        if [ -z ${FFTW_LIB} ]; then
		echo "Error! fftw_lib or fftw_dir must be supplied."
		exit 1
	fi

	if [ -z ${FFTW_INCLUDE} ]; then
		echo "Error! fftw_include or fftw_dir must be supplied."
                exit 1
        fi
fi

if [ -z ${FFTW_INCLUDE} ]; then
        FFTW_INCLUDE=${FFTW_DIR}/include
fi

if [ -z ${FFTW_LIB} ]; then
        FFTW_LIB=${FFTW_DIR}/lib
fi

if [ -z ${TINKER_DIR} ]; then
        TINKER_DIR=../tinker/source
fi

if [ ! -z ${BUILD} ]; then
	pre_install
else
	if [ ! -z ${INSTALL} ]; then 
		install ${FFTW_INCLUDE} ${FFTW_LIB} ${TINKER_DIR} ${ACC} ${HOST}
	else
		if [ ! -z ${CLEAN} ]; then
			clean
		else
			if [ ! -z ${TEST} ]; then
				test ${FFTW_INCLUDE} ${FFTW_LIB} ${TINKER_DIR} ${ACC} ${HOST}  
			fi
		fi
	fi
fi
