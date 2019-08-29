#!/bin/bash

# Author: Andrew Abi-Mansour
# DOC: 08/29/2019
# Email: andrew.gaam [at] gmail [dot] com
# Description: entry script for compiling/testing tinker.gpu
# Usage example: 
#	./run.sh --install --fftw_dir /usr/lib/x86_64-linux-gnu --tinker_dir ../../tinker/source
#	./run.sh --install --fftw_dir ../src/build-travis/fftw --tinker_dir ../src/build-travis/tinker/source

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
	make fftw_dir=$1 libtinker_dir=$2
	echo "Compilation successful. Now run test cases via: make test"
}

function clean {

	if [ -d "build" ]; then rm -r build && echo "Build dir recursively removed."; fi
}

function test {

	if [ ! -d "build" ]; then install $1 $2; fi

	cd build
	make test fftw_dir=/usr/lib/x86_64-linux-gnu libtinker_dir=../../tinker/source
	echo "Test cases successfully completed."
}


#!/bin/bash

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
    -w|--tinker_dir)
    TINKER_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -z|--fftw_dir)
    FFTW_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ ! -z ${BUILD} ]; then
	pre_install
else
	if [ ! -z ${INSTALL} ]; then 
		echo ${FFTW_DIR} 
		echo ${TINKER_DIR}
		install ${FFTW_DIR} ${TINKER_DIR}
	else
		if [ ! -z ${CLEAN} ]; then
			clean
		else
			if [ ! -z ${TEST} ]; then
				test ${FFTW_DIR} ${TINKER_DIR} 
			fi
		fi
	fi
fi
