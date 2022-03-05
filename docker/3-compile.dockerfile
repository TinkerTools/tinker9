### compile tinker9
### docker build -t temp_build_tinker9 -f 3-compile.dockerfile .

ARG       T9Dir=/home/tinker9

# load the base image
FROM      tinkertools/tinker9-devel:cuda11.2-nvhpc22.2
ARG       T9Dir

# clone tinker9
WORKDIR   /home
RUN       git clone https://github.com/tinkertools/tinker9

# configure tinker9
WORKDIR   $T9Dir
RUN       git submodule update --init
# this line is important;
# however it is not required if built interactively
ENV       CUDA_HOME=/usr/local/cuda
RUN       cmake $T9Dir -B $T9Dir/build \
-DCOMPUTE_CAPABILITY=50,60,70,75,80 \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_INSTALL_PREFIX=$T9Dir/bin
ENV       PATH "$PATH:$T9Dir/bin"

# make tinker9
WORKDIR   $T9Dir/build
RUN       make -j8 VERBOSE=1
RUN       make install
