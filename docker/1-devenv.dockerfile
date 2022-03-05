### create the image for the toolchains

### command
###     docker build -t [REPO/NAME:TAG] -f [DOCKER_FILE] .
### e.g.
###     docker build -t tinkertools/tinker9-devel:cuda11.2-nvhpc22.2 -f 1-devenv.dockerfile .

### don't miss the last dot in the command

# load the base image
FROM      nvidia/cuda:11.2.2-devel-ubuntu20.04

# install packages
RUN       export DEBIAN_FRONTEND=noninteractive
RUN       apt update -y && apt install -y ca-certificates libssl-dev gfortran git cmake

# install nvhpc 22.2 multi-cuda version via apt
RUN       echo 'deb [trusted=yes] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | tee /etc/apt/sources.list.d/nvhpc.list
RUN       apt update -y && apt install -y nvhpc-22-2-cuda-multi

# add nvhpc to PATH
ENV       PATH "$PATH:/opt/nvidia/hpc_sdk/Linux_x86_64/22.2/compilers/bin"
ENV       LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/opt/nvidia/hpc_sdk/Linux_x86_64/22.2/compilers/lib"
