### create the image for the runtime

### e.g.
###    docker build -t tinkertools/tinker9-runtime:cuda11.2-nvhpc22.2 -f 2-runenv.dockerfile .

# load the local image for the toolchains
FROM      tinkertools/tinker9-devel:cuda11.2-nvhpc22.2 as tool_image

# load the image for the runtime and install packages
FROM      nvidia/cuda:11.2.2-runtime-ubuntu20.04 as runtime_image
# git and cmake are no longer necessary
RUN       export DEBIAN_FRONTEND=noninteractive
RUN       apt update -y && apt install -y ca-certificates libssl-dev gfortran

ARG       LibPath=/opt/nvidia/hpc_sdk/Linux_x86_64/22.2/compilers/lib
# COPY    --from=SRC_IMAGE  SRC          DEST
COPY      --from=tool_image "${LibPath}" "${LibPath}"
ENV       LD_LIBRARY_PATH "$LD_LIBRARY_PATH:${LibPath}"
