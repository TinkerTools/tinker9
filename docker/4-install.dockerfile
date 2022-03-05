### create the image for tinker9

### command
###     docker build -t tinkertools/tinker9:YYYYMMDD -f 4-install.dockerfile .

# load the local tinker9 build
FROM      temp_build_tinker9 as build_image

# load the runtime image for tinker9
FROM      tinkertools/tinker9-runtime:cuda11.2-nvhpc22.2 as runtime_image

# COPY    --from=SRC_IMAGE   SRC               DEST
COPY      --from=build_image /home/tinker9/bin /home/tinker9/bin
ENV       PATH "$PATH:/home/tinker9/bin/gpu-m"
