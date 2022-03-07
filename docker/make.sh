#!/bin/bash

YYYYMMDD=$(date '+%Y%m%d')
V1=11.2
V2=22.2

#####
#
# create the image for the toolchains
# docker build -t tinkertools/tinker9-devel:cuda${V1}-nvhpc${V2} -f 1-devenv.dockerfile .
#
# create the image for the runtime
# docker build -t tinkertools/tinker9-runtime:cuda${V1}-nvhpc${V2} -f 2-runenv.dockerfile .
#
#####

# push the runtime image to Docker Hub
 echo Dummy push runtime
# sleep 5s

# compile tinker9 in a temporary image
docker build -t temp_build_tinker9 -f 3-compile.dockerfile .
sleep 5s

# create the image for tinker9
docker build -t tinkertools/tinker9:$YYYYMMDD -f 4-install.dockerfile .

# push the image to Docker Hub
echo Dummy push tinker9
# docker image push -a tinkertools/tinker9
sleep 5s

# delete the temporary image
docker rmi temp_build_tinker9


################################################################################


# clean all of the old containers
# docker container prune

# create an alias for an image
# docker image tag IMAGE ALIAS

# login
# docker login
