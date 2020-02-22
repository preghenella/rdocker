#! /usr/bin/env bash

docker run --rm -it \
       --mount type=bind,source=`pwd`,target=/home/fastuser \
       --mount type=bind,source=/tmp/.X11-unix,target=/tmp/.X11-unix \
       -e DISPLAY=$DISPLAY \
       preghenella/fastsim $@
