#! /usr/bin/env bash

docker run --rm -it \
       --mount type=bind,source=`pwd`,target=/home/fastuser \
       -e DISPLAY=host.docker.internal:0 \
       preghenella/fastsim $@
