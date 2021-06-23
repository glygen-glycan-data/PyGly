#!/bin/sh
USER="`id -u`:`id -g`"
CUR=`pwd`
CUR=`readlink -f $CUR`
DOCKER="docker run -u $USER -v $CUR:/data/:Z -v /data2/projects/GlyGen/PyGly/smw/gptwiki/:/data2/projects/GlyGen/PyGly/smw/gptwiki/:Z --rm"
READLINK=`readlink -f "$0"`
BASE=`basename "$0" .sh`
READLINKBASE=`basename "$READLINK" .sh`
docker pull openswath/openswath:latest
if [ "$BASE" == "$READLINKBASE" ]; then
  $DOCKER -i -t openswath/openswath:latest                                                                                   
else
  $DOCKER openswath/openswath:latest $BASE "$@"
fi                                                                                                                           
