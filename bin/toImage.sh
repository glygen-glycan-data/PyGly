#!/bin/sh
THEDIR=`dirname $0`
java -cp $THEDIR/pygly/GlycoCT2ImageBundle.jar GlycoCT2Image out $1 -
