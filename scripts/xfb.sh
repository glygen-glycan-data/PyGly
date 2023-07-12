#!/bin/sh
Xvfb :1 &
XSCR=$!
export DISPLAY=localhost:1.0
