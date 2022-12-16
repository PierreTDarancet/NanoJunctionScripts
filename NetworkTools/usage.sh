#!/bin/bash
NPROC=`grep processor /proc/cpuinfo | wc -l`
USAGE=`ps -e -o pcpu --sort pcpu | awk  'NR!=1' | paste -sd+ | bc`
echo "$USAGE/$NPROC" | bc

