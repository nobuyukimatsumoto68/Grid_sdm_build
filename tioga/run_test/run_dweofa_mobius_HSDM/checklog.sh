#!/bin/bash
# parse texts from log to test the run


log=$1

grep Plaquette ${log} | awk -F: '{print $4, $5}' | uniq
grep Polyakov ${log} | awk -F: '{print $4, $5}'
grep NERSC ${log} | awk -F: '{print $4}'
