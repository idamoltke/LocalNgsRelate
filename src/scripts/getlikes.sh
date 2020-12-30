#!/bin/bash

FILENAME=$1
grep Obtained ${1} | awk -F'[ =:]' '{print $6 "\t" $17}'  | sort -k2 
