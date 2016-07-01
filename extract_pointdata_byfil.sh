#!/bin/bash

## Extracts point data from the WFDEI SWdown NetCDF dataset and 
## writes time series into a text file.
## Arguments:
## 1. path to file
## 2. variable name
## 3. longitude dimension name
## 4. latitude dimension name
## 5. longitude value of site (point)
## 6. latitude value of site (point)

# echo "extracting from $1, output is in out.txt"
ncks -s '%13.9f\n' -C -H -d $3,$5 -d $4,$6 -v $2 $1 >$HOME/tmp/tmp.txt
sed '/^$/d' $HOME/tmp/tmp.txt >$HOME/tmp/out.txt
