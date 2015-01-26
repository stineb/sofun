#!/bin/bash

## This script downloads a list of files from a remote host (dropbox in this case).
## The file list is constructed using 'fprint' to concatenate strings with a 
## numeric counter variable representing months in this case.
## b.stocker@imperial.ac.uk

## Create a list of filenames using 'fprint'
rm filelist.txt
touch filelist.txt
let i=1
while [[ $i -lt 13 ]]; do
  printf "https://dl.dropboxusercontent.com/u/26553053/Rainf_daily_WFDEI_CRU_2000%02d.nc.gz \n"  $i >> filelist.txt
  let i=$i+1
done

## Download files in this list using 'wget'
wget -i filelist.txt

