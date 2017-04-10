#!/bin/bash

list=`cat ./sitelist166_fluxnet2015.txt`

# list="DE-Seh"

for isit in ${list}
do
  for iyr in {2000..2016}
  do 
    echo $iyr
    fil=`ls /alphadata01/bstocker/data/modis_fpar_fluxnet_cutouts_FROMTREVOR/downloadedData/${isit}/*Start${iyr}*asc`
    dir=/alphadata01/bstocker/data/modis_fpar_fluxnet_cutouts_FROMTREVOR/downloadedData/${isit}/${iyr}
    mkdir ${dir}
    mv ${fil} ${dir}
  done
done
