  #!/bin/bash

## Extract values from NetCDF file for one point (gridcell) corresponding to mylon and mylat.
## This script is written for processing WATCH forcing data, there each file contains data
## for one month. Loop over monthly files and write all daily values into a single file for 
## this year (corresponding to 'myyear').
## b.stocker@imperial.ac.uk

sitename='CH-Oe1'
myyear=2002
mylon=7.73
mylat=47.3

dotemp=true
doprec=true
doccov=true
dovap=true

datadir="/alphadata01/bstocker/data"

## Temperature
if $dotemp
  then  
  rm dtemp_tmp_${myyear}.txt
  touch dtemp_tmp_${myyear}.txt

  for idx in `ls -1 ${datadir}/watch_wfdei/Tair_daily_WFDEI_${myyear}*.nc`
  do
    echo "extracting from $idx"
    ncks -s '%13.9f\n' -C -H -d lon,${mylon},${mylon} -d lat,${mylat},${mylat} -v Tair $idx >>dtemp_tmp_${myyear}.txt
  done  

  ## Get rid of blank lines in dtemp and dprec files
  sed '/^$/d' dtemp_tmp_${myyear}.txt > dtemp_tmp2_${myyear}.txt

  ## K -> deg C
  awk '{ printf "%13.9f\n", $1 - 273.15 }' dtemp_tmp2_${myyear}.txt > dtemp_${sitename}_${myyear}.txt

  echo "point data written to dtemp_${sitename}_${myyear}.txt"

fi


## Precipitation
if $doprec
  then
  rm dprec_tmp_${myyear}.txt
  touch dprec_tmp_${myyear}.txt

  for idx in `ls -1 ${datadir}/watch_wfdei/Rainf_daily_WFDEI_CRU_${myyear}*.nc`
  do
    echo "extracting from $idx"
    ncks -s '%13.9f\n' -C -H -d lon,${mylon},${mylon} -d lat,${mylat},${mylat} -v Rainf $idx >>dprec_tmp_${myyear}.txt
  done  

  ## Get rid of blank lines in dtemp and dprec filesxt > mccov_tmp2_${myyear}.txt
  sed '/^$/d' dprec_tmp_${myyear}.txt > dprec_tmp2_${myyear}.txt

  ## kg/m2/s -> mm
  awk '{ printf "%13.9f\n", $1 *24.0*60.0*60.0 }' dprec_tmp2_${myyear}.txt > dprec_${sitename}_${myyear}.txt

  echo "point data written to dprec_${sitename}_${myyear}.txt"

fi


## Cloud cover
if $doccov
  then
  let istart=(${myyear}-1901)*12
  let iend=(${myyear}-1901)*12+11

  echo "extracting from ${datadir}/cru/cru_ts_3.21/cru_ts3.21.1901.2012.cld.dat.nc"
  ncks -s '%13.9f\n' -C -H -d lon,${mylon},${mylon} -d lat,${mylat},${mylat} -d time,1188,1199 -v cld "${datadir}/cru/cru_ts_3.21/cru_ts3.21.1901.2012.cld.dat.nc" >mccov_tmp_${myyear}.txt

  ## Get rid of blank lines in dtemp and dprec filesxt > mccov_tmp2_${myyear}.txt
  sed '/^$/d' mccov_tmp_${myyear}.txt > mccov_tmp2_${myyear}.txt

  ## % cloud covered -> sunshine fraction {0,1}
  awk '{printf "%13.9f\n",(100.0-$1)/100}' mccov_tmp2_${myyear}.txt > mfsun_${sitename}_${myyear}.txt

  echo "point data written to mfsun_${sitename}_${myyear}.txt"

fi


## Vapour pressure
if $dovap
  then
  let istart=108
  let iend=119

  echo "extracting from ${datadir}/cru/cru_ts_3.22/cru_ts3.22.1991.2000.vap.dat.nc"
  ncks -s '%13.9f\n' -C -H -d lon,${mylon},${mylon} -d lat,${mylat},${mylat} -d time,${istart},${iend} -v vap "${datadir}/cru/cru_ts_3.22/cru_ts3.22.1991.2000.vap.dat.nc" >mvapr_tmp_${myyear}.txt

  # Get rid of blank lines in dtemp and dprec files
  sed '/^$/d' mvapr_tmp_${myyear}.txt > mvapr_tmp2_${myyear}.txt

  ## Vapour pressure (hPa)
  awk '{printf "%13.9f\n",$1}' mvapr_tmp2_${myyear}.txt > mvapr_${sitename}_${myyear}.txt

  echo "point data written to mvapr_${sitename}_${myyear}.txt"

fi

## delete temporary files
rm *tmp*
