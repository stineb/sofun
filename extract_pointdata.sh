#!/bin/bash

## Extract values from NetCDF file for one point (gridcell) corresponding to mylon and mylat.
## This script is written for processing WATCH forcing data, there each file contains data
## for one month. Loop over monthly files and write all daily values into a single file for 
## this year (corresponding to 'myyear').
## b.stocker@imperial.ac.uk

myyear='2000'
mylon='7.73'
mylat='47.3'

rm dtemp_tmp_${myyear}.txt
touch dtemp_tmp_${myyear}.txt

rm dprec_tmp_${myyear}.txt
touch dprec_tmp_${myyear}.txt

## Temperature
for idx in `find Tair_daily_WFDEI_${myyear}*.nc`
do
  ncks -s '%13.9f\n' -C -H -d lon,${mylon},${mylon} -d lat,${mylat},${mylat} -v Tair $idx >>dtemp_tmp_${myyear}.txt
done    

## Precipitation
for idx in `find Rainf_daily_WFDEI_CRU_${myyear}*.nc`
do
  ncks -s '%13.9f\n' -C -H -d lon,${mylon},${mylon} -d lat,${mylat},${mylat} -v Rainf $idx >>dprec_tmp_${myyear}.txt
done    

## Cloud cover
ncks -s '%13.9f\n' -C -H -d lon,${mylon},${mylon} -d lat,${mylat},${mylat} -v cld cru_ts3.21.1901.2012.cld.dat.nc >mccov_tmp_${myyear}.txt

## Get rid of blank lines in dtemp and dprec files
sed '/^$/d' dtemp_tmp_${myyear}.txt > dtemp_CH-Oe1_${myyear}.txt
sed '/^$/d' dprec_tmp_${myyear}.txt > dprec_CH-Oe1_${myyear}.txt

