#!/bin/bash

runname="s0_fapar3g_global"
outyears=`seq 1982 2011`
let dt=1

cd output_nc

# ## sum up to annual values
# for yr in ${outyears}
# do
# 	cdo yearsum ${runname}.${yr}.d.gpp.nc ${runname}.${yr}.a.gpp.nc
# done

# ## output is in daily average values, multiply by output periodicity
# cdo -O mulc,dt ${runname}.${yr}.a.gpp.nc ${runname}.${yr}.a.gpp.nc

## stack annual files along time axis
cdo mergetime ${runname}.*.a.gpp.nc ${runname}.a.gpp.nc
cdo mergetime ${runname}.*.a.aet.nc ${runname}.a.aet.nc
cdo mergetime ${runname}.*.a.pet.nc ${runname}.a.pet.nc
cdo mergetime ${runname}.*.a.alpha.nc ${runname}.a.alpha.nc

cd ..