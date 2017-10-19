#!/bin/bash

runname="EXAMPLE_global"
outyears=`seq 2005 2010`
let dt=1

cd output_nc

## sum up to annual values
for yr in ${outyears}
do
	cdo yearsum ${runname}.${yr}.d.gpp.nc ${runname}.${yr}.a.gpp.nc
done

# ## output is in daily average values, multiply by output periodicity
# cdo -O mulc,dt ${runname}.${yr}.a.gpp.nc ${runname}.${yr}.a.gpp.nc

cd ..