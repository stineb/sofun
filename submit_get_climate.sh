#!/bin/sh
#PBS -l walltime=99:00:00,mem=1gb
## This tells the batch manager to limit the walltime for the job to XX hours, YY minutes and ZZ second
## and use PP gb of memory.

module load R
module load nco
module load netcdf

R CMD BATCH $HOME/sofun/getin/get_climate.R $HOME/sofun/getin/get_climate.Rout
