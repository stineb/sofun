#!/bin/sh
#PBS -l walltime=01:00:00,mem=1gb
## This tells the batch manager to limit the walltime for the job to XX hours, YY minutes and ZZ second
## and use PP gb of memory.

module load R
module load nco
module load netcdf
## This jobs requires the Intel math kernel so we must load it at run time.

R CMD BATCH $HOME/sofun/utils_sofun/prepare_input/get_climate_sofun.R 
## This tells the batch manager to execute the program lazy from the examples
## directory of the users home directory.