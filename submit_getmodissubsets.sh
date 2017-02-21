#!/bin/sh
#PBS -l walltime=01:00:00,mem=1gb
## This tells the batch manager to limit the walltime for the job to XX hours, YY minutes and ZZ second
## and use PP gb of memory.

module load R
module unload liblzma

## This jobs requires the Intel math kernel so we must load it at run time.

R CMD BATCH --no-save --no-restore $HOME/sofun/getin/get_modissubset_tseries.R

## This tells the batch manager to execute the program lazy from the examples
## directory of the users home directory.
