#!/bin/bash

## argument 1: SOFUN job name

bsub -N -o ~/hpc_log/$1_%J.out -e ~/hpc_log/$1_%J.err -n 1 -W 1:00 -R "rusage[mem=5000]" -J $1 "echo $1 | ./rungpmodel"
