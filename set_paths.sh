#!/bin/bash

# arguments:
# 1. original path
# 2. user path

perl -i -pe's/${1}/${2}/g' ./src/*f90
perl -i -pe's/${1}/${2}/g' ./params/parameters_Allocation.nml