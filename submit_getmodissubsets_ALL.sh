#!/bin/bash

for idx in `cat test.txt`
#for idx in `cat sitelist_doprofile_fluxnet2015.txt`
do 
    cp submit_getmodissubsets_XXXXXX.sh submit_getmodissubsets_${idx}.sh
    sed -i -- s/XXXXX/${idx}/g submit_getmodissubsets_${idx}.sh
  	echo "submitting getmodissubset script for site ${idx}"
  	qsub submit_getmodissubsets_${idx}.sh	
done
