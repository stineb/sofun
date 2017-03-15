#!/bin/bash
for i in {1..10} 
do 
  for i in {1..10} 
  do 
    for i in {1..10} 
    do 
      nohup nice R CMD BATCH get_modissubsets.R get_modissubsets.Rout
      echo "sleeping 1 min. zzz"
      sleep 60
    done
    echo "sleeping 5 min. zzzzzzzz"
    sleep 300
  done
  echo "sleeping 30 min. zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
  sleep 1800
done
