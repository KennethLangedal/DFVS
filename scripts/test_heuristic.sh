#!/bin/bash

cost=0
tmax=60s

for file in $(find data/heuristic/ -type f | sort); do
    echo -n $file $(head -n 1 $file | awk '{ print $1 }')" "
    timeout -s SIGTERM $tmax ./build/dfvs < $file > data/solution/tmp
    if [ $? -eq 0 ] 
    then
       cost=$(wc -l < data/solution/tmp)
       solved_total=$(( $solved_total + 1 ))
       echo -n $cost " "
       ./verifier/verifier $file data/solution/tmp
    else
       echo "TLE"
    fi
done