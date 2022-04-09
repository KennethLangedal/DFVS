#!/bin/bash

tried_total=0
solved_total=0
tmax=20s

for file in $(find data/exact/ -type f | sort); do
    echo -n $file $(head -n 1 $file | awk '{ print $1 }')" "
    timeout $tmax ./build/dfvs < $file > data/solution/tmp
    if [ $? -eq 0 ] 
    then
       solved_total=$(( $solved_total + 1 ))
       echo -n "Opt found "
       echo -n $(wc -l < data/solution/tmp) " "
       ./verifier/verifier $file data/solution/tmp
    else
       echo "TLE"
    fi
    tried_total=$(( $tried_total + 1 ))
done

echo "----"
echo $solved_total / $tried_total