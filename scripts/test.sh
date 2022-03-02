#!/bin/bash

left_total=0
tried_total=0
solved_total=0

for file in $(find data/exact/ -type f | sort); do
    echo -n $file $(head -n 1 $file | awk '{ print $1 }')" "
    n_left=$(./build/dfvs < $file)
    if [ "$n_left" = "Too Large" ]
    then
        echo $n_left
        continue
    else
        tried_total=$(( $tried_total + 1 ))
    fi
    left_total=$(( $left_total + $n_left ))
    echo -n $n_left" "
    if [ $n_left -eq 0 ] 
    then
        solved_total=$(( $solved_total + 1 ))
        echo -n "Opt found "
        ./verifier/verifier $file data/solution/tmp
    else
        echo
    fi
done

echo "----"
echo $left_total
echo $solved_total / $tried_total