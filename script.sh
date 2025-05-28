#!/bin/bash

YEAR="2024"

run_list=$(awk '
/#32Ar/ {f=1}
/#33Ar/ {f=2}
/^#/ && !/#32Ar|#33Ar/ {f=0}
f==1 || f==2 {
    if ($0 !~ /^#/) print
}' Grouper/Config_Files/${YEAR}/Runs_${YEAR}.txt | tr ' ' '\n')

runs=($run_list)

for run in "${runs[@]}"
do
    cd Grouper
    Grouper $run 
    cd -
done

cd Grouper
Merger
