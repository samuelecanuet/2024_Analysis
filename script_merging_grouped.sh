#!/bin/bash

YEAR="2024"



# run_list=$(awk '
# /#32Ar/ {f=1}
# /#33Ar/ {f=2}
# /^#/ && !/#32Ar|#33Ar/ {f=0}
# f==1 || f==2 {
#     if ($0 !~ /^#/) print
# }' Grouper/Config_Files/${YEAR}/Runs_${YEAR}.txt | tr ' ' '\n')

# only if 32Ar
run_list=$(awk '
/#32Ar / {f=1}
/^#/ && !/#32Ar / {f=0}
f==1 {
    if ($0 !~ /^#/) print
}' Grouper/Config_Files/${YEAR}/Runs_${YEAR}.txt | tr ' ' '\n')


runs=($run_list)


command="hadd -T -f /mnt/hgfs/shared-2/${YEAR}_DATA/DETECTOR_DATA/GROUPED/merged_${YEAR}_grouped.root "
###### GROUPING RUN ######
for run in "${runs[@]}"
do
    command+="/mnt/hgfs/shared-2/${YEAR}_DATA/DETECTOR_DATA/GROUPED/run_${run}_*_grouped.root "
done

# Execute the hadd command
eval $command