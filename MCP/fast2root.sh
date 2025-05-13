#!/bin/bash

Runs=("001" "005" "006" "007" "008" "009" "010" "011" "027" "028" "029" "030" "031" "032" "033" "034" "035" "036" "037" "038" "039" "041" "042" "077" "079" "080" "081" "082")

for run in "${Runs[@]}"
do
    files=$(find /run/media/local1/Disque_Dur/2025_DATA/MCP_DATA/ -type f \( -name "run_${run}_*.fast" -o -name "run_${run}_*_0001.fast" \))

    for file in $files; do
        cd /home/local1/Documents/Fast2Root
        Fast2Root "$file" -o /mnt/hgfs/shared-2/2025_DATA/MCP_DATA/
        cd -
        break
    done
done
