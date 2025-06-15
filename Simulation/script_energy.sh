#!/bin/bash
min=0.1
max=7.0

rsync -av --progress --ignore-existing lecanuet@borlin333:/home/lecanuet/../../../data333/lecanuet/Result/"*MeV.root" /run/media/local1/DATANEX/Samuel-G4/06-03/
# for i in $(seq $min 0.1 $max); do
#     Reader2 $i
# done
