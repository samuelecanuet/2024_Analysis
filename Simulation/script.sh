#!/bin/bash
min=-3.0
max=3.0

for i in $(seq 2.5 0.5 $max); do
    for j in $(seq -1.5 0.5 $max); do
       for t in $(seq $min 0.5 $max); do
        	echo "Position: ($i, $j, $t)"

             # extract fractional parts
            frac_i=${i#*.}
            frac_j=${j#*.}
            frac_t=${t#*.}

            # if *all three* are integers (i.e. .0), skip this iteration
            # if [[ $frac_i -eq 0 && $frac_j -eq 0 && $frac_t -eq 0 ]]; then
            #     continue
            # fi  

            if [ $(echo "$i >= -2.0 && $i <= 2.0" | bc) -eq 1 ] && \
               [ $(echo "$j >= -2.0 && $j <= 2.0" | bc) -eq 1 ] && \
               [ $(echo "$t >= -2.0 && $t <= 2.0" | bc) -eq 1 ]; then
                continue
            fi
          
            rsync -av --progress --ignore-existing lecanuet@borlin333:/home/lecanuet/../../../data333/lecanuet/Result/"*_x*y${i}_z${j}_theta${t}*_CVP1.root" /run/media/local1/DATANEX/Samuel-G4/new/new/
	    done
    done
done
