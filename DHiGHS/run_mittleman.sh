#!/bin/bash


i=0

for f in ../mittleman_symmetric_mps/*; do 
    if [ $i -lt 135 ]
    then
        ((i++))
        continue
    fi
    ./build/bin/highs "$f" 1
    ./build/bin/highs "$f" 0
done

# while [ $i -le 2 ]
# do
#   echo Number: $i
#   ((i++))
# done
