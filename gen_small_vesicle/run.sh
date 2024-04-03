#!/bin/bash

##path="/media/weria/Disk_1/work/DTS/Vesicle_Caimirforce_Ranjit/src/"
path="/Users/weria/Documents/DTSv2.0/src/dts_src/"

if [ $1 -gt 10 ]; then
    echo "The seed number is good "
   "$path"DTS -in input.dts -top out.tsi  -seed $1


else
    echo "seed is bad or has not been provided "
fi


#"$path"DTS -in input.dts -top top.top  -seed 98231

