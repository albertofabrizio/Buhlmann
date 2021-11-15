#!/bin/bash

usage() { echo -e "Usage: $0 -d 40 -p 210 -s 24 -m 40 -o 0.21 -e 0.0 -g 0.75 \n -d Depth [m] \n -p Twinset total pressure [bar] \n -s Twinset total volume [l] \n -m Minimum allowed twinset pressure [bar] \n -o Percentage of Oxygen in Bottom Gas -e Percentage of Helium in Bottom Gas" 1>&2; exit 1; }

while getopts d:p:s:m:o:e:g:h: flag
do
    case "${flag}" in
        d) depth=${OPTARG};;
        p) pressure=${OPTARG};;
        s) size=${OPTARG};;
        m) min=${OPTARG};;
        o) oxo=${OPTARG};;
        e) helium=${OPTARG};;
        g) gl=${OPTARG};; 
        *) usage ;;
    esac
done

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    exit 1
fi

if [ -z $depth ]; then
    echo "No target depth provided"
    exit 1
fi
if [ -z $pressure ]; then
    echo "No tank pressure provided"
    exit 1
fi
if [ -z $size ]; then
    echo "No tank size provided"
    exit 1
fi
if [ -z $min ]; then
    echo "No minimum safety margin provided"
    exit 1
fi

if [ -z $oxo ]; then
    echo "No oxygen provided, using default 0.21"
    oxo="0.21"
fi

if [ -z $helium ]; then
    echo "No helium provided, using default 0.0"
    helium="0.0"
fi

if [ -z $gl ]; then
    echo "No gradient factor low provided, using default 0.75"
    gl="0.75"
fi






liters=$[ $pressure * $size  ]
liters_safe=$[ $liters - $min * $size ]
echo "Total volume of air in twinset: $liters l"

echo "Computing maximum time at depth $depth m, assuming loss of decompression stage and minimum allowed tank pressure $min bar. Bottom gas $oxo/$helium."
i=1
while :
do
    python ../deco_plan.py --depth $depth --time $i  --fo2 $oxo --fhe $helium --glow $gl > tmp
    tot_gas=`tail -1 tmp | awk '{print int($3)}'`
    
    if [[ $tot_gas -gt $liters_safe ]]; then
        break
    fi

    answer="Gas consumption after $i minutes at depth : $tot_gas [l]"
    remain_gas=$[ $liters - $tot_gas ]
    remain_pressure=$(echo "scale=2; $remain_gas/$size" | bc)
    
    i=$[$i+1]

done

echo $answer

echo "Estimated remaining gas upon surfacing $remain_gas l or $remain_pressure bar."

if [ ! -z tmp ]; then
    rm tmp
fi
