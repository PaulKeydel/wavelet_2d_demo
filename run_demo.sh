#!/bin/bash

if [ ! -f ./astronaut.bin ]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

make

if [ "$#" -eq 2 ]; then
    stepsize="$(echo "2^$1" | bc -l)"
    partdepth="$2"

    fname="astronaut_${stepsize}_${partdepth}.svg"
    ./comp_demo astronaut.bin 512 512 $stepsize $partdepth
    python3 visualize.py 512 512 $fname

    fname="radial_${stepsize}_${partdepth}.svg"
    ./comp_demo radial1024.bin 1024 1024 $stepsize $partdepth
    python3 visualize.py 1024 1024 $fname
else
    ./comp_demo astronaut.bin 512 512 8 1
    python3 visualize.py 512 512 demo_stepsize8.svg

    ./comp_demo astronaut.bin 512 512 32 1
    python3 visualize.py 512 512 demo_stepsize32.svg

    ./comp_demo radial1024.bin 1024 1024 2 1
    python3 visualize.py 1024 1024 demo_partdepth1.svg

    ./comp_demo radial1024.bin 1024 1024 2 3
    python3 visualize.py 1024 1024 demo_partdepth3.svg
fi