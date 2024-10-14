#!/bin/bash

if [ ! -f ./astronaut.bin ]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

make

./comp_demo astronaut.bin 512 512 1
python3 visualize.py 512 512 demo_stepsize1.svg

./comp_demo astronaut.bin 512 512 8
python3 visualize.py 512 512 demo_stepsize8.svg

./comp_demo astronaut.bin 512 512 32
python3 visualize.py 512 512 demo_stepsize32.svg
