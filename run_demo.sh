#!/bin/bash

if [ ! -f ./astronaut.bin ]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

make
./comp_demo astronaut.bin 512 512
python3 visualize.py 512 512
