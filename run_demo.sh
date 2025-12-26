#!/bin/bash

if [ ! -f ./astronaut.bin ]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

make

#visualize intra-coding: prediction, transformation, quantization
./comp_demo astronaut.bin 512 512 4 8 3
python3 visualize.py 512 512 demo_stepsize8.svg

./comp_demo astronaut.bin 512 512 4 32 3
python3 visualize.py 512 512 demo_stepsize32.svg

./comp_demo radial1024.bin 1024 1024 4 2 1
python3 visualize.py 1024 1024 demo_partdepth1.svg

./comp_demo radial1024.bin 1024 1024 4 2 3
python3 visualize.py 1024 1024 demo_partdepth3.svg

#visualize RD-points for a specific parameter space (predMode x quantStepSize x quadSplitDepth)
./eval_RD.py RD_demo.svg