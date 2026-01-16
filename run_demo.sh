#!/bin/bash

if ! [[ -f ./astronaut.bin && -f ./camera.bin ]]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

make

read -p "Only start the Jupyter demo notebook? [y/n] " -n 1 choice
echo
if [ "$choice" = "y" ]; then 
    jupyter-notebook visualize_nb.ipynb
elif [ "$choice" = "n" ]; then
    echo "Create SVG image files..."

    #visualize intra-coding: prediction, transformation, quantization
    ./comp_demo astronaut.bin 512 512 4 8 3
    python3 visualize.py 512 512 3 demo_astronaut_pred4_qs8_depth3.svg

    ./comp_demo astronaut.bin 512 512 4 32 3
    python3 visualize.py 512 512 3 demo_astronaut_pred4_qs32_depth3.svg

    ./comp_demo camera.bin 512 512 4 8 2
    python3 visualize.py 512 512 2 demo_camera_pred4_qs8_depth2.svg

    ./comp_demo camera.bin 512 512 4 8 5
    python3 visualize.py 512 512 5 demo_camera_pred4_qs8_depth5.svg

    #visualize RD-points for a specific parameter space (predMode x quantStepSize x quadSplitDepth)
    ./eval_RD.py RD_demo.svg
else
    echo "Invalid choice!"
fi
