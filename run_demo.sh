#!/bin/bash

if ! [[ -f ./astronaut.bin && -f ./camera.bin ]]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

make

read -p "Only start the Jupyter demo notebook? [y/n] " -n 1 choice
echo
if [ "$choice" = "y" ]; then 
    echo "Running Jupyter..."
    jupyter-notebook visualize_nb.ipynb
elif [ "$choice" = "n" ]; then
    #visualize intra-coding (i.e. prediction, transformation, quantization) and RD-optimization
    echo "Create SVG image files..."
    python3 visualize.py
else
    echo "Invalid choice!"
fi
