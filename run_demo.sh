#!/bin/bash

if ! [[ -f ./astronaut.bin && -f ./camera.bin ]]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

cd comp_demo
#make clean
make all
cd ..

echo "What do you want to do?"
echo "  (j)         Start Jupyter demo"
echo "  (v)         Generate figure files"
echo "  (c)         Evaluation of R-D-costs"
echo "  (b)         Validation test for bitstream"
echo "  (any other) Just do make"
read -p "Key? " -n 1 choice
echo
if [ "$choice" = "j" ]; then 
    echo "Running Jupyter..."
    open /Applications/Firefox.app/Contents/MacOS/firefox
    jupyter-notebook visualize_nb.ipynb
elif [ "$choice" = "v" ]; then
    #visualize intra-coding (i.e. prediction, transformation, quantization) and RD-optimization
    echo "Create SVG image files..."
    python3 visualize.py
elif [ "$choice" = "c" ]; then
    #print summary of RD analysis
    echo "Display summary..."
    ./eval_RD.py
elif [ "$choice" = "b" ]; then
    #decode and validate bitstream
    echo "Check if decoder matches encoder..."
    python3 decode.py
else
    echo "Invalid choice!"
fi
