#!/bin/bash

code_length_dist() {
  for i in {1..5}
  do
    num="$(grep -E "):\s\b\w{${i}}\b" bitstream.txt | wc -l)"
    echo "Number of Huffman-coded coefficients that have length "$i":"
    echo "$num"
  done
}

if ! [[ -f ./astronaut.bin && -f ./camera.bin ]]; then
    echo "Creating demo images..."
    python3 create_orig_img.py
fi

make

read -p "Start Jupyter demo (j) OR generate visuals (v) OR evaluate costs (c)? [j/v/c] " -n 1 choice
echo
if [ "$choice" = "j" ]; then 
    echo "Running Jupyter..."
    jupyter-notebook visualize_nb.ipynb
elif [ "$choice" = "v" ]; then
    #visualize intra-coding (i.e. prediction, transformation, quantization) and RD-optimization
    echo "Create SVG image files..."
    python3 visualize.py
    code_length_dist
elif [ "$choice" = "c" ]; then
    #print summary of RD analysis
    echo "Display summary..."
    ./eval_RD.py
else
    echo "Invalid choice!"
fi
