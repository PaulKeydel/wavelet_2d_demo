# wavelet_2d_demo
Two implementations of the wavelet transform (CDF 9/7) and a small demonstration of how wavelets can be used in image compression.
The demonstration is based on a biorthogonal 9/7 wavelet transform that is implemented by using both a classical convolution method and a lifting scheme. Please note that the width and height of the signal must be a power of 2. The first half part of the output signal then contains the approximation coefficients and the second half part contains the wavelet coefficients (the detail coefficients). As lossy image compression formats do, these coefficients will be quantized to produce many zeros in the detail coefficients.

**The principle behind image compression**
First, the original image is splitted into blocks (in our case a single quad split which yields 4 blocks of same size). The sub-block structure makes it possible to apply the compression block by block:
1. Predict the current sub-block by using the blocks before. The first block is predicted with a constant value.
2. Calculate the block residuum, i.e. original minus prediction.
3. Transform the residuum using the CDF9/7 trafo.
4. Quantize the resulting coefficients with a fixed stepsize. A lot of zeros will appear (hopefully).

![Compression demo](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/Figure_1.svg)

**How to use this project?**

If you are working on an UNIX system, just run `./run_demo.sh`. Otherwise follow these three steps:

1. Run `python3 create_orig_img.py`. This python script will generate two different black-white images, `astronaut.bin` and `art8x8.bin`.

2. Compile and run the C file. For example, use astronaut.bin (512x512 pixel) as input image: `./w97 astronaut.bin 512 512`. The executable will perform the prediction, the transform and the quantization. The residual, the quantized coeffs and the reconstructed image are stored in new bin files.

3. Visualize all bin files by running the python script `python3 visualize.py 512 512`. You can see the same figure as above.