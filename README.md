# wavelet_2d_demo
Two implementations of the wavelet transform (CDF 9/7) and a small demonstration of how wavelets can be used in image compression.
The demonstration is based on a biorthogonal 9/7 wavelet transform that is implemented by using both a classical convolution method and a lifting scheme. Please note that the width and height of the signal must be a power of 2. The first half part of the output signal then contains the approximation coefficients and the second half part contains the wavelet coefficients (the detail coefficients). As lossy image compression formats do, these coefficients will be quantized to produce many zeros in the detail coefficients.

**The principle behind image compression**

First, the original image is splitted into blocks of same size (biggest block size is 256x256 pixel but blocks can be smaller). The sub-block structure makes it possible to apply the compression block by block:
1. Predict the current sub-block by using the blocks before. The first block is predicted with a constant value.
2. Calculate the block residuum, i.e. original minus prediction.
3. Transform the residuum using the CDF9/7 trafo.
4. Quantize the resulting coefficients with a fixed stepsize. A lot of zeros will appear (hopefully).


***Example: The effect of partitioning***
![Compression demo blk-128](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/demo_camera_pred4_qs8_depth1.svg)

![Compression demo blk-32](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/demo_camera_pred4_qs8_depth3.svg)

***Example: The effect of quantization***
![Compression demo qp-3](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/demo_astronaut_pred4_qs8_depth3.svg)

![Compression demo qp-5](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/demo_astronaut_pred4_qs32_depth3.svg)

**How to use this project?**

If you are working on an UNIX system, just run `./run_demo.sh`. Otherwise follow these three steps:

1. Run `python3 create_orig_img.py`. This python script will generate two different black-white images, `astronaut.bin` and `radial1024.bin`.

2. Build the C project using the makefile (configured for gcc, maybe you'll need to change the compiler). After building the executable will compress a given binary source file and expects 5 parameters: the source file itself, width and height of the raw image, qunatization parameter QP and the partitioning depth. For example, if we want to compress the astronaut image (512x512 pixel) with step-size 8 (QP=3) and block-size 128: `./comp_demo astronaut.bin 512 512 8 1`. The executable will perform the prediction, the transform and the quantization. The residual, the quantized coeffs and the reconstructed image are stored in new bin files.

3. Visualize all bin files by running the python script `python3 visualize.py 512 512`. You can see the same figure as above.