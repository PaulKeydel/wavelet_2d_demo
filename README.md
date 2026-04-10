# wavelet_2d_demo
This project was started to illustrate the idea behind lossy image compression and it shows how so-called Wavelets can be used in image compression.
The demonstration is based on a biorthogonal 9/7 wavelet transform that is implemented by using both a classical convolution method and a lifting scheme. Please note that the width and height of the signal must be a power of 2. The first half part of the output signal then contains the approximation coefficients and the second half part contains the wavelet coefficients (the detail coefficients). As lossy image compression formats do, these coefficients will be quantized to produce many zeros in the detail coefficients.

**The principle behind image compression**

First, the original image is splitted into blocks of same size (biggest block size is 256x256 pixel but blocks can be smaller). The sub-block structure makes it possible to apply the compression block by block:
1. Predict the current sub-block by using the blocks before. Prediction can be done horizontally, vertically or diagonally. The first block in the upper left corner is always predicted with a constant value.
![Prediction demo](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_prediction.svg)
2. Transform the residuum (i.e. original minus prediction) using the CDF9/7 wavelet transform. The resulting coefficients will then be quantized with an arbitrary stepsize in order to achieve gainful bitrates.
![Trafo demo](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_transform.svg)
3. Encode the qunatized coefficients using Huffman coding. Note that modern compression algorithms use the more efficient arithmetic coding instead of Huffman.
![Encoding demo](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_encoding.svg)
4. Use the obtaining distortion and the length of resulting bitstream to calculate the costs (known as RD-optimization). On the basis of RD-costs one can determine the best compression configuration.
![RD demo](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_RD_costs.svg)


***Example: The effect of partitioning***
![Compression demo blk-128](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_camera_pred4_qs8_depth2.svg)

![Compression demo blk-32](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_camera_pred4_qs8_depth5.svg)

***Example: The effect of quantization***
![Compression demo qp-3](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_astronaut_pred4_qs8_depth3.svg)

![Compression demo qp-5](https://github.com/PaulKeydel/wavelet_2d_demo/blob/main/visuals/demo_astronaut_pred4_qs32_depth3.svg)

**How to use this project?**

If you are working on an UNIX system, just run `./run_demo.sh`. Otherwise follow these three steps:

1. Run `python3 create_orig_img.py`. This python script will generate two different black-white images, `astronaut.bin` and `camera.bin`.

2. Build the C project using the makefile (configured for gcc, maybe you'll need to change the compiler). After building the executable will compress a given binary source file and expects 6 parameters: the source file itself, width and height of the raw image, the prediction mode, the qunatization step-width and the partitioning depth. For example, if we want to compress the astronaut image (512x512 pixel) with pred-mode 4, step-size 8 and block-size 128, we simply run `./comp_demo astronaut.bin 512 512 4 8 1`. The residual, the quantized coeffs, the reconstructed image and its encoding will be stored as binary/text files.

3. Visualize all bin files by running the python script `python3 visualize.py`. All visuals can be found in the `visuals` directory.