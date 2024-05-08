# wavelet_2d_demo
Two implementations of the wavelet transform (CDF 9/7) and a small demonstration of how wavelets can be used in image compression.
The demonstration is based on a biorthogonal 9/7 wavelet transform that is implemented by using both a classical convolution method and a lifting scheme. Please note that the width and height of the signal must be a power of 2. The first half part of the output signal then contains the approximation coefficients and the second half part contains the wavelet coefficients (the detail coefficients). As lossy image compression formats do, these coefficients will be quantized to produce many zeros in the detail coefficients.

**How to use this project?**

1. Run `python3 create_orig_img.py`. This python script will generate two different black-white images, `astronaut.bin` and `art8x8.bin`.

2. Compile and run the C file. For example, use astronaut.bin (512x512 pixel) as input image: `./w97 astronaut.bin 512 512`. The executable will perform the wavelet transform and the quantization. The wavelet coeffs, the quantized coeffs and the reconstructed image are stored in new bin files.

3. Visualize all bin files by running the python script `python3 visualize.py 512 512`. You can see the original (uncompressed) image, both versions of the wavelet coefficients and the effect of quantization within the reconstructed image.