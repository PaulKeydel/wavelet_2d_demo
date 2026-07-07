#include <stdio.h>

//config
#define USE_TAUBMANN    0
#define MAX_BLOCK_SIZE  128
#define LL_FIXED_LENGTH 0
#define TRANSFORM_SKIP  0

//helper functions
int calcBitdepth(int* x, int n);
double calcLambda(int stepSize);
int clipLR(int val, int min, int max);
void blkcpy(int* src, int* dst, int width, int height, int stride);
//add or subtract prediction
void predict(int predMode, int* reco, int* dst, int* resi, int width, int height, int stride, bool hasLeft, bool hasTop, int bitdepth);
//transform predicted block and reconstruct original coefficients
void transform(int* src, int width, int height, int stride, int bitdepthIn);
void inv_transform(int* src, int width, int height, int stride, int bitdepthOut);
//quantize and dequantize
void quantize(int* src, int width, int height, int stride, int quantsize);
void dequantize(int* src, int width, int height, int stride, int quantsize);
//calculate rate and distortion (MSE)
unsigned long coded_bits(int* x, int width, int height, int stride, int bitdepth, int Qp);
unsigned long mse_dist(int* src, int* reco, int width, int height, int stride);
//Huffman-coding of quantized coefficients, write bitstream to file
void encode_fixlen(int n, int len, char* bitsOut);
void encode_huffman(int n, char* bitsOut);
void encode_unit(int* quant, int bitdepth, int stride, int stepSize, int partDepth, int predMode, FILE* fptrEncTxt, FILE* fptrEncBin);
//apply all and compress input image
void compress_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int bitdepth, int stride, int stepSize, int partDepth, int predMode, bool topMargin, bool leftMargin, unsigned long* bits, unsigned long* dist);