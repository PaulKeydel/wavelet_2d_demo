#include <stdio.h>
#include <stdbool.h>

typedef unsigned char uchar;

//config
#define IMG_BITDEPTH     8
#define MAX_BLOCK_SIZE   128
#define USE_TAUBMANN     0
//encoder related config
#define SEARCH_SPLT_DEPTH 4
#define SEARCH_CUT_DEPTH  3

//prediction modes
#define PRED_CONST   0
#define PRED_HOR     1
#define PRED_VER     2
#define PRED_DIAG    3
#define NUM_PREDS    4

//macros for bit manipulation on unsigned char array
#define setBit(arr,k,val) ( arr[(k/8)] = (val == 1) ? arr[(k/8)] | (1 << (7 - (k%8))) : arr[(k/8)] & ~(1 << (7 - (k%8))) )
#define checkBit(arr,k)   ( (arr[(k/8)] >> (7 - (k%8))) & 1 )

//helper functions
int calcBitdepth(int* x, int n);
double calcLambda(int stepSize);
int clipLR(int val, int min, int max);
void blkcpy(int* src, int* dst, int width, int height, int stride);
void array_to_file(const char* fname, const void* data, int typeSize, int len);
void array_from_file(const char* fname, void* data, int typeSize, int len);
//add or subtract prediction
void predict(int predMode, int* reco, int* dst, int* resi, int width, int height, int stride, bool hasLeft, bool hasTop, int bitdepth);
//transform predicted block and reconstruct original coefficients
void transform(int* src, int width, int height, int stride, int bitdepthIn);
void inv_transform(int* src, int width, int height, int stride, int bitdepthOut);
//quantize and dequantize
void quantize(int* src, int width, int height, int stride, int quantsize);
void dequantize(int* src, int width, int height, int stride, int quantsize);
//calculate rate and distortion (MSE)
unsigned long rd_unit_bits(int* x, int stride, int partDepth, int cutCoefs);
unsigned long mse_dist(int* src, int* reco, int width, int height, int stride);
//encoding and decoding
unsigned encode_fixlen(int n, int len, uchar* bitsOut, int bitPos);
unsigned decode_fixlen(int len, uchar* bitsIn, int bitPos, int* n);
unsigned encode_huffman(int n, uchar* bitsOut, int bitPos);
unsigned decode_huffman(uchar* bitsIn, int bitPos, int* n);
unsigned encode_coding_params(int width, int height, int stepSize, uchar* bitsOut, int bitPos);
unsigned decode_coding_params(uchar* bitsIn, int bitPos, int* width, int* height, int* stepSize);
unsigned encode_unit(int* quant, int stride, int partDepth, int predMode, int cutCoefs, uchar* binStream, int bitPos);
unsigned decode_unit(uchar* binStream, int bitPos, int* quant, int stride, int* partDepth, int* predMode, int* cutCoefs);
//compress and reconstruct input image
void cut_detail_coefs(int* src, int width, int height, int stride, int partDepth, int cutCoefs);
void comp_reco_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int bitdepth, int stride, int stepSize, int partDepth, int predMode, int cutCoefs, bool topMargin, bool leftMargin);
void compress_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int stepSize, int partDepth, int predMode, int cutCoefs, bool topMargin, bool leftMargin, unsigned long* rdBits, unsigned long* rdDist);
void reconstruct_unit(int* quant, int* reco, int stride, int stepSize, int partDepth, int predMode, int cutCoefs, bool topMargin, bool leftMargin);