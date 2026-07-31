#include <stdio.h>
#include <stdbool.h>

typedef unsigned char uchar;

//config
#define IMG_BITDEPTH     8
#define MAX_BLOCK_SIZE   128
#define USE_WAVE_LIFTING 1
//encoder related config
#define ENC_MAX_DEPTH    3
#define ENC_SPECIAL_CUT  1

//prediction modes
#define PRED_CONST   0
#define PRED_HOR     1
#define PRED_VER     2
#define PRED_DIAG    3
#define NUM_PREDS    4

//coef cutting modes
#define CUT_OFF      0
#define CUT_HH       1
#define CUT_LH_HL_HH 2
#define NUM_CUTTINGS 3

//2-dim trafo subband types
#define SUBBAND_LL   0
#define SUBBAND_LH   1
#define SUBBAND_HL   2
#define SUBBAND_HH   3

//macros for bit manipulation on unsigned char array
#define setBit(arr,k,val) ( arr[(k/8)] = (val == 1) ? arr[(k/8)] | (1 << (7 - (k%8))) : arr[(k/8)] & ~(1 << (7 - (k%8))) )
#define checkBit(arr,k)   ( (arr[(k/8)] >> (7 - (k%8))) & 1 )

//helper functions
int calcBitdepth(int* x, int n);
double calcLambda(int stepSize);
unsigned long mse_dist(int* src, int* reco, int width, int height, int stride);
int clipLR(int val, int min, int max);
int w2D_subband(int row, int col, int width, int height, int stride);
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
//encoding and decoding
unsigned encode_fixlen(int n, int len, uchar* bitsOut, unsigned bitPos);
unsigned decode_fixlen(int len, uchar* bitsIn, unsigned bitPos, int* n);
unsigned encode_huffman(int n, uchar* bitsOut, unsigned bitPos);
unsigned decode_huffman(uchar* bitsIn, unsigned bitPos, int* n);
unsigned encode_coding_params(int width, int height, int stepSize, uchar* bitsOut, unsigned bitPos);
unsigned decode_coding_params(uchar* bitsIn, unsigned bitPos, int* width, int* height, int* stepSize);
unsigned encode_unit(int* quant, int stride, int partDepth, int* predModes, int* cutModes, uchar* binStream, unsigned bitPos);
unsigned decode_unit(uchar* binStream, unsigned bitPos, int* quant, int stride, int* partDepth, int* predModes, int* cutModes);
//compress and reconstruct input image
void cut_detail_coefs(int* src, int width, int height, int stride, int cutMode);
unsigned rd_est_bits(int* x, int width, int height, int stride, int cutMode);
void comp_subblk(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int stride, int bitdepth, int quantSize, int predMode, int cutMode, bool topMargin, bool leftMargin);
void reco_subblk(int* quant, int* reco, int width, int height, int stride, int bitdepth, int quantSize, int predMode, bool topMargin, bool leftMargin);
void rd_search_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int bitdepth, int quantSize, bool topMargin, bool leftMargin, double lambda, int* bestDepth, int* bestPreds, int* bestCuttings);
void comp_reco_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int bitdepth, int quantSize, int partDepth, int* predModes, int* cutModes, bool topMargin, bool leftMargin);
void compress_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int stepSize, bool topMargin, bool leftMargin, double lambda, uchar* binStream, unsigned* bitPos);
void reconstruct_unit(int* quant, int* reco, int stride, int stepSize, bool topMargin, bool leftMargin, uchar* binStream, unsigned* bitPos);