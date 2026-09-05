#include <stdio.h>
#include <stdbool.h>

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


typedef unsigned char uchar;
typedef unsigned long ulong;

//calculate the Lagrange multiplier from quant size
double calcLambda(int stepSize);
//few helper functions
int calcBitdepth(int* src, int n);
ulong dist_ssd(int* src, int* ref, int width, int height, int stride);
void clipLR(int* src, int width, int height, int stride, int min, int max);
void array_to_file(const char* fname, const void* data, int typeSize, int len);
void array_from_file(const char* fname, void* data, int typeSize, int len);
//partitioning
int partition(int idx, int width, int height, int gridWidth, int stride, bool* hasLeft, bool* hasTop);
//add or subtract prediction
void predict(int predMode, int* reco, int* pred, int* resi, int width, int height, int stride, bool hasLeft, bool hasTop);
//transform predicted block and reconstruct original coefficients
void transform(int* src, int width, int height, int stride);
void inv_transform(int* src, int width, int height, int stride);
//quantize and dequantize
void quantize(int* src, int width, int height, int stride, int quantsize, int cutMode);
void dequantize(int* src, int width, int height, int stride, int quantsize);
//encoding and decoding
unsigned encode_coding_params(int width, int height, int stepSize, uchar* bitsOut, unsigned bitPos);
unsigned decode_coding_params(uchar* bitsIn, unsigned bitPos, int* width, int* height, int* stepSize);
unsigned encode_unit(int* quant, int stride, int partDepth, int* predModes, int* cutModes, uchar* binStream, unsigned bitPos);
unsigned decode_unit(uchar* binStream, unsigned bitPos, int* quant, int stride, int* partDepth, int* predModes, int* cutModes);
//compress and reconstruct input image
unsigned rd_est_bits(int* x, int width, int height, int stride, int cutMode);
void comp_reco_subblk(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int stride, int bitdepth, int quantSize, int predMode, int cutMode, bool topMargin, bool leftMargin);
void comp_reco_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int bitdepth, int quantSize, int partDepth, int* predModes, int* cutModes, bool topMargin, bool leftMargin);
double rd_search_subblk(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int stride, int bitdepth, int quantSize, bool topMargin, bool leftMargin, double lambda, int partDepth, int* bestPred, int* bestCut);
void rd_search_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int bitdepth, int quantSize, bool topMargin, bool leftMargin, double lambda, int* bestDepth, int* bestPreds, int* bestCuts);
void compress_image(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int stepSize, double lambda, uchar* binStream, unsigned* bitPos);
void reconstruct_image(int* quant, int* reco, int* width, int* height, int* stepSize, uchar* binStream, unsigned* bitPos);