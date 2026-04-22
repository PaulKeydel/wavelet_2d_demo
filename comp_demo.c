#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "w97.h"

#define USE_TAUBMANN    0
#define MAX_BLOCK_SIZE  256
#define LL_FIXED_LENGTH 0

int clipLR(int val, int min, int max)
{
  if (val > max)
  {
    return max;
  }
  if (val < min)
  {
    return min;
  }
  return val;
}

int calcBitdepth(int* x, int n)
{
  int max = 0;
  int min = 0;
  for (int i = 0; i < n; i++)
  {
    if (x[i] > max) max = x[i];
    if (x[i] < min) min = x[i];
  }
  return (int)floor(log2((double)(max - min))) + 1;
}

void predict(int predMode, int* reco, int* dst, int* resi, int width, int height, int stride, bool hasLeft, bool hasTop, int bitdepth)
{
  int defaultPred = 1 << (bitdepth - 1);
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int pred = 0;
      if (predMode == 0)
      {
        pred = defaultPred;
      }
      else if (predMode == 1)
      {
        pred = hasLeft ? reco[rowidx * stride - 1] : defaultPred;
      }
      else if (predMode == 2)
      {
        pred = hasTop ? reco[colidx - stride] : defaultPred;
      }
      else if (predMode == 3)
      {
        pred = (hasLeft && hasTop) ? ((reco[colidx - stride] + reco[rowidx * stride - 1]) >> 1) : defaultPred;
      }
      //store prediction signal
      if (dst != NULL)
      {
        dst[rowidx * stride + colidx] = pred;
      }
      //subtract or add prediction value
      if (resi != NULL)
      {
        resi[rowidx * stride + colidx] = clipLR(resi[rowidx * stride + colidx] - pred, -(1 << bitdepth) + 1, (1 << bitdepth) - 1);
      }
      else
      {
        reco[rowidx * stride + colidx] = clipLR(reco[rowidx * stride + colidx] + pred, 0, (1 << bitdepth) - 1);
      }
    }
  }
}

void transform(int* src, int width, int height, int stride, int bitdepthIn)
{
#if !USE_TAUBMANN
  //filter set sqrt(2):sqrt(2) normalization
  double h_syn[7] = { -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 };
  double g_ana[7] = { -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 };
  double h_ana[9] = {  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827 };
  double g_syn[9] = { -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827 };
#else
  //Taubmann:
  double h_syn[7] = { -0.091270, -0.057542,  0.591270, 1.115086,  0.591270, -0.057542, -0.091270 };
  double g_ana[7] = { 0.045635, -0.028771, -0.295635, 0.557543, -0.295635, -0.028771,  0.045635 };
  double h_ana[9] = { 0.026748, -0.016864, -0.078223,  0.266864, 0.602949,  0.266864, -0.078223, -0.016864, 0.026748 };
  double g_syn[9] = { 0.053496,  0.033728, -0.156446, -0.533728, 1.205898, -0.533728, -0.156446,  0.033728, 0.053496 };
#endif
  double* dsrc = (double*)malloc(width * height * sizeof(double));
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dsrc[rowidx * width + colidx] = (double)src[rowidx * stride + colidx];
    }
  }
#if USE_TAUBMANN
  convWT_2d(h_ana, 9, g_ana, 7, dsrc, width, height, width);
  int clipMin = -(1 << (bitdepthIn - 1)) + 1;
  int clipMax = (1 << (bitdepthIn - 1)) - 1;
#else
  lwt97_2d(dsrc, width, height, width);
  int clipMin = -(1 << bitdepthIn) + 1;
  int clipMax = (1 << bitdepthIn) - 1;
#endif
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      src[rowidx * stride + colidx] = clipLR((int)round(dsrc[rowidx * width + colidx]), clipMin, clipMax);
    }
  }
  free(dsrc);
}

void inv_transform(int* src, int width, int height, int stride, int bitdepthOut)
{
#if !USE_TAUBMANN
  //filter set sqrt(2):sqrt(2) normalization
  double h_syn[7] = { -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 };
  double g_ana[7] = { -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 };
  double h_ana[9] = {  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827 };
  double g_syn[9] = { -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827 };
#else
  //Taubmann:
  double h_syn[7] = { -0.091270, -0.057542,  0.591270, 1.115086,  0.591270, -0.057542, -0.091270 };
  double g_ana[7] = { 0.045635, -0.028771, -0.295635, 0.557543, -0.295635, -0.028771,  0.045635 };
  double h_ana[9] = { 0.026748, -0.016864, -0.078223,  0.266864, 0.602949,  0.266864, -0.078223, -0.016864, 0.026748 };
  double g_syn[9] = { 0.053496,  0.033728, -0.156446, -0.533728, 1.205898, -0.533728, -0.156446,  0.033728, 0.053496 };
#endif
  double* dsrc = (double*)malloc(width * height * sizeof(double));
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dsrc[rowidx * width + colidx] = (double)src[rowidx * stride + colidx];
    }
  }
#if USE_TAUBMANN
  invconvWT_2d(h_syn, 7, g_syn, 9, dsrc, width, height, width);
#else
  ilwt97_2d(dsrc, width, height, width);
#endif
  int clipMin = -(1 << (bitdepthOut - 1)) + 1;
  int clipMax = (1 << (bitdepthOut - 1)) - 1;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      src[rowidx * stride + colidx] = clipLR((int)round(dsrc[rowidx * width + colidx]), clipMin, clipMax);
    }
  }
  free(dsrc);
}

void quantize(int* src, int width, int height, int stride, int quantsize)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int sgn = src[rowidx * stride + colidx] < 0 ? -1 : 1;
      src[rowidx * stride + colidx] = sgn * (abs(src[rowidx * stride + colidx]) / quantsize);
    }
  }
}

void dequantize(int* src, int width, int height, int stride, int quantsize)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      src[rowidx * stride + colidx] *= quantsize;
    }
  }
}

unsigned long coded_bits(int* x, int width, int height, int stride, int bitdepth, int Qp)
{
  unsigned long bits = 0UL;
#if LL_FIXED_LENGTH
  //fixed length coding for the LL band
  int trafoBD = bitdepth + 2 - Qp;
  bits = (height / 2) * (width / 2) * (unsigned long)trafoBD;
#endif
  //Huffman coding depending on LL_FIXED_LENGTH: all four subbands or the three details bands only
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      if ((rowidx >= height/2 || colidx >= width/2) || !LL_FIXED_LENGTH)
      {
        int signFlag = 0;
        if (x[rowidx * stride + colidx] != 0)
        {
          signFlag = x[rowidx * stride + colidx] > 0 ? 1 : -1;
        }
        unsigned bitsCurr = (unsigned)(abs(x[rowidx * stride + colidx]) + 1 + abs(signFlag));
        bits += (unsigned long) bitsCurr;
      }
    }
  }
  return bits;
}

unsigned long mse_dist(int* src, int* reco, int width, int height, int stride)
{
  unsigned long dist = 0UL;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dist += ((src[rowidx * stride + colidx] - reco[rowidx * stride + colidx]) * (src[rowidx * stride + colidx] - reco[rowidx * stride + colidx]));
    }
  }
  return dist;
}

void encode_fixlen8(int* x, int width, int height, int stride, const char* encfile)
{
  FILE* fptr = fopen(encfile, "w");
  char* bitsOut = malloc(8);
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      memset(bitsOut, '0', 8);
      int n = x[rowidx * stride + colidx];
      int i = 0;
      assert(n < 256);
      while (n > 0)
      {
        bitsOut[7 - i] = (n % 2) + '0';
        n = n / 2;
        i++;
      }
      fprintf(fptr, "%s\n", (const char*)bitsOut);
    }
  }
  free(bitsOut);
  fclose(fptr);
}

void encode_huffman(int* x, int width, int height, int stride, const char* encfile)
{
  FILE* fptr = fopen(encfile, "w");
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int n = x[rowidx * stride + colidx];
      unsigned bitlen = (unsigned)(abs(n) + 1 + (n != 0));

      char* bitsOut = malloc(bitlen);
      if (n == 0)
      {
        bitsOut[0] = '0';
      }
      else
      {
        memset(bitsOut, '1', bitlen - 2);
        bitsOut[bitlen - 2] = '0';
        bitsOut[bitlen - 1] = n < 0 ? '0' : '1';
      }
      fprintf(fptr, "%s\n", (const char*)bitsOut);
      free(bitsOut);
    }
  }
  fclose(fptr);
}

void array_to_file(const char* fname, int* data, int len)
{
  FILE *fp = fopen(fname, "wb");
  if(fp == NULL)
  {
    printf("error creating file");
    exit(EXIT_FAILURE);
  }
  fwrite((const void*)data, sizeof(int), len, fp);
  fclose(fp);
}

void array_from_file(const char* fname, int* data, int len)
{
  FILE *fp = fopen(fname, "rb");
  if(fp == NULL)
  {
    printf("error opening file");
    exit(EXIT_FAILURE);
  }
  fread((void*)data, sizeof(int), len, fp);
  fclose(fp);
}

void blkcpy(int* src, int* dst, int width, int height, int stride)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    memcpy(dst + rowidx * stride, src + rowidx * stride, width * sizeof(int));
  }
}

int main(int argc, char **argv)
{
  if (argc != 7)
  {
    printf("Not enough or missing parameter!\n");
    printf("Usage: comp_demo <image_file> <width> <height> <prediction mode> <quant step-size> <partitioning depth>\n");
    return -1;
  }

  //declarations
  int width     = atoi(argv[2]);
  int height    = atoi(argv[3]);
  int predMode  = atoi(argv[4]);
  int stepSize  = atoi(argv[5]);
  int partDepth = atoi(argv[6]);
  int* x        = (int*)malloc(width * height * sizeof(int));
  int* pred     = (int*)malloc(width * height * sizeof(int));
  int* resi     = (int*)malloc(width * height * sizeof(int));
  int* trafo    = (int*)malloc(width * height * sizeof(int));
  int* quant    = (int*)malloc(width * height * sizeof(int));
  int* reco     = (int*)malloc(width * height * sizeof(int));

  //load data and copy to orig.bin
  array_from_file(argv[1], x, width * height);
  array_to_file("orig.bin", x, width * height);

  //estimate bit-depth and QP value
  int bitdepth = calcBitdepth(x, width * height);
  int QP       = calcBitdepth(&stepSize, 1) - 1;
  printf("processing input: bit-depth input image=%d, QP=%d, partitioning depth=%d\n", bitdepth, QP, partDepth);

  //adjust quantization step-size due to different transform scaling
#if USE_TAUBMANN
  int quantSize = stepSize >> 1;
#else
  int quantSize = stepSize;
#endif

  //rate-distortion parameter
  unsigned long bits = 0UL;
  unsigned long dist = 0UL;

  //partitioning
  int blkWidth  = MAX_BLOCK_SIZE >> partDepth;
  int blkHeight = MAX_BLOCK_SIZE >> partDepth;
  int blkStride = width / blkWidth;
  int numBlocks = (width / blkWidth) * (height / blkHeight);

  //quad split and for each block do first compression and then decompression
  for (int subblk = 0; subblk < numBlocks; subblk++)
  {
    bool leftMargin = true;
    bool topMargin = true;
    if (subblk == 0)
    {
      leftMargin = false;
      topMargin = false;
    }
    if (subblk != 0 && subblk / blkStride == 0) topMargin = false;
    if (subblk != 0 && subblk % blkStride == 0) leftMargin = false;
    int blkMode = -1;
    if (predMode < 4)
    {
      blkMode = predMode;
    }
    else if (predMode == 4)
    {
      blkMode = 2 * (int)topMargin + (int)leftMargin;
    }

    int offset = (subblk / blkStride) * blkHeight * width + (subblk % blkStride) * blkWidth;
    int* currOrig = x + offset;
    int* currPred = pred + offset;
    int* currResi = resi + offset;
    int* currTrafo = trafo + offset;
    int* currQuant = quant + offset;
    int* currReco = reco + offset;
    //compression: prediction, 9/7 transformtion and quatization
    blkcpy(currOrig, currResi, blkWidth, blkHeight, width);
    predict(blkMode, currReco, currPred, currResi, blkWidth, blkHeight, width, leftMargin, topMargin, bitdepth);
    blkcpy(currResi, currTrafo, blkWidth, blkHeight, width);
    transform(currTrafo, blkWidth, blkHeight, width, bitdepth + 1);
    blkcpy(currTrafo, currQuant, blkWidth, blkHeight, width);
    quantize(currQuant, blkWidth, blkHeight, width, quantSize);
    bits += coded_bits(currQuant, blkWidth, blkHeight, width, bitdepth, QP);
    blkcpy(currQuant, currReco, blkWidth, blkHeight, width);
    //decompression
    dequantize(currReco, blkWidth, blkHeight, width, quantSize);
    inv_transform(currReco, blkWidth, blkHeight, width, bitdepth + 1);
    predict(blkMode, currReco, NULL, NULL, blkWidth, blkHeight, width, leftMargin, topMargin, bitdepth);
    dist += mse_dist(currOrig, currReco, blkWidth, blkHeight, width);
  }

  //check quantization output in terms of bitdepth
  assert(calcBitdepth(quant, width * height) <= bitdepth + 2 - QP);

  //safe everything to files
  array_to_file("pred.bin", pred, width * height);
  array_to_file("resi.bin", resi, width * height);
  array_to_file("coeffs.bin", quant, width * height);
  array_to_file("reco.bin", reco, width * height);

  //encoding
  encode_huffman(quant, width, height, width, "enc_comp.txt");
  encode_fixlen8(x, width, height, width, "enc_orig.txt");

  printf("Relative distortion (MSE): %f\n", (double)dist / (double)(width * height));
  printf("Average symbol length (Bits): %f\n", (double)bits / (double)(width * height));
  printf("Compression rate: %f\n", 1.0 - (double)bits / (double)(bitdepth * width * height));

  free(x);
  free(pred);
  free(resi);
  free(trafo);
  free(quant);
  free(reco);
  return 0;
}