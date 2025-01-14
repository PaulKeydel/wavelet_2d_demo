#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "w97.h"

#define USE_TAUBMANN   0
#define MAX_BLOCK_SIZE 256

double clipLR(double val, int bitdepth, int shift)
{
  assert(shift < bitdepth);
  double min = 0;
  double max = (1 << bitdepth) - 1;
  double shift_val = (1 << shift) - 1;
  if (val + shift_val < min) return min - shift_val;
  if (val + shift_val > max) return max - shift_val;
  return val;
}

void predict(int predMode, double* reco, double* dst, int width, int height, int stride, int bitdepth)
{
  int defaultPred = 1 << (bitdepth - 1);
  int maxPixel = (1 << bitdepth) - 1;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      double pred = 0;
      if (predMode == 0)
      {
        pred = defaultPred;
      }
      else if (predMode == 1)
      {
        pred = reco[rowidx * stride - 1];
      }
      else if (predMode == 2)
      {
        pred = reco[colidx - stride];
      }
      else if (predMode == 3)
      {
        pred = 0.5 * (reco[colidx - stride] + reco[rowidx * stride - 1]);
      }
      //subtract or add prediction value
      if (dst != NULL)
      {
        dst[rowidx * stride + colidx] = clipLR(dst[rowidx * stride + colidx] - pred, bitdepth + 1 , bitdepth);
      }
      else
      {
        reco[rowidx * stride + colidx] = clipLR(reco[rowidx * stride + colidx] + pred, bitdepth, 0);
      }
    }
  }
}

void transform(double* x, int width, int height, int stride)
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
#if USE_TAUBMANN
  convWT_2d(h_ana, 9, g_ana, 7, x, width, height, stride);
#else
  lwt97_2d(x, width, height, stride);
#endif
}

void inv_transform(double* x, int width, int height, int stride)
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
#if USE_TAUBMANN
  invconvWT_2d(h_syn, 7, g_syn, 9, x, width, height, stride);
#else
  ilwt97_2d(x, width, height, stride);
#endif
}

int estBitdepth(double* x, int width, int height, int stride)
{
  double max = fabs(x[0]);
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      if (fabs(x[rowidx * stride + colidx]) > max) max = fabs(x[rowidx * stride + colidx]);
    }
  }
  return (int)ceil(log2(max));
}

void quantize(double* x, int width, int height, int stride, int bitdepth, int Qp)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int sgn = x[rowidx * stride + colidx] < 0 ? -1 : 1;
      x[rowidx * stride + colidx] = sgn * floor(fabs(x[rowidx * stride + colidx]) / (1 << Qp));
#if USE_TAUBMANN
      x[rowidx * stride + colidx] = clipLR(x[rowidx * stride + colidx], bitdepth + 1 - Qp, bitdepth - Qp);
#else
      x[rowidx * stride + colidx] = clipLR(x[rowidx * stride + colidx], bitdepth + 2 - Qp, bitdepth + 1 - Qp);
#endif
    }
  }
}

void dequantize(double* x, int width, int height, int stride, int Qp)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      x[rowidx * stride + colidx] *= (1 << Qp);
    }
  }
}

unsigned long coded_bits(double* x, int width, int height, int stride, int bitdepth, int Qp)
{
  unsigned long bits = 0UL;
#if USE_TAUBMANN
  int trafoBD = bitdepth + 1 - Qp;
#else
  int trafoBD = bitdepth + 2 - Qp;
#endif
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      if (rowidx < height/2 && colidx < width/2)
      {
        //fixed length coding for the LL band
        bits += (unsigned long)trafoBD;
      }
      else
      {
        //Huffman coding for all three details bands
        bits += (unsigned long)(x[rowidx * stride + colidx] + 1);
      }
    }
  }
  return bits;
}

unsigned long mse_dist(double* src, double* reco, int width, int height, int stride)
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

void dbls_to_file(const char* fname, double* data, int len)
{
  int* int_res = (int*)malloc(len * sizeof(int));
  for (int i = 0; i < len; i++) 
  {
    int_res[i] = (int)round(data[i]);
  }
  FILE *fp = fopen(fname, "wb");
  if(fp == NULL)
  {
    printf("error creating file");
  }
  fwrite((const void*)int_res, sizeof(int), len, fp);
  fclose(fp);
  free(int_res);
}

void dbls_from_file(const char* fname, double* data, int len)
{
  int* int_res = (int*)malloc(len * sizeof(int));
  FILE *fp = fopen(fname, "rb");
  if(fp == NULL)
  {
    printf("error opening file");
  }
  fread((void*)int_res, sizeof(int), len, fp);
  fclose(fp);
  for (int i = 0; i < len; i++) 
  {
    data[i] = (double)int_res[i];
  }
  free(int_res);
}

void blkcpy(double* src, double* dst, int width, int height, int stride)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    memcpy(dst + rowidx * stride, src + rowidx * stride, width * sizeof(double));
  }
}

int main(int argc, char **argv)
{
  if (argc != 6)
  {
    printf("Not enough or missing parameter!\n");
    printf("Usage: comp_demo <image_file> <width> <height> <quant step-size> <partitioning depth>\n");
    return -1;
  }

  //declarations
  int width     = atoi(argv[2]);
  int height    = atoi(argv[3]);
  int QP        = (int)round(log2(atoi(argv[4])));
  int partDepth = atoi(argv[5]);
  double* x     = (double*)malloc(width * height * sizeof(double));
  double* resi  = (double*)malloc(width * height * sizeof(double));
  double* trafo = (double*)malloc(width * height * sizeof(double));
  double* reco  = (double*)malloc(width * height * sizeof(double));

  //load data and copy to orig.bin
  dbls_from_file(argv[1], x, width * height);
  dbls_to_file("orig.bin", x, width * height);

  //estimate bit-depth
  int bitdepth = estBitdepth(x, width, height, width);
  printf("processing input: bit-depth input image=%d, QP=%d\n", bitdepth, QP);

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
    int predMode = 3;
    if (subblk == 0) predMode = 0;
    if (subblk / blkStride == 0) predMode = 1;
    if (subblk % blkStride == 0) predMode = 2;

    int offset = (subblk / blkStride) * blkHeight * width + (subblk % blkStride) * blkWidth;
    double* currOrig = x + offset;
    double* currResi = resi + offset;
    double* currTrafo = trafo + offset;
    double* currReco = reco + offset;
    //compression: prediction, 9/7 transformtion and quatization
    blkcpy(currOrig, currResi, blkWidth, blkHeight, width);
    predict(predMode, currReco, currResi, blkWidth, blkHeight, width, bitdepth);
    blkcpy(currResi, currTrafo, blkWidth, blkHeight, width);
    transform(currTrafo, blkWidth, blkHeight, width);
    quantize(currTrafo, blkWidth, blkHeight, width, bitdepth, QP);
    bits += coded_bits(currTrafo, blkWidth, blkHeight, width, bitdepth, QP);
    blkcpy(currTrafo, currReco, blkWidth, blkHeight, width);
    //decompression
    dequantize(currReco, blkWidth, blkHeight, width, QP);
    inv_transform(currReco, blkWidth, blkHeight, width);
    predict(predMode, currReco, NULL, blkWidth, blkHeight, width, bitdepth);
    dist += mse_dist(currOrig, currReco, blkWidth, blkHeight, width);
  }
  //safe everything to files
  dbls_to_file("resi.bin", resi, width * height);
  dbls_to_file("coeffs.bin", trafo, width * height);
  dbls_to_file("reco.bin", reco, width * height);

  printf("Relative distortion (MSE): %f\n", (double)dist / (double)(width * height));
  printf("Compression rate: %f\n", 1.0 - (double)bits / (double)(bitdepth * width * height));
  
  free(x);
  free(resi);
  free(trafo);
  free(reco);
  return 0;
}