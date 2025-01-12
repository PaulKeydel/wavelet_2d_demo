#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "w97.h"

#define USE_TAUBMANN 0

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

void predict(int subblk, double* reco, double* dst, int width, int height, int stride, int bitdepth)
{
  int defaultPred = 1 << (bitdepth - 1);
  int maxPixel = (1 << bitdepth) - 1;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      double pred = 0;
      if (subblk == 0)
      {
        pred = defaultPred;
      }
      else if (subblk == 1)
      {
        pred = reco[rowidx * stride - 1];
      }
      else if (subblk == 2)
      {
        pred = reco[colidx - stride];
      }
      else if (subblk == 3)
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

void transform(double* x, int width, int height, int stride, int shift)
{
#if !USE_TAUBMANN
  //filter set sqrt(2):sqrt(2) normalization
  double h_syn[7] = { -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 };
  double g_ana[7] = { -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 };
  double h_ana[9] = {  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827 };
  double g_syn[9] = { -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827 };
#else
  Taubmann:
  double h_syn[7] = { -0.091270, -0.057542,  0.591270, 1.115086,  0.591270, -0.057542, -0.091270 };
  double g_ana[7] = { 0.045635, -0.028771, -0.295635, 0.557543, -0.295635, -0.028771,  0.045635 };
  double h_ana[9] = { 0.026748, -0.016864, -0.078223,  0.266864, 0.602949,  0.266864, -0.078223, -0.016864, 0.026748 };
  double g_syn[9] = { 0.053496,  0.033728, -0.156446, -0.533728, 1.205898, -0.533728, -0.156446,  0.033728, 0.053496 };
#endif
  for (int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {
      x[i * stride + j] += (double)((1 << shift) - 1);
    }
  }
  lwt97_2d(x, width, height, stride);
}

void inv_transform(double* x, int width, int height, int stride, int shift)
{
#if !USE_TAUBMANN
  //filter set sqrt(2):sqrt(2) normalization
  double h_syn[7] = { -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 };
  double g_ana[7] = { -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 };
  double h_ana[9] = {  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827 };
  double g_syn[9] = { -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827 };
#else
  Taubmann:
  double h_syn[7] = { -0.091270, -0.057542,  0.591270, 1.115086,  0.591270, -0.057542, -0.091270 };
  double g_ana[7] = { 0.045635, -0.028771, -0.295635, 0.557543, -0.295635, -0.028771,  0.045635 };
  double h_ana[9] = { 0.026748, -0.016864, -0.078223,  0.266864, 0.602949,  0.266864, -0.078223, -0.016864, 0.026748 };
  double g_syn[9] = { 0.053496,  0.033728, -0.156446, -0.533728, 1.205898, -0.533728, -0.156446,  0.033728, 0.053496 };
#endif
  invconvWT_2d(h_syn, 7, g_syn, 9, x, width, height, stride);
  for (int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {
      x[i * stride + j] -= (double)((1 << shift) - 1);
    }
  }
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
      x[rowidx * stride + colidx] = clipLR(x[rowidx * stride + colidx], bitdepth + 2 - Qp, 0);
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
  if (argc != 5)
  {
    printf("Not enough parameter!\n");
    printf("Usage: comp_demo <image_file> <width> <height> <quant step-size>\n");
    return -1;
  }

  //declarations
  int width = atoi(argv[2]);
  int height = atoi(argv[3]);
  int QP = (int)round(log2(atoi(argv[4])));
  double* x = (double*)malloc(width * height * sizeof(double));
  double* resi = (double*)malloc(width * height * sizeof(double));
  double* trafo = (double*)malloc(width * height * sizeof(double));
  double* reco = (double*)malloc(width * height * sizeof(double));

  //load data and copy to orig.bin
  dbls_from_file(argv[1], x, width * height);
  dbls_to_file("orig.bin", x, width * height);

  //estimate bit-depth
  int bitdepth = estBitdepth(x, width, height, width);
  printf("processing input: bit-depth input image=%d, QP=%d\n", bitdepth, QP);

  //quad split and for each block do first compression and then decompression
  for (int blk = 0; blk < 4; blk++)
  {
    //depth-1 partitioning
    int offset = ((blk % 2) * (width / 2)) + ((blk / 2) * (height * width / 2));
    for (int subblk = blk; subblk < (blk == 0 ? 4 : (blk + 1)); subblk++)
    {
      int suboffset = offset;
      int blkratio = 2;
      if (blk == 0)
      {
        //depth-2 partitioning
        suboffset = ((subblk % 2) * (width / 4)) + ((subblk / 2) * (height * width / 4));
        blkratio = 4;
      }
      double* currOrig = x + suboffset;
      double* currResi = resi + suboffset;
      double* currTrafo = trafo + suboffset;
      double* currReco = reco + suboffset;
      //compression: prediction, 9/7 transformtion and quatization
      blkcpy(currOrig, currResi, width/blkratio, height/blkratio, width);
      predict(subblk, currReco, currResi, width/blkratio, height/blkratio, width, bitdepth);
      blkcpy(currResi, currTrafo, width/blkratio, height/blkratio, width);
      transform(currTrafo, width/blkratio, height/blkratio, width, bitdepth);
      quantize(currTrafo, width/blkratio, height/blkratio, width, bitdepth, QP);
      blkcpy(currTrafo, currReco, width/blkratio, height/blkratio, width);
      //decompression
      dequantize(currReco, width/blkratio, height/blkratio, width, QP);
      inv_transform(currReco, width/blkratio, height/blkratio, width, bitdepth);
      predict(subblk, currReco, NULL, width/blkratio, height/blkratio, width, bitdepth);
    }
  }
  //safe everything to files
  dbls_to_file("resi.bin", resi, width * height);
  dbls_to_file("coeffs.bin", trafo, width * height);
  dbls_to_file("reco.bin", reco, width * height);
  
  free(x);
  free(resi);
  free(trafo);
  free(reco);
  return 0;
}