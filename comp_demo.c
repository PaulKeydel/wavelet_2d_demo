#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "w97.h"

void predict(int blk, double* reco, double* dst, int width, int height, int stride)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int pred = 0;
      if (blk == 0)
      {
        pred = 32;
      }
      else if (blk == 1)
      {
        pred = reco[rowidx * stride - 1];
      }
      else if (blk == 2)
      {
        pred = reco[colidx - stride];
      }
      else if (blk == 3)
      {
        pred = 0.5 * (reco[colidx - stride] + reco[rowidx * stride - 1]);
      }
      //subtract or add prediction value
      if (dst != NULL)
      {
        dst[rowidx * stride + colidx] -= pred;
      }
      else
      {
        reco[rowidx * stride + colidx] += pred;
      }
    }
  }
}

int estQuantStepSize(double* x, int width, int height, int stride)
{
  double min = x[0];
  double max = x[0];
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      if (x[rowidx * stride + colidx] > max) max = x[rowidx * stride + colidx];
      if (x[rowidx * stride + colidx] < min) min = x[rowidx * stride + colidx];
    }
  }
  int bitdepth = (int)ceil(log2(max - min));
  return (int)pow(2, bitdepth / 2 + 1);
}

void quantize(double* x, int width, int height, int stride, int stepsize)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int sgn = x[rowidx * stride + colidx] < 0 ? -1 : 1;
      x[rowidx * stride + colidx] = sgn * floor(fabs(x[rowidx * stride + colidx]) / stepsize);
    }
  }
}

void dequantize(double* x, int width, int height, int stride, int stepsize)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      x[rowidx * stride + colidx] *= stepsize;
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
  if (argc != 4)
  {
    printf("Not enough parameter!\n");
    printf("Usage: w97 <image_file> <width> <height>\n");
    return -1;
  }

  //declarations
  int width = atoi(argv[2]);
  int height = atoi(argv[3]);
  double* x = (double*)malloc(width * height * sizeof(double));
  double* resi = (double*)malloc(width * height * sizeof(double));
  double* trafo = (double*)malloc(width * height * sizeof(double));
  double* reco = (double*)malloc(width * height * sizeof(double));

  //load data and copy to orig.bin
  dbls_from_file(argv[1], x, width * height);
  dbls_to_file("orig.bin", x, width * height);
    
  //filter set sqrt(2):sqrt(2) normalization
  double h_syn[7] = { -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 };
  double g_ana[7] = { -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 };
  double h_ana[9] = {  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827};
  double g_syn[9] = { -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827};
    
  //Taubmann:
  //h_syn: -0.091270, -0.057542,  0.591270, 1.115086,  0.591270, -0.057542, -0.091270
  //g_ana:  0.045635, -0.028771, -0.295635, 0.557543, -0.295635, -0.028771,  0.045635
  //h_ana: 0.026748, -0.016864, -0.078223,  0.266864, 0.602949,  0.266864, -0.078223, -0.016864, 0.026748
  //g_syn: 0.053496,  0.033728, -0.156446, -0.533728, 1.205898, -0.533728, -0.156446,  0.033728, 0.053496

  //estimate step size for quantization depending on input bitdepth
  int stepsize0 = estQuantStepSize(x, width, height, width);
  printf("Quantization step size: %d\n", stepsize0);

  //quad split and for each block do first compression and then decompression
  for (int blk = 0; blk < 4; blk++)
  {
    //block partitioning
    int offset = ((blk % 2) * (width / 2)) + ((blk / 2) * (height * width / 2));
    double* currOrig = x + offset;
    double* currResi = resi + offset;
    double* currTrafo = trafo + offset;
    double* currReco = reco + offset;
    //compression: prediction, 9/7 transformtion and quatization
    blkcpy(currOrig, currResi, width/2, height/2, width);
    predict(blk, currReco, currResi, width/2, height/2, width);
    blkcpy(currResi, currTrafo, width/2, height/2, width);
    lwt97_2d(currTrafo, width/2, height/2, width);
    quantize(currTrafo, width/2, height/2, width, stepsize0);
    blkcpy(currTrafo, currReco, width/2, height/2, width);
    //decompression
    dequantize(currReco, width/2, height/2, width, stepsize0);
    invconvWT_2d(h_syn, 7, g_syn, 9, currReco, width/2, height/2, width);
    predict(blk, currReco, NULL, width/2, height/2, width);
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