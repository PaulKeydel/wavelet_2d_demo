#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "lib_enc_dec.h"


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

void compress(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int bitdepth, int stepSize, double lambda, unsigned long* totalBits, unsigned long* totalDist, FILE* fptrEncTxt, FILE* fptrEncBin)
{
  int numUnitsX = width / MAX_BLOCK_SIZE;
  int numUnitsY = height / MAX_BLOCK_SIZE;
  int numUnits  = numUnitsX * numUnitsY;
  for (int ui = 0; ui < numUnits; ui++)
  {
    bool leftMargin = true;
    bool topMargin = true;
    if (ui == 0)
    {
      leftMargin = false;
      topMargin = false;
    }
    if (ui != 0 && ui / numUnitsX == 0) topMargin = false;
    if (ui != 0 && ui % numUnitsX == 0) leftMargin = false;

    int offset = (ui / numUnitsX) * MAX_BLOCK_SIZE * width + (ui % numUnitsX) * MAX_BLOCK_SIZE;
    double bestCost = MAXFLOAT;
    int bestPred    = 0;
    int bestDepth   = 0;
    for (int partDepth = 0; partDepth < 5; partDepth++)
    {
      for (int predMode = 0; predMode < 4; predMode++)
      {
        unsigned long blkBits = 0UL;
        unsigned long blkDist = 0UL;
        compress_unit(
          x + offset,
          pred + offset,
          resi + offset,
          trafo + offset,
          quant + offset,
          reco + offset,
          bitdepth,
          width,
          stepSize,
          partDepth,
          predMode,
          topMargin,
          leftMargin,
          &blkBits,
          &blkDist
        );
        double blkCost = (double)blkDist + lambda * (double)blkBits;
        if (blkCost < bestCost)
        {
          bestDepth = partDepth;
          bestPred = predMode;
          bestCost = blkCost;
        }
      }
    }
    compress_unit(
      x + offset,
      pred + offset,
      resi + offset,
      trafo + offset,
      quant + offset,
      reco + offset,
      bitdepth,
      width,
      stepSize,
      bestDepth,
      bestPred,
      topMargin,
      leftMargin,
      totalBits,
      totalDist
    );
    encode_unit(quant + offset, bitdepth, width, stepSize, bestDepth, bestPred, fptrEncTxt, fptrEncBin);
  }
}


int main(int argc, char **argv)
{
  if (argc != 5 && argc != 6)
  {
    printf("Not enough or missing parameter!\n");
    printf("Usage: comp_demo <image_file> <width> <height> <quant step-size> <optional: lambda>\n");
    return -1;
  }

  //declarations
  int width     = atoi(argv[2]);
  int height    = atoi(argv[3]);
  int stepSize  = atoi(argv[4]);
  double lambda = argc == 6 ? atof(argv[5]) : calcLambda(stepSize);
  int* x        = (int*)malloc(width * height * sizeof(int));
  int* pred     = (int*)malloc(width * height * sizeof(int));
  int* resi     = (int*)malloc(width * height * sizeof(int));
  int* trafo    = (int*)malloc(width * height * sizeof(int));
  int* quant    = (int*)malloc(width * height * sizeof(int));
  int* reco     = (int*)malloc(width * height * sizeof(int));

  //load data and copy to orig.bin
  array_from_file(argv[1], x, width * height);
  array_to_file("outputs/orig.bin", x, width * height);

  //estimate bit-depth and QP value
  int bitdepth = calcBitdepth(x, width * height);
  int QP       = calcBitdepth(&stepSize, 1) - 1;
  printf("processing input: bit-depth input image=%d, QP=%d, Lambda=%f\n", bitdepth, QP, lambda);

  //rate-distortion parameter
  unsigned long totalBits = 0UL;
  unsigned long totalDist = 0UL;

  //output text file for resulting bitstream
  FILE* fptrEncTxt = fopen("outputs/bitstream.txt", "w");
  FILE* fptrEncBin = fopen("outputs/bitstream.bin", "wb");

  //find RD-optimized compression mode and do compression + encoding
  compress(
    x,
    pred,
    resi,
    trafo,
    quant,
    reco,
    width,
    height,
    bitdepth,
    stepSize,
    lambda,
    &totalBits,
    &totalDist,
    fptrEncTxt,
    fptrEncBin
  );
  //check quantization output in terms of bitdepth
  assert(calcBitdepth(quant, width * height) <= bitdepth + 1 - QP);

  //safe everything to files
  array_to_file("outputs/pred.bin", pred, width * height);
  array_to_file("outputs/resi.bin", resi, width * height);
  array_to_file("outputs/coeffs.bin", quant, width * height);
  array_to_file("outputs/reco.bin", reco, width * height);

  //free file buffer
  fclose(fptrEncTxt);
  fclose(fptrEncBin);

  printf("Relative distortion (MSE): %f\n", (double)totalDist / (double)(width * height));
  printf("Average symbol length (Bits): %f\n", (double)totalBits / (double)(width * height));
  printf("Compression rate: %f\n", 1.0 - (double)totalBits / (double)(bitdepth * width * height));

  free(x);
  free(pred);
  free(resi);
  free(trafo);
  free(quant);
  free(reco);
  return 0;
}