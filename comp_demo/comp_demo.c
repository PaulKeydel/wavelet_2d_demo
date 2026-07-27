#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "lib_enc_dec.h"


void compress(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int stepSize, double lambda, uchar* byteStream, unsigned* binLen)
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

    compress_unit(x + offset, pred + offset, resi + offset, trafo + offset, quant + offset, reco + offset, width, stepSize, topMargin, leftMargin, lambda, byteStream, binLen);
  }
  assert(*binLen % 8 == 0);
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
  array_from_file(argv[1], (void*)x, sizeof(int), width * height);
  array_to_file("outputs_enc/orig.bin", (const void*)x, sizeof(int), width * height);

  //estimate bit-depth and QP value
  int bitdepth = calcBitdepth(x, width * height);
  int QP       = calcBitdepth(&stepSize, 1) - 1;

  if (bitdepth > IMG_BITDEPTH)
  {
    printf("Actual bitdepth of input image is greater than configured bitdepth!\n");
    return -1;
  }
  printf("processing input: bit-depth input image=%d, QP=%d, Lambda=%f\n", bitdepth, QP, lambda);

  //buffer for resulting bitstream
  int maxQuant      = (1 << (bitdepth + 1 - QP)) - 1;
  int numUnits      = (width / MAX_BLOCK_SIZE) * (height / MAX_BLOCK_SIZE);
  int numSubUnits   = (1 << ENC_MAX_DEPTH) * (1 << ENC_MAX_DEPTH);
  int numBytes      = (32 + 3 * numUnits + 4 * numSubUnits + (maxQuant + 2) * MAX_BLOCK_SIZE * MAX_BLOCK_SIZE) / 8;
  uchar* byteStream = (uchar*)malloc(numBytes);
  memset(byteStream, 0, numBytes);

  //find RD-optimized compression mode and do compression + encoding
  unsigned binLen = encode_coding_params(width, height, stepSize, byteStream, 0);
  compress(x, pred, resi, trafo, quant, reco, width, height, stepSize, lambda, byteStream, &binLen);
  //check quantization output in terms of bitdepth
  assert(calcBitdepth(quant, width * height) <= bitdepth + 1 - QP);

  //safe everything to files
  array_to_file("outputs_enc/pred.bin", (const void*)pred, sizeof(int), width * height);
  array_to_file("outputs_enc/resi.bin", (const void*)resi, sizeof(int), width * height);
  array_to_file("outputs_enc/coeffs.bin", (const void*)quant, sizeof(int), width * height);
  array_to_file("outputs_enc/reco.bin", (const void*)reco, sizeof(int), width * height);
  array_to_file("outputs_enc/bitstream.bin", (const void*)byteStream, sizeof(uchar), binLen / 8);

  //rate-distortion summary where we neglect fixlen-coded parameters
  double totalDist = (double)mse_dist(x, reco, width, height, width);
  double totalBits = (double)binLen;
  printf("Relative distortion (MSE): %f\n", totalDist / (double)(width * height));
  printf("Average symbol length (Bits): %f\n", totalBits / (double)(width * height));
  printf("Compression rate: %f\n", 1.0 - totalBits / (double)(bitdepth * width * height));

  free(byteStream);
  free(x);
  free(pred);
  free(resi);
  free(trafo);
  free(quant);
  free(reco);
  return 0;
}