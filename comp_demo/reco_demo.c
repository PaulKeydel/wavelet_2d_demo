#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lib_enc_dec.h"


unsigned get_fsize(const char* fname)
{
  FILE* fp = fopen(fname, "rb");
  if(fp == NULL)
  {
    printf("Error opening file");
    exit(EXIT_FAILURE);
  }
  fseek(fp, 0L, SEEK_END);
  long sz = ftell(fp);
  fseek(fp, 0L, SEEK_SET);
  return (unsigned)sz;
}

void reconstruct(int* quant, int* reco, int width, int height, int stepSize, uchar* byteStream, unsigned* binLen)
{
  int bestPred    = 0;
  int bestDepth   = 0;
  int bestCutting = 0;

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

    *binLen += decode_unit(byteStream, *binLen, quant + offset, width, &bestDepth, &bestPred, &bestCutting);
    //printf("ui: %d\n", ui);
    //printf("  cutting: %d\n", bestCutting);
    //printf("  pred: %d\n", bestPred);
    reconstruct_unit(quant + offset, reco + offset, width, stepSize, bestDepth, bestPred, bestCutting, topMargin, leftMargin);
  }
}


int main(int argc, char **argv)
{
  if (argc != 2)
  {
    printf("Not enough or missing parameter!\n");
    printf("Usage: reco_demo <bitstream-file>\n");
    return -1;
  }
  int width;
  int height;
  int stepSize;

  unsigned numBytes = get_fsize(argv[1]);

  uchar* byteStream = (uchar*)malloc(numBytes);
  array_from_file(argv[1], (void*)byteStream, sizeof(uchar), numBytes);

  unsigned binLen = decode_coding_params(byteStream, 0, &width, &height, &stepSize);

  int* quant = (int*)malloc(width * height * sizeof(int));
  int* reco  = (int*)malloc(width * height * sizeof(int));

  reconstruct(quant, reco, width, height, stepSize, byteStream, &binLen);

  array_to_file("outputs_dec/coeffs.bin", (const void*)quant, sizeof(int), width * height);
  array_to_file("outputs_dec/reco.bin", (const void*)reco, sizeof(int), width * height);

  free(byteStream);
  free(quant);
  free(reco);
  return 0;
}