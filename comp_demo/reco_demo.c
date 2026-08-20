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

  unsigned binLen = 0U;

  const int numPixels4K = 3840 * 2160;
  int* quant = (int*)malloc(numPixels4K * sizeof(int));
  int* reco  = (int*)malloc(numPixels4K * sizeof(int));

  reconstruct_image(quant, reco, &width, &height, &stepSize, byteStream, &binLen);

  array_to_file("outputs_dec/coeffs.bin", (const void*)quant, sizeof(int), width * height);
  array_to_file("outputs_dec/reco.bin", (const void*)reco, sizeof(int), width * height);

  free(byteStream);
  free(quant);
  free(reco);
  return 0;
}