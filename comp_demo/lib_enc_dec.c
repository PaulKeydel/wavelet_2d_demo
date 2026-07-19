#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lib_enc_dec.h"
#include "w97.h"


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

double calcLambda(int stepSize)
{
  //double res = 3.392 * stepSize * stepSize - 24.93 * stepSize + 53.88;
  double res = 2.536 * stepSize * stepSize - 21.27 * stepSize + 76.61;
  //double res = 2.898 * stepSize * stepSize - 29.36 * stepSize + 113.5;
  res = res < 5 ? 5 : res;
  return res;
}

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

void blkcpy(int* src, int* dst, int width, int height, int stride)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    memcpy(dst + rowidx * stride, src + rowidx * stride, width * sizeof(int));
  }
}

void array_to_file(const char* fname, const void* data, int typeSize, int len)
{
  FILE* fp = fopen(fname, "wb");
  if(fp == NULL)
  {
    printf("error creating file");
    exit(EXIT_FAILURE);
  }
  fwrite(data, typeSize, len, fp);
  fclose(fp);
}

void array_from_file(const char* fname, void* data, int typeSize, int len)
{
  FILE* fp = fopen(fname, "rb");
  if(fp == NULL)
  {
    printf("Error opening file");
    exit(EXIT_FAILURE);
  }
  fread(data, typeSize, len, fp);
  fclose(fp);
}

void predict(int predMode, int* reco, int* dst, int* resi, int width, int height, int stride, bool hasLeft, bool hasTop, int bitdepth)
{
  int defaultPred = 1 << (bitdepth - 1);
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int pred = 0;
      if (predMode == PRED_CONST)
      {
        pred = defaultPred;
      }
      else if (predMode == PRED_HOR)
      {
        pred = hasLeft ? reco[rowidx * stride - 1] : defaultPred;
      }
      else if (predMode == PRED_VER)
      {
        pred = hasTop ? reco[colidx - stride] : defaultPred;
      }
      else if (predMode == PRED_DIAG)
      {
        pred = ((hasLeft ? reco[rowidx * stride - 1] : defaultPred) + (hasTop ? reco[colidx - stride] : defaultPred)) >> 1;
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
  double* dsrc = (double*)malloc(width * height * sizeof(double));
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dsrc[rowidx * width + colidx] = (double)src[rowidx * stride + colidx];
    }
  }
#if USE_TAUBMANN
  convWT_2d(&FilterTaubmann, dsrc, width, height, width);
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
  double* dsrc = (double*)malloc(width * height * sizeof(double));
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dsrc[rowidx * width + colidx] = (double)src[rowidx * stride + colidx];
    }
  }
#if USE_TAUBMANN
  invconvWT_2d(&FilterTaubmann, dsrc, width, height, width);
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

unsigned long rd_unit_bits(int* x, int stride, int partDepth, int cutCoefs)
{
  unsigned long bits = 0UL;
  int blkSize = MAX_BLOCK_SIZE >> partDepth;
  for (int rowidx = 0; rowidx < MAX_BLOCK_SIZE; rowidx++)
  {
    for (int colidx = 0; colidx < MAX_BLOCK_SIZE; colidx++)
    {
      //Huffman coding all relevant subbands
      if (cutCoefs == 0 || subband(rowidx % blkSize, colidx % blkSize, blkSize, blkSize, stride) != cutCoefs)
      {
        int absVal = abs(x[rowidx * stride + colidx]);
        bits += (unsigned long)(absVal == 0 ? 1 : absVal + 2);
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

unsigned encode_fixlen(int n, int len, uchar* bitsOut, int bitPos)
{
  assert(len >= calcBitdepth(&n, 1));
  for (int i = 0; i < len; i++)
  {
    setBit(bitsOut, (bitPos + i), 0);
  }
  int i = 0;
  while (n > 0)
  {
    setBit(bitsOut, (bitPos + len - i - 1), (n % 2));
    n = n / 2;
    i++;
  }
  return len;
}

unsigned decode_fixlen(int len, uchar* bitsIn, int bitPos, int* n)
{
  *n = 0;
  for (int i = 0; i < len; i++)
  {
    *n += checkBit(bitsIn, (bitPos + len - 1 - i)) * (1 << i);
  }
  return len;
}

unsigned encode_huffman(int n, uchar* bitsOut, int bitPos)
{
  unsigned bitlen = (unsigned)(abs(n) + 1 + (n != 0));
  if (n == 0)
  {
    setBit(bitsOut, bitPos, 0);
  }
  else
  {
    for (int i = 0; i < bitlen - 2; i++)
    {
      setBit(bitsOut, (bitPos + i), 1);
    }
    setBit(bitsOut, (bitPos + bitlen - 2), 0);
    setBit(bitsOut, (bitPos + bitlen - 1), (n < 0 ? 0 : 1));
  }
  return bitlen;
}

unsigned decode_huffman(uchar* bitsIn, int bitPos, int* n)
{
  int i = 0;
  while (checkBit(bitsIn, (bitPos + i)) != 0)
  {
    i++;
  }
  *n = i;
  if (i > 0)
  {
    *n *= (checkBit(bitsIn, (bitPos + i + 1)) == 0 ? -1 : 1);
  }
  unsigned bitlen = (i == 0 ? 1U : (unsigned)i + 2U);
  return bitlen;
}

unsigned encode_coding_params(int width, int height, int stepSize, uchar* bitsOut, int bitPos)
{
  encode_fixlen(width, 12, bitsOut, bitPos);
  encode_fixlen(height, 12, bitsOut, bitPos + 12);
  encode_fixlen(stepSize, 8, bitsOut, bitPos + 24);
  return 32;
}

unsigned decode_coding_params(uchar* bitsIn, int bitPos, int* width, int* height, int* stepSize)
{
  decode_fixlen(12, bitsIn, bitPos, width);
  decode_fixlen(12, bitsIn, bitPos + 12, height);
  decode_fixlen(8, bitsIn, bitPos + 24, stepSize);
  return 32;
}

unsigned encode_unit(int* quant, int stride, int partDepth, int predMode, int cutCoefs, uchar* binStream, int bitPos)
{
  unsigned binLen = 0;
  binLen += encode_fixlen(partDepth, 3, binStream, bitPos + binLen);
  binLen += encode_fixlen(predMode, 2, binStream, bitPos + binLen);
  binLen += encode_fixlen(cutCoefs, 2, binStream, bitPos + binLen);

  int blkSize = MAX_BLOCK_SIZE >> partDepth;
  for (int rowidx = 0; rowidx < MAX_BLOCK_SIZE; rowidx++)
  {
    for (int colidx = 0; colidx < MAX_BLOCK_SIZE; colidx++)
    {
      if (cutCoefs == 0 || subband(rowidx % blkSize, colidx % blkSize, blkSize, blkSize, stride) != cutCoefs)
      {
        int symbol = *(quant + rowidx * stride + colidx);
        binLen += encode_huffman(symbol, binStream, bitPos + binLen);
      }
    }
  }
  if (binLen % 8 > 0)
  {
    binLen += (8 - (binLen % 8));
  }
  return binLen;
}

unsigned decode_unit(uchar* binStream, int bitPos, int* quant, int stride, int* partDepth, int* predMode, int* cutCoefs)
{
  unsigned binLen = 0;
  binLen += decode_fixlen(3, binStream, bitPos + binLen, partDepth);
  binLen += decode_fixlen(2, binStream, bitPos + binLen, predMode);
  binLen += decode_fixlen(2, binStream, bitPos + binLen, cutCoefs);

  int blkSize = MAX_BLOCK_SIZE >> *partDepth;
  for (int rowidx = 0; rowidx < MAX_BLOCK_SIZE; rowidx++)
  {
    for (int colidx = 0; colidx < MAX_BLOCK_SIZE; colidx++)
    {
      if (*cutCoefs && subband(rowidx % blkSize, colidx % blkSize, blkSize, blkSize, stride) == *cutCoefs)
      {
        *(quant + rowidx * stride + colidx) = 0;
      }
      else
      {
        binLen += decode_huffman(binStream, bitPos + binLen, quant + rowidx * stride + colidx);
      }
    }
  }
  if (binLen % 8 != 0)
  {
    binLen += 8 - (binLen % 8);
  }
  assert(binLen % 8 == 0);
  return binLen;
}

void cut_detail_coefs(int* src, int width, int height, int stride, int cutCoefs)
{
  if (cutCoefs == 0)
  {
    return;
  }
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      if (subband(rowidx, colidx, width, height, stride) == cutCoefs)
      {
        src[rowidx * stride + colidx] = 0;
      }
    }
  }
}

void comp_reco_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int bitdepth, int stride, int stepSize, int partDepth, int predMode, int cutCoefs, bool topMargin, bool leftMargin)
{
  //adjust quantization step-size (due to different transform scalings)
#if USE_TAUBMANN
  int quantSize = stepSize;
#else
  int quantSize = stepSize << 1;
#endif
  //split and for each subblock do first compression and then decompression
  int blkWidth  = MAX_BLOCK_SIZE >> partDepth;
  int blkHeight = MAX_BLOCK_SIZE >> partDepth;
  int blkStride = 1 << partDepth;
  int blkNum    = (1 << partDepth) * (1 << partDepth);
  for (int subblk = 0; subblk < blkNum; subblk++)
  {
    bool hasLeft = true;
    bool hasTop = true;
    if (subblk == 0)
    {
      hasLeft = leftMargin;
      hasTop = topMargin;
    }
    else if (subblk / blkStride == 0)
    {
      hasTop = topMargin;
    }
    else if (subblk % blkStride == 0)
    {
      hasLeft = leftMargin;
    }

    int offset = (subblk / blkStride) * blkHeight * stride + (subblk % blkStride) * blkWidth;
    int* currQuant = quant + offset;
    int* currReco = reco + offset;
    //compression: prediction, 9/7 transformtion and quatization
    if (x != NULL && pred != NULL && resi != NULL && trafo != NULL)
    {
      int* currOrig = x + offset;
      int* currPred = pred + offset;
      int* currResi = resi + offset;
      int* currTrafo = trafo + offset;
      blkcpy(currOrig, currResi, blkWidth, blkHeight, stride);
      predict(predMode, currReco, currPred, currResi, blkWidth, blkHeight, stride, hasLeft, hasTop, bitdepth);
      blkcpy(currResi, currTrafo, blkWidth, blkHeight, stride);
      transform(currTrafo, blkWidth, blkHeight, stride, bitdepth + 1);
      cut_detail_coefs(currTrafo, blkWidth, blkHeight, stride, cutCoefs);
      blkcpy(currTrafo, currQuant, blkWidth, blkHeight, stride);
      quantize(currQuant, blkWidth, blkHeight, stride, quantSize);
    }
    //decompression
    blkcpy(currQuant, currReco, blkWidth, blkHeight, stride);
    dequantize(currReco, blkWidth, blkHeight, stride, quantSize);
    inv_transform(currReco, blkWidth, blkHeight, stride, bitdepth + 1);
    predict(predMode, currReco, NULL, NULL, blkWidth, blkHeight, stride, hasLeft, hasTop, bitdepth);
  }
}

void compress_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int stepSize, int partDepth, int predMode, int cutCoefs, bool topMargin, bool leftMargin, unsigned long* rdBits, unsigned long* rdDist)
{
  comp_reco_unit(x, pred, resi, trafo, quant, reco, IMG_BITDEPTH, stride, stepSize, partDepth, predMode, cutCoefs, topMargin, leftMargin);
  //estimate bits per unit and calculate MSE
  *rdBits += rd_unit_bits(quant, stride, partDepth, cutCoefs);
  *rdDist += mse_dist(x, reco, MAX_BLOCK_SIZE, MAX_BLOCK_SIZE, stride);
}

void reconstruct_unit(int* quant, int* reco, int stride, int stepSize, int partDepth, int predMode, int cutCoefs, bool topMargin, bool leftMargin)
{
  comp_reco_unit(NULL, NULL, NULL, NULL, quant, reco, IMG_BITDEPTH, stride, stepSize, partDepth, predMode, cutCoefs, topMargin, leftMargin);
}