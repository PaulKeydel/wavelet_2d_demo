#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lib_enc_dec.h"
#include "w97.h"
#include "lib_bin_coding.h"


double calcLambda(int stepSize)
{
  double lambda = 2 * (stepSize - 4) * (stepSize - 4) + 28;
  lambda = stepSize < 4 ? 28 : lambda;
  return lambda;
}

int calcBitdepth(int* src, int n)
{
  int max = src[0];
  int min = src[0];
  for (int i = 1; i < n; i++)
  {
    if (src[i] > max) max = src[i];
    if (src[i] < min) min = src[i];
  }
  int ref = (n == 1 ? abs(*src) : (max - min));
  return (int)floor(log2((double)ref)) + 1;
}

ulong dist_ssd(int* src, int* ref, int width, int height, int stride)
{
  ulong dist = 0UL;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dist += ((src[rowidx * stride + colidx] - ref[rowidx * stride + colidx]) * (src[rowidx * stride + colidx] - ref[rowidx * stride + colidx]));
    }
  }
  return dist;
}

void clipLR(int* src, int width, int height, int stride, int min, int max)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int val = src[rowidx * stride + colidx];
      if (val > max)
      {
        val = max;
      }
      if (val < min)
      {
        val = min;
      }
      src[rowidx * stride + colidx] = val;
    }
  }
}

int subband_at_pos(int row, int col, int width, int height, int stride)
{
  if ((row < height / 2) && (col < width / 2))
  {
    return SUBBAND_LL;
  }
  else if ((row < height / 2) && (col >= width / 2))
  {
    return SUBBAND_LH;
  }
  else if ((row >= height / 2) && (col < width / 2))
  {
    return SUBBAND_HL;
  }
  else
  {
    return SUBBAND_HH;
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

void predict(int predMode, int* reco, int* pred, int* resi, int width, int height, int stride, bool hasLeft, bool hasTop)
{
  int defaultPred = 1 << (IMG_BITDEPTH - 1);
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int predVal = defaultPred;
      if (predMode == PRED_HOR)
      {
        predVal = hasLeft ? reco[rowidx * stride - 1] : defaultPred;
      }
      else if (predMode == PRED_VER)
      {
        predVal = hasTop ? reco[colidx - stride] : defaultPred;
      }
      else if (predMode == PRED_DIAG)
      {
        predVal = ((hasLeft ? reco[rowidx * stride - 1] : defaultPred) + (hasTop ? reco[colidx - stride] : defaultPred)) >> 1;
      }
      //store prediction signal and calculate residual
      if (pred != NULL && resi != NULL)
      {
        pred[rowidx * stride + colidx] = predVal;
        resi[rowidx * stride + colidx] = resi[rowidx * stride + colidx] - predVal;
      }
      //update current reconstruction value
      else
      {
        reco[rowidx * stride + colidx] += predVal;
      }
    }
  }
}

void transform(int* src, int width, int height, int stride)
{
  double* dsrc = (double*)malloc(width * height * sizeof(double));
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dsrc[rowidx * width + colidx] = (double)src[rowidx * stride + colidx];
    }
  }
#if USE_WAVE_LIFTING
  lwt97_2d(dsrc, width, height, width);
#else
  convWT_2d(dsrc, width, height, width);
#endif
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      src[rowidx * stride + colidx] = (int)round(dsrc[rowidx * width + colidx]);
    }
  }
  free(dsrc);
}

void inv_transform(int* src, int width, int height, int stride)
{
  double* dsrc = (double*)malloc(width * height * sizeof(double));
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      dsrc[rowidx * width + colidx] = (double)src[rowidx * stride + colidx];
    }
  }
#if USE_WAVE_LIFTING
  ilwt97_2d(dsrc, width, height, width);
#else
  invconvWT_2d(dsrc, width, height, width);
#endif
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      src[rowidx * stride + colidx] = (int)round(dsrc[rowidx * width + colidx]);
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

unsigned encode_coding_params(int width, int height, int stepSize, uchar* bitsOut, unsigned bitPos)
{
  encode_fixlen(width, 12, bitsOut, bitPos);
  encode_fixlen(height, 12, bitsOut, bitPos + 12);
  encode_fixlen(stepSize, 8, bitsOut, bitPos + 24);
  return 32U;
}

unsigned decode_coding_params(uchar* bitsIn, unsigned bitPos, int* width, int* height, int* stepSize)
{
  decode_fixlen(12, bitsIn, bitPos, width);
  decode_fixlen(12, bitsIn, bitPos + 12, height);
  decode_fixlen(8, bitsIn, bitPos + 24, stepSize);
  return 32U;
}

unsigned encode_unit(int* quant, int stride, int partDepth, int* predModes, int* cutModes, uchar* binStream, unsigned bitPos)
{
  unsigned binLen = 0;
  int blkWidth  = MAX_BLOCK_SIZE >> partDepth;
  int blkHeight = MAX_BLOCK_SIZE >> partDepth;
  int blkStride = 1 << partDepth;
  int blkNum    = (1 << partDepth) * (1 << partDepth);
  int blkSize   = MAX_BLOCK_SIZE >> partDepth;

  binLen += encode_fixlen(partDepth, 3, binStream, bitPos + binLen);
  for (int subblk = 0; subblk < blkNum; subblk++)
  {
    binLen += encode_fixlen(predModes[subblk], 2, binStream, bitPos + binLen);
    binLen += encode_fixlen(cutModes[subblk], 2, binStream, bitPos + binLen);
    int row = (subblk / blkStride) * blkHeight;
    int col = (subblk % blkStride) * blkWidth;
    for (int rowidx = row; rowidx < (row + blkHeight); rowidx++)
    {
      for (int colidx = col; colidx < (col + blkWidth); colidx++)
      {
        int currSubband = subband_at_pos(rowidx % blkSize, colidx % blkSize, blkSize, blkSize, stride);
        if (cutModes[subblk] == CUT_HH && currSubband == SUBBAND_HH)
        {
          continue;
        }
        else if (cutModes[subblk] == CUT_LH_HL_HH && currSubband != SUBBAND_LL)
        {
          continue;
        }
        else
        {
          int symbol = *(quant + rowidx * stride + colidx);
          binLen += encode_huffman(symbol, binStream, bitPos + binLen);
        }
      }
    }
  }
  if (binLen % 8 > 0)
  {
    binLen += (8 - (binLen % 8));
  }
  return binLen;
}

unsigned decode_unit(uchar* binStream, unsigned bitPos, int* quant, int stride, int* partDepth, int* predModes, int* cutModes)
{
  unsigned binLen = 0;
  binLen += decode_fixlen(3, binStream, bitPos + binLen, partDepth);

  int blkWidth  = MAX_BLOCK_SIZE >> *partDepth;
  int blkHeight = MAX_BLOCK_SIZE >> *partDepth;
  int blkStride = 1 << *partDepth;
  int blkNum    = (1 << *partDepth) * (1 << *partDepth);
  int blkSize   = MAX_BLOCK_SIZE >> *partDepth;

  for (int subblk = 0; subblk < blkNum; subblk++)
  {
    binLen += decode_fixlen(2, binStream, bitPos + binLen, predModes + subblk);
    binLen += decode_fixlen(2, binStream, bitPos + binLen, cutModes + subblk);
    int row = (subblk / blkStride) * blkHeight;
    int col = (subblk % blkStride) * blkWidth;
    for (int rowidx = row; rowidx < (row + blkHeight); rowidx++)
    {
      for (int colidx = col; colidx < (col + blkWidth); colidx++)
      {
        int currSubband = subband_at_pos(rowidx % blkSize, colidx % blkSize, blkSize, blkSize, stride);
        if (cutModes[subblk] == CUT_HH && currSubband == SUBBAND_HH)
        {
          *(quant + rowidx * stride + colidx) = 0;
        }
        else if (cutModes[subblk] == CUT_LH_HL_HH && currSubband != SUBBAND_LL)
        {
          *(quant + rowidx * stride + colidx) = 0;
        }
        else
        {
          binLen += decode_huffman(binStream, bitPos + binLen, quant + rowidx * stride + colidx);
        }
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

void cut_detail_coefs(int* src, int width, int height, int stride, int cutMode)
{
  if (cutMode == CUT_OFF)
  {
    return;
  }
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int currSubband = subband_at_pos(rowidx, colidx, width, height, stride);
      if (cutMode == CUT_HH && currSubband == SUBBAND_HH)
      {
        src[rowidx * stride + colidx] = 0;
      }
      else if (cutMode == CUT_LH_HL_HH && currSubband != SUBBAND_LL)
      {
        src[rowidx * stride + colidx] = 0;
      }
    }
  }
}

unsigned rd_est_bits(int* x, int width, int height, int stride, int cutMode)
{
  unsigned bits = 4U; //4 bits for coding pred mode und cutting mode
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int currSubband = subband_at_pos(rowidx, colidx, width, height, stride);
      //Huffman coding all relevant subbands
      if (cutMode == CUT_HH && currSubband == SUBBAND_HH)
      {
        continue;
      }
      else if (cutMode == CUT_LH_HL_HH && currSubband != SUBBAND_LL)
      {
        continue;
      }
      else
      {
        int absVal = abs(x[rowidx * stride + colidx]);
        bits += (unsigned)(absVal == 0 ? 1 : absVal + 2);
      }
    }
  }
  return bits;
}

void comp_subblk(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int stride, int bitdepth, int quantSize, int predMode, int cutMode, bool topMargin, bool leftMargin)
{
  int clipMin = -(1 << bitdepth) + 1;
  int clipMax = (1 << bitdepth) - 1;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    memcpy(resi + rowidx * stride, x + rowidx * stride, width * sizeof(int));
  }
  predict(predMode, reco, pred, resi, width, height, stride, leftMargin, topMargin);
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    memcpy(trafo + rowidx * stride, resi + rowidx * stride, width * sizeof(int));
  }
  transform(trafo, width, height, stride);
  clipLR(trafo, width, height, stride, clipMin * TRAFO_SCALE_2D, clipMax * TRAFO_SCALE_2D);
  cut_detail_coefs(trafo, width, height, stride, cutMode);
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    memcpy(quant + rowidx * stride, trafo + rowidx * stride, width * sizeof(int));
  }
  quantize(quant, width, height, stride, quantSize);
}

void reco_subblk(int* quant, int* reco, int width, int height, int stride, int bitdepth, int quantSize, int predMode, bool topMargin, bool leftMargin)
{
  int clipMin = -(1 << bitdepth) + 1;
  int clipMax = (1 << bitdepth) - 1;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    memcpy(reco + rowidx * stride, quant + rowidx * stride, width * sizeof(int));
  }
  dequantize(reco, width, height, stride, quantSize);
  inv_transform(reco, width, height, stride);
  clipLR(reco, width, height, stride, clipMin, clipMax);
  predict(predMode, reco, NULL, NULL, width, height, stride, leftMargin, topMargin);
  clipLR(reco, width, height, stride, 0, clipMax);
}

void rd_search_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int bitdepth, int quantSize, bool topMargin, bool leftMargin, double lambda, int* bestDepth, int* bestPreds, int* bestCuttings)
{
  double bestCost = MAXFLOAT;
  for (int partDepth = 0; partDepth <= ENC_MAX_DEPTH; partDepth++)
  {
    int blkWidth  = MAX_BLOCK_SIZE >> partDepth;
    int blkHeight = MAX_BLOCK_SIZE >> partDepth;
    int blkStride = 1 << partDepth;
    int blkNum    = (1 << partDepth) * (1 << partDepth);
    double sumCosts = 0.0;
    int* bestPredsInLevel    = malloc(blkNum * sizeof(int));
    int* bestCuttingsInLevel = malloc(blkNum * sizeof(int));
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
      int* currOrig = x + offset;
      int* currPred = pred + offset;
      int* currResi = resi + offset;
      int* currTrafo = trafo + offset;
      int* currQuant = quant + offset;
      int* currReco = reco + offset;

      double bestCostSubblk = MAXFLOAT;
      for (int predMode = 0; predMode < NUM_PREDS; predMode++)
      {
        for (int cutMode = 0; cutMode < NUM_CUTTINGS; cutMode++)
        {
#if ENC_SPECIAL_CUT
          if (partDepth == ENC_MAX_DEPTH && cutMode == CUT_LH_HL_HH)
          {
            continue;
          }
#endif
          comp_subblk(currOrig, currPred, currResi, currTrafo, currQuant, currReco, blkWidth, blkHeight, stride, bitdepth, quantSize, predMode, cutMode, hasTop, hasLeft);
          reco_subblk(currQuant, currReco, blkWidth, blkHeight, stride, bitdepth, quantSize, predMode, hasTop, hasLeft);

          //estimate bits per unit and calculate MSE
          unsigned int blkBits = rd_est_bits(currQuant, blkWidth, blkHeight, stride, cutMode);
          unsigned long blkDist = dist_ssd(currOrig, currReco, blkWidth, blkHeight, stride);
          double blkCost = (double)blkDist + lambda * (double)blkBits;
          if (blkCost < bestCostSubblk)
          {
            *(bestPredsInLevel + subblk) = predMode;
            *(bestCuttingsInLevel + subblk) = cutMode;
            bestCostSubblk = blkCost;
          }
        }
      }
      sumCosts += bestCostSubblk;
    }
    if (sumCosts < bestCost)
    {
      bestCost = sumCosts;
      *bestDepth = partDepth;
      memcpy(bestPreds, bestPredsInLevel, blkNum * sizeof(int));
      memcpy(bestCuttings, bestCuttingsInLevel, blkNum * sizeof(int));
    }
    free(bestPredsInLevel);
    free(bestCuttingsInLevel);
  }
}

void comp_reco_unit(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int stride, int bitdepth, int quantSize, int partDepth, int* predModes, int* cutModes, bool topMargin, bool leftMargin)
{
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

    int predMode = predModes[subblk];
    int cutMode = cutModes[subblk];
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
      comp_subblk(currOrig, currPred, currResi, currTrafo, currQuant, currReco, blkWidth, blkHeight, stride, bitdepth, quantSize, predMode, cutMode, hasTop, hasLeft);
    }
    //decompression
    reco_subblk(currQuant, currReco, blkWidth, blkHeight, stride, bitdepth, quantSize, predMode, hasTop, hasLeft);
  }
}

void compress_image(int* x, int* pred, int* resi, int* trafo, int* quant, int* reco, int width, int height, int stepSize, double lambda, uchar* binStream, unsigned* bitPos)
{
  //check minimum trafo block size
  assert(MAX_BLOCK_SIZE >> ENC_MAX_DEPTH >= 16);

  //encode width, height and stepsize
  *bitPos = encode_coding_params(width, height, stepSize, binStream, 0);

  //adjust quantization step-size regarding transform scaling
  int quantSize = stepSize * TRAFO_SCALE_2D;

  int maxNumUnits = (1 << ENC_MAX_DEPTH) * (1 << ENC_MAX_DEPTH);
  int numUnitsX   = width / MAX_BLOCK_SIZE;
  int numUnitsY   = height / MAX_BLOCK_SIZE;
  int numUnits    = numUnitsX * numUnitsY;
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

    int bestDepth     = -1;
    int* bestPreds    = (int*)malloc(maxNumUnits * sizeof(int));
    int* bestCuttings = (int*)malloc(maxNumUnits * sizeof(int));

    rd_search_unit(x + offset, pred + offset, resi + offset, trafo + offset, quant + offset, reco + offset, width, IMG_BITDEPTH, quantSize, topMargin, leftMargin, lambda, &bestDepth, bestPreds, bestCuttings);

    comp_reco_unit(x + offset, pred + offset, resi + offset, trafo + offset, quant + offset, reco + offset, width, IMG_BITDEPTH, quantSize, bestDepth, bestPreds, bestCuttings, topMargin, leftMargin);

    *bitPos += encode_unit(quant + offset, width, bestDepth, bestPreds, bestCuttings, binStream, *bitPos);

    free(bestPreds);
    free(bestCuttings);
  }
  assert(*bitPos % 8 == 0);
}

void reconstruct_image(int* quant, int* reco, int* width, int* height, int* stepSize, uchar* binStream, unsigned* bitPos)
{
  //decode width, height and stepsize
  *bitPos = decode_coding_params(binStream, 0, width, height, stepSize);

  //adjust quantization step-size regarding transform scaling
  int quantSize = *stepSize * TRAFO_SCALE_2D;

  int maxNumUnits = (1 << ENC_MAX_DEPTH) * (1 << ENC_MAX_DEPTH);
  int numUnitsX   = *width / MAX_BLOCK_SIZE;
  int numUnitsY   = *height / MAX_BLOCK_SIZE;
  int numUnits    = numUnitsX * numUnitsY;
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

    int offset = (ui / numUnitsX) * MAX_BLOCK_SIZE * *width + (ui % numUnitsX) * MAX_BLOCK_SIZE;

    int partDepth  = -1;
    int* predModes = (int*)malloc(maxNumUnits * sizeof(int));
    int* cutModes  = (int*)malloc(maxNumUnits * sizeof(int));

    *bitPos += decode_unit(binStream, *bitPos, quant + offset, *width, &partDepth, predModes, cutModes);

    comp_reco_unit(NULL, NULL, NULL, NULL, quant + offset, reco + offset, *width, IMG_BITDEPTH, quantSize, partDepth, predModes, cutModes, topMargin, leftMargin);

    free(predModes);
    free(cutModes);
  }
}