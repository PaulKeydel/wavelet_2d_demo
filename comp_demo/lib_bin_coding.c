#include <stdlib.h>
#include "lib_bin_coding.h"


//macros for bit manipulation on unsigned char array
#define setBit(arr,k,val) ( arr[(k)/8] = (val == 1) ? arr[(k)/8] | (1 << (7 - (k)%8)) : arr[(k)/8] & ~(1 << (7 - (k)%8)) )
#define checkBit(arr,k)   ( (arr[(k)/8] >> (7 - (k)%8)) & 1 )


unsigned encode_fixlen(int n, int len, uchar* bitsOut, unsigned bitPos)
{
  for (int i = 0; i < len; i++)
  {
    setBit(bitsOut, bitPos + i, 0);
  }
  int i = 0;
  while (n > 0)
  {
    setBit(bitsOut, bitPos + len - i - 1, (n % 2));
    n = n / 2;
    i++;
  }
  return len;
}

unsigned decode_fixlen(int len, uchar* bitsIn, unsigned bitPos, int* n)
{
  *n = 0;
  for (int i = 0; i < len; i++)
  {
    *n += checkBit(bitsIn, bitPos + len - 1 - i) * (1 << i);
  }
  return len;
}

unsigned encode_huffman(int n, uchar* bitsOut, unsigned bitPos)
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
      setBit(bitsOut, bitPos + i, 1);
    }
    setBit(bitsOut, bitPos + bitlen - 2, 0);
    setBit(bitsOut, bitPos + bitlen - 1, (n < 0 ? 0 : 1));
  }
  return bitlen;
}

unsigned decode_huffman(uchar* bitsIn, unsigned bitPos, int* n)
{
  *n = 0;
  while (checkBit(bitsIn, bitPos + *n) != 0)
  {
    *n = *n + 1;
  }
  unsigned bitlen = (*n == 0 ? 1U : *n + 2U);
  if (*n > 0)
  {
    *n *= (checkBit(bitsIn, bitPos + bitlen - 1) == 0 ? -1 : 1);
  }
  return bitlen;
}

unsigned encode_blk_huffman(int* src, int width, int height, int stride, uchar* bitsOut, unsigned bitPos)
{
  unsigned binLen = 0U;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int symbol = *(src + rowidx * stride + colidx);
      binLen += encode_huffman(symbol, bitsOut, bitPos + binLen);
    }
  }
  return binLen;
}

unsigned decode_blk_huffman(uchar* bitsIn, unsigned bitPos, int* src, int width, int height, int stride)
{
  unsigned binLen = 0U;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      binLen += decode_huffman(bitsIn, bitPos + binLen, src + rowidx * stride + colidx);
    }
  }
  return binLen;
}