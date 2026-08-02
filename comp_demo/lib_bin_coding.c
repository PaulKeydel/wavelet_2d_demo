#include <stdlib.h>
#include "lib_bin_coding.h"


unsigned encode_fixlen(int n, int len, uchar* bitsOut, unsigned bitPos)
{
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

unsigned decode_fixlen(int len, uchar* bitsIn, unsigned bitPos, int* n)
{
  *n = 0;
  for (int i = 0; i < len; i++)
  {
    *n += checkBit(bitsIn, (bitPos + len - 1 - i)) * (1 << i);
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
      setBit(bitsOut, (bitPos + i), 1);
    }
    setBit(bitsOut, (bitPos + bitlen - 2), 0);
    setBit(bitsOut, (bitPos + bitlen - 1), (n < 0 ? 0 : 1));
  }
  return bitlen;
}

unsigned decode_huffman(uchar* bitsIn, unsigned bitPos, int* n)
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