//typedef for binary stream
typedef unsigned char uchar;

//macros for bit manipulation on unsigned char array
#define setBit(arr,k,val) ( arr[(k/8)] = (val == 1) ? arr[(k/8)] | (1 << (7 - (k%8))) : arr[(k/8)] & ~(1 << (7 - (k%8))) )
#define checkBit(arr,k)   ( (arr[(k/8)] >> (7 - (k%8))) & 1 )

//fix-length coding
unsigned encode_fixlen(int n, int len, uchar* bitsOut, unsigned bitPos);
unsigned decode_fixlen(int len, uchar* bitsIn, unsigned bitPos, int* n);
//Huffman coding
unsigned encode_huffman(int n, uchar* bitsOut, unsigned bitPos);
unsigned decode_huffman(uchar* bitsIn, unsigned bitPos, int* n);