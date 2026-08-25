typedef unsigned char uchar;

//fix-length coding
unsigned encode_fixlen(int n, int len, uchar* bitsOut, unsigned bitPos);
unsigned decode_fixlen(int len, uchar* bitsIn, unsigned bitPos, int* n);
//Huffman coding
unsigned encode_huffman(int n, uchar* bitsOut, unsigned bitPos);
unsigned decode_huffman(uchar* bitsIn, unsigned bitPos, int* n);