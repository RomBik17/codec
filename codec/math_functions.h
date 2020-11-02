
#include "video_constants.h"

#define CONST_BITS 13
#define ONE	((int) 1)
#define CONST_SCALE (ONE << CONST_BITS)
#define FIX(x)	((int) ((x) * CONST_SCALE))
#define DESCALE(x) ((int)((x) / CONST_SCALE))

double PSNR(char* startFrame, char* dumpedFrame);
void quanting(int* block, int blockSize, int Q);
void dequanting(int* block, int blockSize, int Q);
void IDCT(int* block, int blockSize);
void DCT(int* block, int blockSize);
void SAD(int* currentBlock, int* predictedBlock, int* sadBlock, int blockSize);
void reverseSAD(int* currentBlock, int* predictedBlock, int* sadBlock, int blockSize);
