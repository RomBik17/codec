
#include "decoder.h"
using namespace std;

void reconstructedBlockRead(int* reconstructedBlock, int* decodedFrame, int* previousFrame,
    int* modeMatrix, int blockSize, int y, int x)
{
    int wBlock;
    if (((int)(LENGTH / blockSize)) * blockSize == LENGTH) wBlock = LENGTH / blockSize;
    else wBlock = LENGTH / blockSize + 1;
    int mode = modeMatrix[(y / blockSize) * wBlock + (x / blockSize)];

    //by checking mode read the reconstructed block
    int recY;
    int tempX, tempY;
    if (mode == 0)
    {
        if (x >= blockSize)
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                if (y + i >= HEIGHT) tempY = (HEIGHT - 1) * LENGTH;
                else tempY = (y + i) * LENGTH;
                for (int j = 0; j < blockSize; ++j)
                {
                    if (x - blockSize + j >= LENGTH) tempX = LENGTH - 1;
                    else tempX = x - blockSize + j;
                    reconstructedBlock[recY + j] = decodedFrame[tempY + tempX];
                }
            }
        }
        else
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                for (int j = 0; j < blockSize; ++j)
                {
                    reconstructedBlock[recY + j] = 128;
                }
            }
        }
    }
    else if (mode == 1)
    {
        if (y >= blockSize)
        {
            int tempX, tempY;
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                if (y - blockSize + i >= HEIGHT) tempY = (HEIGHT - 1) * LENGTH;
                else tempY = (y - blockSize + i) * LENGTH;
                for (int j = 0; j < blockSize; ++j)
                {
                    if (x + j >= LENGTH) tempX = LENGTH - 1;
                    else tempX = x + j;
                    reconstructedBlock[recY + j] = decodedFrame[tempY + tempX];
                }
            }
        }
        else
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                for (int j = 0; j < blockSize; ++j)
                {
                    reconstructedBlock[recY + j] = 128;
                }
            }
        }
    }
    else if (mode == 2)
    {
        for (int i = 0; i < blockSize; ++i)
        {
            recY = i * blockSize;
            if (y + i >= HEIGHT) tempY = (HEIGHT - 1) * LENGTH;
            else tempY = (y + i) * LENGTH;
            for (int j = 0; j < blockSize; ++j)
            {
                if (x + j >= LENGTH) tempX = LENGTH - 1;
                else tempX = x + j;
                reconstructedBlock[recY + j] = previousFrame[tempY + tempX];
            }
        }
    }
    else if (mode == 3)
    {
        if (x >= blockSize)
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                if (y + i >= HEIGHT) tempY = (HEIGHT - 1) * LENGTH;
                else tempY = (y + i) * LENGTH;
                for (int j = 0; j < blockSize; ++j)
                {
                    if (x - blockSize + j >= LENGTH) tempX = LENGTH - 1;
                    else tempX = x - blockSize + j;
                    reconstructedBlock[recY + j] = previousFrame[tempY + tempX];
                }
            }
        }
        else
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                for (int j = 0; j < blockSize; ++j)
                {
                    reconstructedBlock[recY + j] = 128;
                }
            }
        }
    }
    else if (mode == 4)
    {
        if (x < LENGTH - blockSize)
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                if (y + i >= HEIGHT) tempY = (HEIGHT - 1) * LENGTH;
                else tempY = (y + i) * LENGTH;
                for (int j = 0; j < blockSize; ++j)
                {
                    if (x + blockSize + j >= LENGTH) tempX = LENGTH - 1;
                    else tempX = x + blockSize + j;
                    reconstructedBlock[recY + j] = previousFrame[tempY + tempX];
                }
            }
        }
        else
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                for (int j = 0; j < blockSize; ++j)
                {
                    reconstructedBlock[recY + j] = 128;
                }
            }
        }
    }
    else if (mode == 5)
    {
        if (y >= blockSize)
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                if (y - blockSize + i >= HEIGHT) tempY = (HEIGHT - 1) * LENGTH;
                else tempY = (y - blockSize + i) * LENGTH;
                for (int j = 0; j < blockSize; ++j)
                {
                    if (x + j >= LENGTH) tempX = LENGTH - 1;
                    else tempX = x + j;
                    reconstructedBlock[recY + j] = previousFrame[tempY + tempX];
                }
            }
        }
        else
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                for (int j = 0; j < blockSize; ++j)
                {
                    reconstructedBlock[recY + j] = 128;
                }
            }
        }
    }
    else if (mode == 6)
    {
        if (y < HEIGHT - blockSize)
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                if (y + blockSize + i >= HEIGHT) tempY = (HEIGHT - 1) * LENGTH;
                else tempY = (y + blockSize + i) * LENGTH;
                for (int j = 0; j < blockSize; ++j)
                {
                    if (x + j >= LENGTH) tempX = LENGTH - 1;
                    else tempX = x + j;
                    reconstructedBlock[recY + j] = previousFrame[tempY + tempX];
                }
            }
        }
        else
        {
            for (int i = 0; i < blockSize; ++i)
            {
                recY = i * blockSize;
                for (int j = 0; j < blockSize; ++j)
                {
                    reconstructedBlock[recY + j] = 128;
                }
            }
        }
    }
}

void blockDecoding(int* frameCoef, int* decodedFrame, int* previousFrame, int* modeMatrix, int blockSize, int Q, int y, int x)
{
    int* block = new int[blockSize * blockSize];
    int* reconstructedBlock = new int[blockSize * blockSize];
    int* residualBlock = new int[blockSize * blockSize];

    int wBlock;
    if (((int)(LENGTH / blockSize)) * blockSize == LENGTH) wBlock = LENGTH;
    else wBlock = (LENGTH / blockSize + 1) * blockSize;
    //read coeficients
    int tempX, tempY;
    for (int i = y; i < (blockSize + y); ++i)
    {
        int r = i * wBlock;
        int v = (i - y) * blockSize - x;
        for (int j = x; j < (blockSize + x); ++j)
        {
            residualBlock[v + j] = frameCoef[r + j];
        }
    }

    //read reconstructed block
    reconstructedBlockRead(reconstructedBlock, decodedFrame, previousFrame, modeMatrix, blockSize, y, x);

    //dequantize, retransform and reconstruct block
    dequanting(residualBlock, blockSize, Q);
    IDCT(residualBlock, blockSize);
    reverseSAD(block, reconstructedBlock, residualBlock, blockSize);

    for (int i = y; i < (blockSize + y); ++i)
    {
        for (int j = x; j < (blockSize + x); ++j)
        {
            int r = (i - y) * blockSize - x;
            int h = i * LENGTH;
            if (i < HEIGHT && j < LENGTH)
            {
                //if block is inside of frame then write the data
                if (block[r + j] < 0) block[r + j] = 0;
                else if (block[r + j] > 255) block[r + j] = 255;
                decodedFrame[h + j] = block[r + j];
            }
        }
    }

    delete[] block;
    delete[] residualBlock;
    delete[] reconstructedBlock;
}

//return no-zero data in frame
void split(int* frameCoef, int* decodedFrame, int* previousFrame, int* modeMatrix, int blockSize, int Q)
{

    for (int i = 0; i < HEIGHT; i += blockSize)
    {
        for (int j = 0; j < LENGTH; j += blockSize)
        {
            blockDecoding(frameCoef, decodedFrame, previousFrame, modeMatrix, blockSize, Q, i, j);
        }
    }
}

void decoder(int blockSize, int Q)
{
    const char* in_string = "coef_file.dat";
    const char* out_string = "decoded_y_file.yuv";
    const char* mode_string = "mode_file.dat";

    ifstream mode_file(mode_string, ios::in | ios::binary);
    ofstream decoded_file(out_string, ios::out | ios::binary);
    ifstream in_file(in_string, ios::in | ios::binary);

    int* decodedFrame = new int[pixels_on_video];
    int* previousFrame = new int[pixels_on_video];


    int* modeMatrix;
    int hBlock, wBlock;
    if (((int)(HEIGHT / blockSize)) * blockSize == HEIGHT) hBlock = HEIGHT / blockSize;
    else hBlock = HEIGHT / blockSize + 1;

    if (((int)(LENGTH / blockSize)) * blockSize == LENGTH) wBlock = LENGTH / blockSize;
    else wBlock = LENGTH / blockSize + 1;

    modeMatrix = new int[hBlock * wBlock];

    int* frameCoef;
    frameCoef = new int[hBlock * blockSize * wBlock * blockSize];

    while (true)
    {
        //read prediction modes
        for (int i = 0; i < hBlock * wBlock; ++i)
        {
            mode_file >> modeMatrix[i];
        }

        if (mode_file.eof()) break;

        //read coeficients
        for (int i = 0; i < hBlock * blockSize * wBlock * blockSize; ++i)
        {
            in_file >> frameCoef[i];
        }

        split(frameCoef, decodedFrame, previousFrame, modeMatrix, blockSize, Q);

        //write decoded frame
        for (int i = 0; i < pixels_on_video; ++i)
        {
            previousFrame[i] = decodedFrame[i];
            char a = (unsigned char)decodedFrame[i];
            decoded_file.write(&a, 1);
        }
    }

    delete[] modeMatrix;
    delete[] decodedFrame;
    delete[] previousFrame;
    delete[] frameCoef;

    mode_file.close();
    decoded_file.close();
    in_file.close();
}
