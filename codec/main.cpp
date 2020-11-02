
#include "encoder.h"
#include "decoder.h"
using namespace std;


int main()
{
    cout << "Choose:\n    1) encode\n    2) decode" << endl;
    int variant;
    cin >> variant;

    cout << "Input the size of block from 2 to 64:" << endl;
    int blockSize;
    cin >> blockSize;
    cout << "Input the size of Q from 2 to 64:" << endl;
    int Q;
    cin >> Q;

    switch (variant)
    {
    case 1:
    {
        encoder(blockSize, Q);
        break;
    }
    case 2:
    {
        decoder(blockSize, Q);
        break;
    }
    default:
        break;
    }

    cout << "Finished";
    return 0;
}