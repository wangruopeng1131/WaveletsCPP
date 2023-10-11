# WaveletsCPP
The C++ version of WaveletPacket.

## Description
This code is translated from Pywavelets.

## Usage
Decompose ECG.
```c++
#include "waveletPackets.h"
#include "wavelets.h"

short maxLevel = 3;
float freqBand = float(256) / float(2 << (maxLevel - 1));
MODE mode = MODE_SYMMETRIC;
int order = 4;
WAVELET_NAME wavelet = SYM;
int data_len = int(data.size());
double *data_ptr = data.data();
WaveletPacket wp(data_ptr, data_len, wavelet, mode, maxLevel, order);
vector<Node*> freqTree = wp.getLevel(maxLevel, "freq");
```

Reconstruct ECG from different frequency band.
```c++
double **reconstruct = new double *[iterFreqs.size()];
std::vector<WaveletPacket> newList;
for (int i = 0; i < iterFreqs.size(); ++i)
{
    WaveletPacket new_wp(nullptr, data_len, wavelet, mode, maxLevel, order);
    newList.push_back(new_wp);
}

// 小波包重构
for (int i = 0; i < iterFreqs.size(); ++i)
{
    std::map<std::string, int> iter = iterFreqs[i];
    int fmin = iter["fmin"];
    int fmax = iter["fmax"];

    for (int j = 0; j < freqTree.size(); ++j)
    {
        float bandMin = float(j) * freqBand;
        float bandMax = bandMin + freqBand;

        if (fmin <= bandMin && fmax >= bandMax)
        {
            newList[i].setitem(freqTree[j]->nodeFullName(), freqTree[j]->getData());
        }
    }
    reconstruct[i] = newList[i].reconstruct();
}
```

## Build
```
mkdir build
cd build 
cmake .. -DCMAKE_BUILD_TYPE=Release
```

## TODO
Wavelet 2D
Wavelet 3D