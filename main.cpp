//
// Created by wangr on 2023/10/8.
//

#include "waveletPackets.h"
#include "wavelets.h"
#include <fstream>

int main()
{
    std::ifstream d(R"(F:\project\wavelet\ecg.txt)");
    std::vector<double> data;
    std::string lines;
    std::size_t pos{};
    while (std::getline(d, lines))
    {
        double temp = std::stod(lines, &pos);
        data.emplace_back(temp);
    }
    d.close();

    short maxLevel = 3;
    float freqBand = float(256) / float(2 << (maxLevel - 1));
    MODE mode = MODE_SYMMETRIC;
    int order = 4;
    WAVELET_NAME wavelet = SYM;
    int data_len = int(data.size());
    double *data_ptr = data.data();
    std::vector<std::map<std::string, int>> iterFreqs = {{std::pair<std::string, int>("fmin", 1),
                                                                 std::pair<std::string, int>("fmax", 5)},
                                                         {std::pair<std::string, int>("fmin", 5),
                                                                 std::pair<std::string, int>("fmax", 20)},
                                                         {std::pair<std::string, int>("fmin", 20),
                                                                 std::pair<std::string, int>("fmax", 50)}};

    WaveletPacket wp(data_ptr, data_len, wavelet, mode, maxLevel, order);
    auto freqTree = wp.getLevel(maxLevel, "freq");

    // 小波包重构变量声明
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
            // 第i个频段的最小频率
            float bandMin = float(j) * freqBand;
            // 第i个频段的最大频率
            float bandMax = bandMin + freqBand;
            // 判断第i个频段是否在要分析的范围内
            if (fmin <= bandMin && fmax >= bandMax)
            {
                newList[i].setitem(freqTree[j]->nodeFullName(), freqTree[j]->getData());
            }
        }
        reconstruct[i] = newList[i].reconstruct();
    }
    return 0;
}