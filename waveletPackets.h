#pragma once

#include <vector>
#include <map>
#include "wavelets.h"
#include "wt.h"
#include "common.h"
#include <iostream>
#include <algorithm>

//小波包分解参考PyWavelets库中小波包分解文件，基本按照文件中Python代码移植过来

template<class node_type>
class BaseNode;

class Node;

class Node2D;

class NodeND;


inline std::vector<std::string> getGraycodeOrder(short level)
{
    std::string x = "a";
    std::string y = "d";
    std::vector<std::string> grayCodeOrder = {x, y};
    std::vector<std::string> grayCodeOrderReverse = {y, x};
    for (int i = 0; i < level - 1; i++)
    {
        std::vector<std::string> temp;
        for (const std::string& j: grayCodeOrder)
        {
            temp.emplace_back(x + j);
        }

        std::vector<std::string> temp2;
        for (const std::string& j: grayCodeOrderReverse)
        {
            temp2.emplace_back(y + j);
        }
        grayCodeOrder.clear();
        grayCodeOrder.insert(grayCodeOrder.end(), temp.begin(), temp.end());
        grayCodeOrder.insert(grayCodeOrder.end(), temp2.begin(), temp2.end());
        grayCodeOrderReverse.assign(grayCodeOrder.begin(), grayCodeOrder.end());
        std::reverse(grayCodeOrderReverse.begin(), grayCodeOrderReverse.end());
    }

    return grayCodeOrder;
}

// 小波分解节点类
template<class node_type>
class BaseNode
{
public:
    BaseNode<node_type> *parent; // 父节点
    DiscreteWavelet *wavelet;
    double *data = nullptr;
    MODE mode;  // 分解方式
    short level = 0;  // 当前层级
    std::string path; // 节点路径
    short maxLevel = 0;  // 最大层级
    size_t data_len = 0;  // 数据长度
    size_t filter_len = 0;  // 小波长度

    // 不赋值
    std::vector<std::string> PARTS;
    short PART_LEN = 1;
    std::map<std::string, node_type *> subNodeList; // 子节点map

    node_type *derivedNode = nullptr;

    BaseNode(BaseNode<node_type> *parent, double *data, size_t data_len, const std::string &name,
             const std::vector<std::string> &subName) :
            parent(parent), data(data), data_len(data_len), PARTS(subName)
    {
        if (parent != nullptr)
        {
            level = parent->level + 1;
            path = parent->path + name;
            maxLevel = parent->maxLevel;
            mode = parent->mode;
            wavelet = parent->wavelet;
            filter_len = wavelet->dec_len;
        }

        for (const std::string &i: subName)
        {
            setNode(nullptr, i);
        }

        derivedNode = static_cast<node_type * >(this);
    }

    ~BaseNode()
    {

        for (auto iter = subNodeList.rbegin(); iter != subNodeList.rend(); iter++)
        {
            if (iter->second != nullptr)
                delete iter->second;
        }
        if (data != nullptr && this->level != 0)
        {
            delete[] data;
        }
    }

    double *reconstruct()
    {
        if (subNodeList["a"] != nullptr || subNodeList["d"] != nullptr)
        {
            if (this->level == 0)
            {
                auto t = derivedNode->_reconstruct();
                this->data = t;
                return t;
            }
            return derivedNode->_reconstruct();
        }
        return data;
    }

    void decompose()
    {
        if (this->level < this->maxLevel)
        {
            derivedNode->_decompose();
        }
        else
        {
            /*throw "达到最大分解层数";*/
        }
    }

    [[nodiscard]] BaseNode<node_type> getSubNode(std::string name) const
    {
        return subNodeList[name];
    }

    [[nodiscard]] std::string nodeName() const
    {
        return path.substr(path.size() - PART_LEN, 1);
    }

    [[nodiscard]] std::string nodeFullName() const
    {
        return path;
    }

    [[nodiscard]] int getLevel() const
    {
        return level;
    }

    node_type *getNode(std::string name)
    {
        return subNodeList[name];
    }

    void setNode(node_type *node, std::string name)
    {
        subNodeList.insert(std::pair<std::string, node_type *>(name, node));
    }

    short getMaxLevel()
    {
        if (maxLevel != 0)
        {
            return maxLevel;
        }
        else if (parent != nullptr)
        {
            maxLevel = parent->getMaxLevel();
        }
        else if (data != nullptr)
        {
            maxLevel = level + dwt_max_level(data_len, filter_len);
        }

        return maxLevel;
    }


    [[nodiscard]] auto *getData() const
    {
        return data;
    };

    node_type *getitem(const std::string& name)
    {

    }

    void _setitem(node_type *nodePointer, std::string name, double *new_data)
    {
        if (nodePointer->level < name.size())
        {
            int pos = nodePointer->level;
            std::string newNodeName = name.substr(pos, 1);
            if (nodePointer->getNode(newNodeName) == nullptr)
            {
                derivedNode->createSubnode(nodePointer, newNodeName, nullptr);
            }
            _setitem(nodePointer->subNodeList[newNodeName], name, new_data);
        }
        else
        {
            nodePointer->data = new_data;
        }
    }

    void setitem(std::string name, double *new_data)
    {
        std::string newNodeName = name.substr(level, 1);
        if (getNode(newNodeName) == nullptr)
        {
            int output = dwt_buffer_length(data_len, wavelet->dec_len, mode);
            auto *newNode = new node_type(this, nullptr, output, newNodeName, PARTS);
            subNodeList[newNodeName] = newNode;
        }
        _setitem(this->subNodeList[newNodeName], name, new_data);
    };
};


class Node : public BaseNode<Node>
{
public:
    std::string A = "a";
    std::string D = "d";
    std::vector<std::string> subName;

    Node(BaseNode<Node> *parent, double *data, size_t data_len, std::string name, std::vector<std::string> subName) :
            BaseNode<Node>(parent, data, data_len, name, subName), subName(subName) {};

    // 分解输入信号
    void walk(std::vector<Node *> *nodeList, bool is_decompose, short nowLevel)
    {
        if (level <= maxLevel && level != nowLevel)
        {
            for (const std::string& name: PARTS)
            {
                Node *subNode = this->getNode(name);
                if (is_decompose && data != nullptr)
                {
                    decompose();
                    subNode = getNode(name);
                    if (subNode != nullptr)
                    {
                        subNode->walk(nodeList, is_decompose, nowLevel);
                    }
                }
            }
        }
        else
        {
            nodeList->push_back(this);
        }
    };

    void _decompose()
    {
        if (data != nullptr)
        {

            auto output = dwt_buffer_length(data_len, wavelet->dec_len, mode);
            auto *data_a = new double[output]();
            auto *data_d = new double[output]();
            int dec_a_res = dec_a(data, data_len, wavelet, data_a, output, mode); // 低通
            int dec_d_res = dec_d(data, data_len, wavelet, data_d, output, mode); // 高通
            if (dec_a_res && dec_d_res)
            {
                /*throw "小波分解错误";*/
            }
            auto *subNodeA = new Node(this, data_a, output, A, subName);
            auto *subNodeD = new Node(this, data_d, output, D, subName);
            subNodeList[A] = subNodeA;
            subNodeList[D] = subNodeD;
        }
    }

    double *_reconstruct()
    {
        auto subNodeA = getNode(A);
        auto subNodaD = getNode(D);
        double *data_a = nullptr;
        double *data_d = nullptr;
        int input_len;
        if (subNodeA != nullptr)
        {
            data_a = subNodeA->reconstruct();
            input_len = subNodeA->data_len;
        }
        if (subNodaD != nullptr)
        {
            data_d = subNodaD->reconstruct();
            input_len = subNodaD->data_len;
        }

        if (subNodeA == nullptr && subNodaD == nullptr)
        {
            /*throw "此节点为叶节点不能重构";*/
        }
        else
        {
            auto output_len = idwt_buffer_length(input_len, wavelet->rec_len, mode);
            auto *output = new double[output_len]();
            int res = idwt(data_a, input_len, data_d, input_len, output, output_len, this->wavelet,
                           this->mode); // 小波逆变换
            return output;
        }
    };

    static void createSubnode(Node *nodePointer, const std::string &name, double *data)
    {
        auto output = dwt_buffer_length(nodePointer->data_len, nodePointer->wavelet->dec_len, nodePointer->mode);
        auto *newNode = new Node(nodePointer, data, output, name, nodePointer->subName);
        nodePointer->subNodeList[name] = newNode;
    }

};

// 小波包分解
class WaveletPacket : public Node
{
public:
    WaveletPacket(double *data, size_t data_len, WAVELET_NAME wavelet, MODE mode, short maxLevel, unsigned int order) :
            Node(nullptr, data, data_len, "", {"a", "d"})
    {

        this->wavelet = discrete_wavelet(wavelet, order);
        this->mode = mode;
        if (data != nullptr)
        {
            this->data_len = data_len;
            if (maxLevel == 0)
            {
                this->maxLevel = dwt_max_level(data_len, this->wavelet->dec_len);
            }
            else
            {
                this->maxLevel = maxLevel;
            }
        }
    }

    // 初始化后使用此函数获取每个节点的数据
    std::vector<Node *> getLevel(short level, const std::string &order)
    {
        std::vector<Node *> nodeList;
        walk(&nodeList, true, level);
        if (order == "natural")
        {
            return nodeList;
        }
        else if (order == "freq")
        {
            std::map<std::string, Node *> name2node;
            for (auto tmp: nodeList)
            {
                name2node.insert(std::pair<std::string, Node *>(tmp->nodeFullName(), tmp));
            }
            std::vector<std::string> graycodeorder = getGraycodeOrder(level);

            nodeList.clear();
            for (const std::string &path: graycodeorder)
            {
                if (name2node.count(path))
                {
                    nodeList.push_back(name2node[path]);
                }
            }
            return nodeList;
        }
    }
};


