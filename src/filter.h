#pragma once

#include <vector>
#include "disfft.h"
#include "math.h"
#include <unordered_map>


class Filter {
public:
    using Node = SplittingTree::Node;
    using NodePtr = SplittingTree::NodePtr;


    std::vector<complex_t> time;
    std::vector<complex_t> frequency;
    ssize_t size;

    Filter(NodePtr node, ssize_t size_) : size(size_) {
        std::vector<bool> path_mask = node->GetRootPathMask(); // optimize
        phase.resize(static_cast<size_t>(log2(size)), 0);
        uint64_t label = node->label;

        for (int i = 0; i + 1 < path_mask.size(); ++i) {
            int j = path_mask.size() - i - 1;
            if (path_mask[j]) {
                phase[i] = std::exp(-2 * PI * I * label / (1 << (i + 1)));
            }
        }

        G[0] = 1;
        for (int i = 0; i < phase.size(); ++i) {
            if (phase[i] != 0) {
                std::unordered_map<int, complex_t> updated_G;
                for (auto it : G) {
                    int index = it.first;
                    updated_G[index] = G[index] / 2 + phase[i] * G[(index + (size >> (i + 1))) % size] / 2;
                    updated_G[(index + size - (size >> (i + 1))) % size] = G[(index + size - (size >> (i + 1))) % size] / 2 + phase[i] * G[index] / 2;
                }
                G.swap(updated_G);
                updated_G.clear();
            }
        }
    }

    complex_t FilterTime(int psi) const {
        return G[psi];
    }

    complex_t FilterFrequency(int psi) const {
        complex_t freq =  1;
        for (size_t i = 0; i < phase.size(); ++i) {
            if (phase[i] != 0) {
                freq *= (1 + phase[i] * std::exp(2 * PI * I * psi / (1 << (i + 1))));
            }
        }
        return freq;
    }

private:

    std::vector<complex_t> phase;
    std::unordered_map<int, complex_t> G;
};
