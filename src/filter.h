#pragma once

#include <vector>
#include "tree.h"
#include "math.h"
#include <unordered_map>

#define _USE_MATH_DEFINES
const double PI = M_PI;
const complex_t I = complex_t(0, 1);


class Filter {
public:
    using Node = SplittingTree::Node;
    using NodePtr = SplittingTree::NodePtr;

    Filter(const NodePtr& node, int64_t size) : size_(size) {
        std::vector<bool> path_mask = node->GetRootPathMask(); // optimize
        phase_.resize(static_cast<size_t>(log2(size)), 0);
        uint64_t label = node->label;

        for (int i = 0; i + 1 < path_mask.size(); ++i) {
            int j = path_mask.size() - i - 1;
            if (path_mask[j]) {
                phase_[i] = std::exp(-2 * PI * I * label / (1 << (i + 1)));
            }
        }

        filter_[0] = 1;
        for (int i = 0; i < phase_.size(); ++i) {
            if (phase_[i] != 0) {
                std::unordered_map<int, complex_t> updated_G;
                for (auto it : filter_) {
                    int index = it.first;
                    updated_G[index] = filter_[index] / 2 + phase_[i] * filter_[(index + (size >> (i + 1))) % size] / 2;
                    updated_G[(index + size - (size >> (i + 1))) % size] = filter_[(index + size - (size >> (i + 1))) % size] / 2 + phase_[i] * filter_[index] / 2;
                }
                filter_.swap(updated_G);
                updated_G.clear();
            }
        }
    }

    complex_t FilterTime(int psi) const {
        return filter_[psi];
    }

    complex_t FilterFrequency(int psi) const {
        complex_t freq =  1;
        for (size_t i = 0; i < phase_.size(); ++i) {
            if (phase_[i] != 0) {
                freq *= (1 + phase_[i] * std::exp(2 * PI * I * psi / (1 << (i + 1))));
            }
        }
        return freq;
    }

private:

    int64_t size_;
    std::vector<complex_t> phase_;
    std::unordered_map<int, complex_t> filter_;
};
