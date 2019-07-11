#pragma once

#include <vector>
#include <unordered_map>
#include "tree.h"


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
                phase_[i] = CalcKernel(-label, 1 << (i + 1));
            }
        }

        filter_[0] = 1;
        for (int i = 0; i < phase_.size(); ++i) {
            if (NonZero(phase_[i])) {
                std::unordered_map<int, complex_t> updated_filter;
                for (auto it : filter_) {
                    int index = it.first;
                    int64_t shift = size >> (i + 1);
                    updated_filter[index] = (filter_[index] + phase_[i] * filter_[(index + shift) % size]) / 2.;
                    updated_filter[(index + size - shift) % size] = (filter_[(index + size - shift) % size] + phase_[i] * filter_[index]) / 2.;
                }
                filter_.swap(updated_filter);
                updated_filter.clear();
            }
        }
    }

    const std::unordered_map<int, complex_t>& FilterTime() const {
        return filter_;
    }

    complex_t FilterFrequency(int psi) const {
        complex_t freq =  1;
        for (size_t i = 0; i < phase_.size(); ++i) {
            if (NonZero(phase_[i])) {
                freq *= (1. + phase_[i] * CalcKernel(psi, 1 << (i + 1))) / 2.;
            }
        }
        return freq;
    }

private:

    int64_t size_;
    std::vector<complex_t> phase_;
    std::unordered_map<int, complex_t> filter_;
};
