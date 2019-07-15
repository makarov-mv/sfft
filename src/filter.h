#pragma once

#include <vector>
#include <unordered_map>
#include "tree.h"


class Filter {
public:
    using Node = SplittingTree::Node;
    using NodePtr = SplittingTree::NodePtr;

    Filter(const NodePtr& node, int64_t size) {
        path_ = node->GetRootPath();
        auto label = node->label;
        label_ = label;

        filter_[0] = 1;
        for (int i = 0; i < static_cast<int>(path_.size()); ++i) {
            std::unordered_map<int64_t, complex_t> updated_filter;
            int64_t shift = size >> (path_[i] + 1);
            auto phase = CalcKernel(-label, 1 << (path_[i] + 1));
            for (auto it : filter_) {
                int64_t index = it.first;
                updated_filter[index] = (FilterValueAtTime(index) + phase * FilterValueAtTime((index + shift) % size)) / 2.;
                updated_filter[(index + size - shift) % size] = (FilterValueAtTime((index + size - shift) % size) + phase * FilterValueAtTime(index)) / 2.;
            }
            filter_.swap(updated_filter);
            updated_filter.clear();
        }
    }

    complex_t FilterValueAtTime(int64_t time) const {
        auto it = filter_.find(time);
        if (it != filter_.end()) {
            return it->second;
        }
        return 0.;
    }

    const std::unordered_map<int64_t, complex_t>& FilterTime() const {
        return filter_;
    }

    complex_t FilterFrequency(int psi) const {
        double g = PI * (psi - label_);
        complex_t freq =  1;
        int j = 0;
        for (int level : path_) {
            for (; j < level; ++j) {
                g /= 2;
            }
            freq *= CalcFrequencyFactor(g);
        }
        return freq;
    }

private:
    complex_t CalcFrequencyFactor(double g) const {
        return {(1 + cos(g)) / 2, sin(g) / 2};
    }

    int64_t label_;
    std::vector<int> path_;
    std::unordered_map<int64_t, complex_t> filter_;
};
