#pragma once

#include <vector>
#include <unordered_map>
#include <iterator>
#include "tree.h"
#include "signal.h"

namespace std {
    template<>
    struct hash<Key> {
    public:
        std::size_t operator()(const Key& key) const {
            return hasher_(key.Flatten());
        }

    private:
        std::hash<int64_t> hasher_;
    };
}

using FrequencyMap = std::unordered_map<Key, complex_t>;
using NodePtr = SplittingTree::NodePtr;
using Node = SplittingTree::Node;

FrequencyMap MapUnion(const FrequencyMap& a, const FrequencyMap& b) {
    FrequencyMap c;
    for (const auto& w : a) {
        c[w.first] += w.second;
    }
    for (const auto& w : b) {
        c[w.first] += w.second;
    }
    return c;
}

class Filter {
public:

    Filter(const SplittingTree& tree, const NodePtr& node, const SignalInfo& info) : label_(info, node->label), info_(info), period_size_(CalcLog(info.SignalWidth())) {
        path_ = tree.GetRootPath(node);
        std::unordered_map<Key, complex_t> filter;
        filter[Key(info, 0)] = 1;
        for (int i = 0; i < static_cast<int>(path_.size()); ++i) {
            std::unordered_map<Key, complex_t> updated_filter;
            int current_period = CalcCurrentPeriod(i);
            int64_t subtree_level = CalcSubtreeLevel(i);
            int64_t shift = info.SignalWidth() >> subtree_level;
            assert(current_period >= 0 && current_period < info.Dimensions());
            auto phase = CalcKernel(-label_[current_period], 1 << subtree_level);
            phase_.push_back(phase);
            for (auto it : filter) {
                Key index = it.first;
                updated_filter[index] = (MapFilterValueAtTime(filter, index) + phase * MapFilterValueAtTime(filter, index.IncreaseAt(current_period, shift))) / 2.;
                updated_filter[index.IncreaseAt(current_period, -shift)] = (MapFilterValueAtTime(filter, index.IncreaseAt(current_period, -shift)) + phase * MapFilterValueAtTime(filter, index)) / 2.;
            }
            filter.swap(updated_filter);
            updated_filter.clear();
        }
        filter_.reserve(filter.size());
        filter_.insert(filter_.end(), std::make_move_iterator(filter.begin()), std::make_move_iterator(filter.end()));
    }

    const std::vector<std::pair<Key, complex_t>>& FilterTime() const {
        return filter_;
    }

    complex_t FilterFrequency(const Key& key) const {
        complex_t freq = 1.;
        for (size_t i = 0; i < path_.size(); ++i) {
            int current_period = CalcCurrentPeriod(i);
            freq *= (1. + phase_[i] * CalcKernel(key[current_period], 1 << CalcSubtreeLevel(i))) / 2.;
        }
        return freq;
    }

    complex_t FilterValueAtTime(const Key& time) const {
        for (const auto& freq : filter_) {
            if (freq.first == time) {
                return freq.second;
            }
        }
        return 0.;
    }


private:
    int CalcSubtreeLevel(int path_pos) const {
        return path_[path_pos] - CalcCurrentPeriod(path_pos) * period_size_ + 1;
    }

    int CalcCurrentPeriod(int path_pos) const {
        return path_[path_pos] / period_size_;
    }

    complex_t MapFilterValueAtTime(const std::unordered_map<Key, complex_t>& filter, const Key& time) const {
        auto it = filter.find(time);
        if (it != filter.end()) {
            return it->second;
        }
        return 0.;
    }

// TODO: Calc frequencies faster
//    complex_t CalcFrequencyFactor(double g) const {
//        return {(1 + cos(g)) / 2, sin(g) / 2};
//    }

    Key label_;
    SignalInfo info_;
    int period_size_;
    std::vector<int> path_;
    std::vector<complex_t> phase_;
    std::vector<std::pair<Key, complex_t>> filter_;
};