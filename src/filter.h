#pragma once

#include <vector>
#include <unordered_map>
#include <map>
#include <iterator>
#include "tree.h"
#include "signal.h"
#include "pair_allocator.h"
#include <memory_resource>

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
    using Allocator = StackAllocator<std::pair<const Key, complex_t>>;
    using Hashmap = std::unordered_map<Key, complex_t, std::hash<Key>, std::equal_to<>, Allocator>;

    Filter(const SplittingTree& tree, const NodePtr& node, const SignalInfo& info) : label_(info, node->label), info_(info) {
        path_ = tree.GetRootPath(node);
        static Allocator main;
        static Allocator secondary;
        main.consolidate();
        secondary.consolidate();
        Allocator* main_ptr = &main;
        Allocator* secondary_ptr = &secondary;
        Hashmap filter_base(main);
        Hashmap updated_filter_base(secondary);
        Hashmap* filter_ptr = &filter_base;
        Hashmap* upd_filter_ptr = &updated_filter_base;
        filter_base[Key(info, 0)] = 1;
        for (int i = 0; i < static_cast<int>(path_.size()); ++i) {
            Hashmap& filter = *filter_ptr;
            Hashmap& updated_filter = *upd_filter_ptr;
            int current_period = CalcCurrentPeriod(i);
            int64_t subtree_level = CalcSubtreeLevel(i);
            int64_t shift = info.SignalWidth(current_period) >> subtree_level;
            assert(current_period >= 0 && current_period < info.Dimensions());
            
            int64_t phase_shift = label_[current_period];
            auto phase = CalcKernel(-phase_shift, 1 << subtree_level);
            phase_shift_.push_back(phase_shift);
            
            current_period_.push_back(current_period);
            SubtreeLevel_.push_back(subtree_level);
            for (auto it : filter) {
                Key index = it.first;
                updated_filter[index] = (MapFilterValueAtTime(filter, index) + phase * MapFilterValueAtTime(filter, index.IncreaseAt(current_period, shift))) / 2.;
                updated_filter[index.IncreaseAt(current_period, -shift)] = (MapFilterValueAtTime(filter, index.IncreaseAt(current_period, -shift)) + phase * MapFilterValueAtTime(filter, index)) / 2.;
            }
            std::swap(filter_ptr, upd_filter_ptr);
            std::swap(main_ptr, secondary_ptr);
            upd_filter_ptr->~Hashmap();
            secondary_ptr->consolidate();
            new (upd_filter_ptr) Hashmap(*secondary_ptr);
        }
        Hashmap& filter = *filter_ptr;
        filter_.reserve(filter.size());
        filter_.insert(filter_.end(), std::make_move_iterator(filter.begin()), std::make_move_iterator(filter.end()));        
    }

    const std::vector<std::pair<Key, complex_t>>& FilterTime() const {
        return filter_;
    }

    complex_t FilterFrequency(const Key& key) const {
        complex_t freq = 1.;
        for (size_t i = 0; i < path_.size(); ++i) {
            //int current_period = CalcCurrentPeriod(i);
            int64_t numer = key[current_period_[i]] - phase_shift_[i] + (1 << SubtreeLevel_[i]);
            int64_t denom = 1 << SubtreeLevel_[i];
            if ((numer & (denom-1)) == (1 << (SubtreeLevel_[i]-1))){
                return 0;
            } else {
                freq *= (1. + CalcKernel(numer, denom)) / 2.;
            }
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
        int cur_level = path_[path_pos];
        for (int i = 0; i < info_.Dimensions(); ++i) {
            if (cur_level >= info_.LogSignalWidth(i)) {
                cur_level -= info_.LogSignalWidth(i);
            } else {
                break;
            }
        }
        return cur_level + 1;
    }

    int CalcCurrentPeriod(int path_pos) const {
        int cur_level = path_[path_pos];
        int cur_period = 0;
        while (cur_period < info_.Dimensions() && cur_level >= info_.LogSignalWidth(cur_period)) {
            cur_level -= info_.LogSignalWidth(cur_period);
            ++cur_period;
        }
        return cur_period;
    }

    complex_t MapFilterValueAtTime(const Hashmap& filter, const Key& time) const {
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
    std::vector<int> path_;
    std::vector<int64_t> phase_shift_;
    std::vector<int> current_period_;
    std::vector<int64_t> SubtreeLevel_;
    std::vector<std::pair<Key, complex_t>> filter_;
};
