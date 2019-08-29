#pragma once

#include <vector>
#include <unordered_map>
#include <iterator>
#include "tree.h"

class SignalInfo {
public:
    SignalInfo(int dimensions, int64_t signal_width):
        dimensions_(dimensions), signal_width_(signal_width), signal_size_(CalcSignalSize(dimensions, signal_width)), log_signal_width_(CalcLog(signal_width)) {
    }

    int Dimensions() const {
        return dimensions_;
    }

    int64_t SignalWidth() const {
        return signal_width_;
    }

    int64_t LogSignalWidth() const {
        return log_signal_width_;
    }

    int64_t SignalSize() const {
        return signal_size_;
    }

    bool operator==(const SignalInfo& other) const {
        return dimensions_ == other.dimensions_ && signal_width_ == other.signal_width_;
    }

private:
    static int64_t CalcSignalSize(int dimensions, int64_t signal_width) {
        int64_t res = 1;
        for (int i = 0; i < dimensions; ++i) {
            res *= signal_width;
        }
        return res;
    }

    int dimensions_;
    int64_t signal_width_;
    int64_t signal_size_;
    int64_t log_signal_width_;
};


class Key {
public:
    explicit Key(const SignalInfo& info) : info_(info), indices_(info.Dimensions(), 0) {
    }

    explicit Key(const SignalInfo& info, std::vector<int64_t> key) : info_(info), indices_(std::move(key)) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
    }

    explicit Key(const SignalInfo& info, const std::initializer_list<int64_t>& key) : info_(info), indices_(info.Dimensions(), 0) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
        std::copy(key.begin(), key.end(), indices_.begin());
    }

    explicit Key(const SignalInfo& info, int64_t flatten) : info_(info), indices_(info.Dimensions(), 0) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
        SetFromFlatten(flatten);
    }

    void SetZero() {
        indices_.assign(indices_.size(), 0);
    }

    // leftmost dimension is highest in the tree
    const int64_t& operator[](int index) const {
        return indices_[index];
    }

    int64_t& operator[](int index) {
        return indices_[index];
    }

    // in flattened form, least significant dimension is highest
    int64_t Flatten() const {
        int64_t res = 0;
        for (int i = info_.Dimensions() - 1; i >= 0; --i) {
            res = (res << info_.LogSignalWidth()) | indices_[i];
        }
        return res;
    }

    bool operator==(const Key& other) const {
        assert(info_ == other.info_);
        return indices_ == other.indices_;
    }

    SignalInfo GetSignalInfo() const {
        return info_;
    }

    Key IncreaseAt(int index, int64_t value) const {
        std::vector<int64_t> new_indices(indices_);
        new_indices[index] += value + info_.SignalWidth();
        new_indices[index] %= info_.SignalWidth();
        return Key{info_, std::move(new_indices)};
    }

    void StoreDifference(const Key& a, const Key& b) {
        int64_t mod = info_.SignalWidth() - 1;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = (a.indices_[i] - b.indices_[i] + info_.SignalWidth()) & mod;
        }
    }

    Key operator-() const {
        std::vector<int64_t> new_indices(info_.Dimensions());
        for (int i = 0; i < info_.Dimensions(); ++i) {
            new_indices[i] = (-indices_[i] + info_.SignalWidth()) % info_.SignalWidth();
        }
        return Key{info_, std::move(new_indices)};
    }

    int64_t operator*(const Key& key) const {
        int64_t result = 0;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            result += key[i] * indices_[i];
        }
        return result;
    }

    void SetFromFlatten(int64_t flat) {
        int64_t mod = info_.SignalWidth() - 1;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = flat & mod;
            flat >>= info_.LogSignalWidth();
        }
    }

private:
    SignalInfo info_;
    std::vector<int64_t> indices_;
};

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