#pragma once

#include <vector>
#include <unordered_map>
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
    explicit Key(const SignalInfo& info) : info_(info), flattened_(0) {
    }

    explicit Key(const SignalInfo& info, std::vector<int64_t> key) : info_(info), flattened_(0) {
        assert(info.Dimensions() == static_cast<int>(key.size()));
        for (auto it = key.rbegin(); it != key.rend(); ++it) {
            flattened_ = flattened_ * info.SignalWidth() + *it;
        }
    }

    explicit Key(const SignalInfo& info, const std::initializer_list<int64_t>& key) : Key(info, std::vector<int64_t>(key)) {
    }

    explicit Key(const SignalInfo& info, int64_t flatten) : info_(info), flattened_(flatten) {
    }

    // leftmost dimension is highest in the tree
    int64_t operator[](int64_t index) const {
        return (flattened_ & ((info_.SignalWidth() - 1) << (index * info_.LogSignalWidth()))) >> (index * info_.LogSignalWidth());
    }

    // in flattened form, least significant dimension is highest
    int64_t Flatten() const {
        return flattened_;
    }

    bool operator==(const Key& other) const {
        assert(info_ == other.info_);
        return flattened_ == other.flattened_;
    }

    SignalInfo GetSignalInfo() const {
        return info_;
    }

    Key IncreaseAt(int64_t index, int64_t value) const {
        int64_t new_flattened = flattened_;
        int64_t old_value = this->operator[](index);
        new_flattened -= old_value << (index * info_.LogSignalWidth());
        old_value += value + info_.SignalWidth();
        old_value %= info_.SignalWidth();
        new_flattened += old_value << (index * info_.LogSignalWidth());
        return Key{info_, new_flattened};
    }

    Key operator-(const Key& key) const {
        int64_t new_flattened = 0;
        int64_t a = flattened_;
        int64_t b = key.flattened_;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            new_flattened += ((
                (a & (info_.SignalWidth() - 1))
                - (b & (info_.SignalWidth() - 1))
                + info_.SignalWidth()) & (info_.SignalWidth() - 1)) << (i * info_.LogSignalWidth());
            a >>= info_.LogSignalWidth();
            b >>= info_.LogSignalWidth();
        }
        return Key{info_, new_flattened};
    }

    Key operator-() const {
        Key zero{info_, 0};
        return zero - *this;
    }

    int64_t operator*(const Key& key) const {
        int64_t result = 0;
        int64_t a = flattened_;
        int64_t b = key.flattened_;
        int64_t mask = info_.SignalWidth() - 1;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            result += (a & mask) * (b & mask);
            a >>= info_.LogSignalWidth();
            b >>= info_.LogSignalWidth();
        }
        return result;
    }

private:

    SignalInfo info_;
    int64_t flattened_;
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
    c.insert(a.begin(), a.end());
    c.insert(b.begin(), b.end());
    return c;
}

class Filter {
public:

    Filter(const SplittingTree& tree, const NodePtr& node, const SignalInfo& info) : label_(info, node->label), info_(info), period_size_(CalcLog(info.SignalWidth())) {
        path_ = tree.GetRootPath(node); //
        filter_[Key(info, 0)] = 1; //
        for (int i = 0; i < static_cast<int>(path_.size()); ++i) {
            std::unordered_map<Key, complex_t> updated_filter; //
            int current_period = CalcCurrentPeriod(i); //
            int64_t subtree_level = CalcSubtreeLevel(i);
            int64_t shift = info.SignalWidth() >> subtree_level; //
            auto phase = CalcKernel(-label_[current_period], 1 << subtree_level); //
            phase_.push_back(phase);
            for (auto it : filter_) {
                Key index = it.first; //
                updated_filter[index] = (FilterValueAtTime(index) + phase * FilterValueAtTime(index.IncreaseAt(current_period, shift))) / 2.; //
                updated_filter[index.IncreaseAt(current_period, -shift)] = (FilterValueAtTime(index.IncreaseAt(current_period, -shift)) + phase * FilterValueAtTime(index)) / 2.; //
            }
            filter_.swap(updated_filter);
            updated_filter.clear();
        }
    }

    const std::unordered_map<Key, complex_t>& FilterTime() const {
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
        auto it = filter_.find(time);
        if (it != filter_.end()) {
            return it->second;
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

// TODO: Calc frequencies faster
//    complex_t CalcFrequencyFactor(double g) const {
//        return {(1 + cos(g)) / 2, sin(g) / 2};
//    }

    Key label_;
    SignalInfo info_;
    int period_size_;
    std::vector<int> path_;
    std::vector<complex_t> phase_;
    std::unordered_map<Key, complex_t> filter_;
};