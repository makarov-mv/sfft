#pragma once

#include "filter.h"
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "random"
#include <array>

class SignalInfo {
public:
    SignalInfo(int dimensions, int64_t signal_width): dimensions_(dimensions), signal_width_(signal_width) {
    }

    int Dimensions() const {
        return dimensions_;
    }

    int64_t SignalWidth() const {
        return signal_width_;
    }

    int64_t SignalSize() const {
        int64_t res = 1;
        int exp = dimensions_;
        int64_t base = signal_width_;
        for (;;) {
            if (exp & 1) {
                res *= base;
            }
            exp >>= 1;
            if (!exp) {
                return res;
            }
            base *= base;
        }
    }

    bool operator==(const SignalInfo& other) const {
        return dimensions_ == other.dimensions_ && signal_width_ == other.signal_width_;
    }

private:
    int dimensions_;
    int64_t signal_width_;
};


class Key {
public:
    explicit Key(const SignalInfo& info) : info_(info), indices_(info.Dimensions(), 0) {
    }

    explicit Key(const SignalInfo& info, std::vector<int64_t> key) : info_(info), indices_(std::move(key)) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
    }

    explicit Key(const SignalInfo& info, const std::initializer_list<int64_t>& key) : info(info_), indices_(info.Dimensions(), 0) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
        std::copy(key.begin(), key.end(), indices_.begin());
    }

    explicit Key(const SignalInfo& info, int64_t flatten) : info_(info), indices_(info.Dimensions(), 0) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
        SetFromFlatten(flatten);
    }

    int64_t& operator[](int index) {
        return indices_[index];
    }

    int64_t Flatten() const {
        int64_t res = 0;
        for (auto ind : indices_) {
            res = res * info_.SignalWidth() + ind;
        }
        return res;
    }

    bool operator==(const Key& other) const {
        assert(info_ == other.info_);
        return indices_ == other.indices_;
    }

    SignalInfo SignalInfo() const {
        return info_;
    }

private:
    void SetFromFlatten(int64_t flat) {
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = flat % info_.SignalWidth();
            flat /= info_.SignalWidth();
        }
        std::reverse(indices_.begin(), indices_.end());
    }

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

class Signal {
public:
    virtual ~Signal() = default;

    virtual complex_t ValueAtTime(const Key& key) const = 0;
};

class DataSignal: public Signal {
public:
    DataSignal(const SignalInfo& info, FrequencyMap v):
            info_(info),
            values_(std::move(v)) {
    }

    complex_t ValueAtTime(const Key& key) const override {
        assert(key.SignalInfo() == info_);
        return values_[key];
    }

private:
    const SignalInfo info_;
    FrequencyMap values_;
};

class IndexGenerator {
public:
    IndexGenerator(const SignalInfo& info, int64_t seed):
        info_(info), rand_gen_(seed), index_gen_(0, info.SignalSize()) {
    }

    Key Next() {
        return {info_, index_gen_(rand_gen_)};
    }

private:
    SignalInfo info_;
    std::mt19937_64 rand_gen_;
    std::uniform_int_distribution<int64_t> index_gen_;
};

//filter.FilterFrequency(const Key& key)
//operator- for Key
//value at time by modulo signal_size
//filter.FilterTime() consists of all no zero freq

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree::NodePtr& cone_node, const SignalInfo& info, int64_t sparsity, IndexGenerator& delta) {
    auto filter = Filter(cone_node, info.SignalWidth(), info.Dimensions());
    int64_t max_iters = 100;//std::max<int64_t>(llround(2 * sparsity * log2(sparsity) * log2(sparsity) * log2(signal_size)), 2);
    for (int64_t i = 0; i < max_iters; ++i) {
        auto time = delta.Next();
        complex_t recovered_at_time = 0;
        for (auto freq: recovered_freq) {
            auto dot = freq.first * time;
            recovered_at_time += CalcKernel(dot, signal_size) * freq.second * filter.FilterFrequency(freq.first) / complex_t(signal_size ^ dimension, 0);
        }
        complex_t filtered_at_time = 0;
        for (auto value: filter.FilterTime()) {
            filtered_at_time += value.second * x.ValueAtTime(time - value.first);
        }
        if (NonZero(filtered_at_time - recovered_at_time)) {
            return true;
        }
    }
    return false;
}

FrequencyMap SparseFFT(const Signal& x, int64_t signal_size, int64_t sparsity) {
    assert(signal_size > 1);
    SplittingTree tree{};
    FrequencyMap recovered_freq;
    IndexGenerator delta{signal_size, 61};

    while (!tree.IsEmpty()) {
        NodePtr node = tree.GetLightestNode();
        if (node->level == CalcLog(signal_size)) {
            auto filter = Filter(node, signal_size);
            complex_t recovered = 0;
            for (auto freq : recovered_freq) {
                recovered += freq.second * filter.FilterFrequency(freq.first);
            }
            complex_t filtered = 0;
            for (auto value : filter.FilterTime()) {
                filtered += value.second * x.ValueAtTime((signal_size - value.first) % signal_size);
            }
            recovered_freq[node->label] += signal_size * 1. * filtered - recovered;
            tree.RemoveNode(node);
        } else {
            node->AddChildren();
            if (!ZeroTest(x, recovered_freq, node->left, signal_size, sparsity, delta)) {
                tree.RemoveNode(node->left);
            }
            if (!ZeroTest(x, recovered_freq, node->right, signal_size, sparsity, delta)) {
                tree.RemoveNode(node->right);
            }
        }
    }
    return recovered_freq;
}