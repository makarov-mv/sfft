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

    const Key operator-(const Key& key) const {
        std::vector<int64_t> new_indices(info_.Dimensions());
        for (size_t i = 0; i < info_.Dimensions(); ++i) {
            new_indices[i] = (indices_[i] - key[i] + info_.SignalWidth()) % info_.SignalWidth();
        }
        return {info, std::move(new_indices)};
    }

    const Ket operator-() const {
        std::vector<int64_t> new_indices(info_.Dimensions());
        for (size_t i = 0; i < info_.Dimensions(); ++i) {
            new_indices[i] = (-indices_[i] + info_.SignalWidth()) % info_.SignalWidth();
        }
        return {info, std::move(new_indices)};
    }

    double operator*(const Key& key) const {
        int64_t result = 0;
        for (size_t i = 0; i < info_.Dimensions(); ++i) {
            result += key[i] * indices_[i];
        }
        return result;
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
//filter.FilterTime() consists of all no zero freq return map(key, complex_t)

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree::NodePtr& cone_node, const SignalInfo& info, int64_t sparsity, IndexGenerator& delta) {
    auto filter = Filter(cone_node, info.SignalWidth(), info.Dimensions());
    int64_t max_iters = std::max<int64_t>(llround(2 * sparsity * log2(sparsity) * log2(sparsity) * log2(info.SignalSize())), 2); // check
    for (int64_t i = 0; i < max_iters; ++i) {
        auto time = delta.Next();
        complex_t recovered_at_time = 0;
        for (const auto& freq: recovered_freq) {
            recovered_at_time += CalcKernel(freq.first * time, info.SignalWidth()) * freq.second * filter.FilterFrequency(freq.first);
        }
        recovered_at_time /= static_cast<double>(info.SignalSize());
        complex_t filtered_at_time = 0;
        for (const auto& value: filter.FilterTime()) {
            filtered_at_time += value.second * x.ValueAtTime(time - value.first);
        }
        if (NonZero(filtered_at_time - recovered_at_time)) {
            return true;
        }
    }
    return false;
}

FrequencyMap SparseFFT(const Signal& x, const SignalInfo& info, int64_t sparsity) {
    assert(info.SignalSize() > 1);
    SplittingTree tree{};
    FrequencyMap recovered_freq;
    IndexGenerator delta{info, 61};

    while (!tree.IsEmpty()) {
        NodePtr node = tree.GetLightestNode();
        if (node->level == info.Dimensions() * CalcLog(info.SignalSize())) {
            auto filter = Filter(node, info.SignalWidth(), info.Dimensions());
            complex_t recovered = 0;
            for (const auto& freq : recovered_freq) {
                recovered += freq.second * filter.FilterFrequency(freq.first);
            }
            complex_t filtered = 0;
            for (const auto& value : filter.FilterTime()) {
                filtered += value.second * x.ValueAtTime(-value.first);
            }
            recovered_freq[Key(info, node->label)] += static_cast<double>(info.SignalSize()) * filtered - recovered;
            tree.RemoveNode(node);
        } else {
            node->AddChildren();
            if (!ZeroTest(x, recovered_freq, node->left, info, sparsity, delta)) {
                tree.RemoveNode(node->left);
            }
            if (!ZeroTest(x, recovered_freq, node->right, info, sparsity, delta)) {
                tree.RemoveNode(node->right);
            }
        }
    }
    return recovered_freq;
}