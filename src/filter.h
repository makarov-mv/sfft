#pragma once

#include <vector>
#include <unordered_map>
#include "tree.h"

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


class Filter {
public:

    Filter(const NodePtr& node, const SignalInfo& info) {
        std::vector<bool> path_mask = node->GetRootPathMask(); // optimize
        phase_.resize(CalcLog(size), 0);
        auto label = node->label;

        for (int i = 0; i + 1 < static_cast<int>(path_mask.size()); ++i) {
            if (path_mask[i]) {
                phase_[i] = CalcKernel(-label, 1 << (i + 1));
            }
        }

        filter_[0] = 1;
        for (int i = 0; i < static_cast<int>(phase_.size()); ++i) {
            if (NonZero(phase_[i])) {
                std::unordered_map<int64_t, complex_t> updated_filter;
                int64_t shift = size >> (i + 1);
                for (auto it : filter_) {
                    int64_t index = it.first;
                    updated_filter[index] = (FilterValueAtTime(index) + phase_[i] * FilterValueAtTime((index + shift) % size)) / 2.;
                    updated_filter[(index + size - shift) % size] = (FilterValueAtTime((index + size - shift) % size) + phase_[i] * FilterValueAtTime(index)) / 2.;
                }
                filter_.swap(updated_filter);
                updated_filter.clear();
            }
        }
    }

    complex_t FilterValueAtTime(int64_t time) const {
        auto it = filter_.find(time);
        if (it != filter_.end()) {
            return it->second;
        }
        return 0.;
    }

    const std::unordered_map<Key, complex_t>& FilterTime() const {
        return filter_;
    }

    complex_t FilterFrequency(const Key& key) const {
        complex_t freq =  1;
        for (size_t i = 0; i < phase_.size(); ++i) {
            if (NonZero(phase_[i])) {
                freq *= (1. + phase_[i] * CalcKernel(psi, 1 << (i + 1))) / 2.;
            }
        }
        return freq;
    }

private:

    std::vector<complex_t> phase_;
    std::unordered_map<Key, complex_t> filter_;
};
