#pragma once

#include "filter.h"
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "random"
#include <array>

// можно вообще не хранить индексы и хранить только хеш

//struct Key {
//    Key(int16_t dimension) : dimension_(dimension) {
//        indices_ = malloc(dimension_ * sizeof(int64_t));
//    }
//
//    Key(const std::vector<int64_t>& key) : dimension_(key.size()){
//        indices_ = malloc(dimension_ * sizeof(int64_t));
//        memcpy(indices_, key.data(), dimension_ * sizeof(int64_t));
//    }
//
//    Key(const std::initializer_list<int64_t>& key) = 0;
//
//    ~Key() {
//        free(indices_);
//    }
//
//    int64_t& operator[] (int16_t index) {
//        return indices_[index];
//    }
//
//    int16_t size() const {
//        return dimension_;
//    }
//
//
//    int16_t dimension_;
//    int64_t* indices_;
//}; // do not delete

template <class T>
inline void hash_combine(size_t& s, const T& value) {
    std::hash<T> hash;
    s ^= hash(value) + 0x9e3779b9 + (s << 6) + (s >> 2);
}

//namespace std {
//    template<>
//    struct hash<Key>
//    {
//        size_t operator()(const Key& key) const {
//            size_t result = 0;
//            for (int16_t i = 0; i < key.size(); ++i) {
//                hash_combine(result, key[i]);
//            }
//            return result;
//        }
//    };
//} // do not delete


// I have tried to use initializer_list, but there is no random access operators in it
// TODO rewrite with initializer_list if needed

namespace std {
    template<T>
    struct hash<std::initializer_list<T>>
    {
        size_t operator()(const std::initializer_list<T>& key) const {
            size_t result = 0;
            for (auto i : key) {
                hash_combine(result, i);
            }
            return result;
        }
    };
}

namespace std {
    template<class T>
    struct equal_to<std::initializer_list<T>>
    {
        bool operator()(const std::initializer_list<T>& l, const std::initializer_list<T>& r) const {
            auto hash = std::hash<std::initializer_list<T>>();
            return hash(l) == hash(r);
        }
    };
}


using Key = std::initializer_list<int64_t>;
using FrequencyMap = std::unordered_map<Key, complex_t>;
using NodePtr = SplittingTree::NodePtr;
using Node = SplittingTree::Node;

class Signal {
public:
    virtual ~Signal() = default;

    virtual complex_t ValueAtTime(const Key& key) const = 0;
};

//TODO rewrite with some suitable tensor library


// Зачем нужна эта структура??
class DataSignal: public Signal {
public:
    DataSignal(int64_t signal_size, int16_t dimension, FrequencyMap&& v):
            signal_size_(signal_size),
            dimension_(dimension),
            values_(std::move(v)) {
    }

    complex_t ValueAtTime(const Key& key) const override {
        return values_[key];
    }

private:
    const int64_t signal_size_;
    const int16_t dimension_;
    FrequencyMap values_;
};

class IndexGenerator {
public:
    IndexGenerator(int64_t signal_size, int16_t dimension, int64_t seed): rand_gen_(seed), index_gen_(0, signal_size), dimension_(dimension) {
    }

    //how not to copy
    Key Next() {
        return index_gen_(rand_gen_);
    }

private:
    std::mt19937_64 rand_gen_;
    std::uniform_int_distribution<int64_t> index_gen_;
    int16_t dimension_;
};

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree::NodePtr& cone_node, int64_t signal_size, int64_t sparsity, IndexGenerator& delta) {
    auto filter = Filter(cone_node, signal_size);
    int64_t max_iters = std::max<int64_t>(llround(2 * sparsity * log2(sparsity) * log2(sparsity) * log2(signal_size)), 2);
    for (int64_t i = 0; i < max_iters; ++i) {
        auto time = delta.Next();
        complex_t recovered_at_time = 0;
        for (auto freq: recovered_freq) {
            recovered_at_time += CalcKernel(freq.first * time, signal_size) * freq.second * filter.FilterFrequency(freq.first) / complex_t(signal_size, 0);
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