#pragma once

#include "filter.h"
#include "signal.h"
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "random"
#include "optional"
#include "fftw3.h"
#include "projection_recovery.h"
#include "iostream"
#include <typeinfo>
#include <random>


class IndexGenerator {
public:
    IndexGenerator(const SignalInfo& info, int64_t sparsity, double zero_test_koef, int64_t seed):
        info_(info), rand_gen_(seed), index_gen_(0, info.SignalSize()) {
            
            int64_t table_size = std::max<int64_t>(llround(zero_test_koef * sparsity * log2(info.SignalSize())), 1);
            
            Key time(info);
            for (int64_t iter = 0; iter < table_size; ++iter) {
                Next(time);
                indices_.push_back(time);
                //random_signs_.push_back((index_gen_(rand_gen_)%2 - 0.5)*2.);
            }

    }

    void Next(Key& value) {
        value.SetFromFlatten(index_gen_(rand_gen_));
    }
    
    std::vector<Key> indices_;
    //std::vector<double> random_signs_;

private:
    SignalInfo info_;
    std::mt19937_64 rand_gen_;
    std::uniform_int_distribution<int32_t> index_gen_;
};

struct TransformSettings {
    bool use_preemptive_tests{true};
    double zero_test_koef{1};
    bool use_comb{true};
    bool assume_random_phase{false};
    int random_phase_sparsity_koef{1};
    bool use_projection_recovery{true};
};

class AlignedVector {
public:
    AlignedVector(int n) {
        data_ = nullptr;
        realloc_data(n);
    }

    void expand(int n) {
        if (n > size_) {
            realloc_data(n);
        }
    }

    ~AlignedVector() {
        if (data_ != nullptr) {
            std::free(data_);
        }
    }

    double& operator[](int i) {
        return data_[i];
    }

private:
    int size_;
    double* data_;

    void realloc_data(int n) {
        if (data_ != nullptr) {
            std::free(data_);
            data_ = nullptr;
        }
        if (n % 2 == 1) {
            n += 1;
        }
        data_ = (double *) std::malloc(n * sizeof(double));
        for (int i = 0; i < n; ++i) {
            data_[i] = 0;
        }
        size_ = n;
    }
};

bool ZeroTest(const Signal& x, const FrequencyMap& recovered_freq, const SplittingTree& tree,
              const SplittingTree::NodePtr& cone_node, const SignalInfo& info, int64_t sparsity, IndexGenerator& delta,
              const TransformSettings& settings) {

    auto filter = Filter(tree, cone_node, info);
        
    int64_t max_iters = std::max<int64_t>(llround(settings.zero_test_koef * sparsity * log2(info.SignalSize())), 1);
    std::vector<std::pair<Key, complex_t>> freq_precalc;
    freq_precalc.reserve(recovered_freq.size());
    for (const auto& freq: recovered_freq) {
        complex_t filter_value = filter.FilterFrequency(freq.first);
        if (NonZero(filter_value)) {
            freq_precalc.emplace_back(freq.first, freq.second * filter_value);
        }
    }
        
    Key diff(info);
    Key time(info);
    

    //static AlignedVector phi(freq_precalc.size());
    //phi.expand(freq_precalc.size());
    
    //static AlignedVector x_kernel(freq_precalc.size());
    //x_kernel.expand(freq_precalc.size());
    //static AlignedVector y_kernel(freq_precalc.size());
    //y_kernel.expand(freq_precalc.size());
    static std::vector<complex_t> recovered_sum;
    recovered_sum.clear();
    recovered_sum.assign(freq_precalc.size(), 0);
    static std::vector<complex_t> filtered_sum;
    filtered_sum.clear();
    filtered_sum.assign(filter.FilterTime().size(), 0);
    double phi_koef = 2 * PI / info.SignalWidth();
    
    complex_t total_sum = 0;
    
    //std::cout << "signal size: " << info.SignalSize() << ", sparsity: " << sparsity << ", max iter: " << max_iters << '\n';
    
    for (int64_t iter = 0; iter < max_iters; ++iter) {
        //delta.Next(time);
        time = delta.indices_[iter];
        //double rand_sign = 1.; //delta.random_signs_[iter];
        
        if (info.IsSmallSignalWidth()) {
            for (int j = 0; j < static_cast<int>(freq_precalc.size()); ++j) {
                recovered_sum[j] += complex_t(GetTableCos(freq_precalc[j].first * time, info.SignalWidth()), GetTableSin(freq_precalc[j].first * time, info.SignalWidth()));
            }
        } else {
            for (int j = 0; j < static_cast<int>(freq_precalc.size()); ++j) {
                recovered_sum[j] += complex_t(cos((freq_precalc[j].first * time) * phi_koef), sin((freq_precalc[j].first * time) * phi_koef));
            }
        }
                
        for (int j =0; j < static_cast<int>(filter.FilterTime().size()); ++j) {
            diff.StoreDifference(time, filter.FilterTime()[j].first);
            filtered_sum[j] += x.ValueAtTime(diff);
        }
        
        if (iter < 1){
            total_sum = 0;

            for (int j = 0; j < static_cast<int>(freq_precalc.size()); ++j) {
                total_sum += recovered_sum[j] * freq_precalc[j].second;
            }
            total_sum /= static_cast<double>(info.SignalSize());
            for (int j =0; j < static_cast<int>(filter.FilterTime().size()); ++j) {
                total_sum -= filtered_sum[j] * filter.FilterTime()[j].second;
            }
            if (NonZero(total_sum)) {
                return true;
            }
        }
    }
    
    total_sum = 0;
    for (int j = 0; j < static_cast<int>(freq_precalc.size()); ++j) {
        total_sum += recovered_sum[j] * freq_precalc[j].second;
    }
    total_sum /= static_cast<double>(info.SignalSize());
    for (int j =0; j < static_cast<int>(filter.FilterTime().size()); ++j) {
        total_sum -= filtered_sum[j] * filter.FilterTime()[j].second;
    }
    if (NonZero(total_sum / sqrt(static_cast<double>(max_iters)))) {
        return true;
    }
    //std::cout << "zerotest filtered value: " << total_sum / sqrt(static_cast<double>(max_iters)) << '\n';
    return false;
}

std::optional<FrequencyMap> SparseFFT(const Signal& x, const SignalInfo& info, int64_t expected_sparsity, SplittingTree* parent_tree,
                       SplittingTree::NodePtr& parent_node, const FrequencyMap& known_freq, IndexGenerator& delta, const TransformSettings& settings) {
    FrequencyMap recovered_freq;
    FrequencyMap total_freq = known_freq;
    SplittingTree tree(parent_tree, parent_node);
    
    while (!tree.IsEmpty() && (settings.assume_random_phase || (tree.LeavesCount() + static_cast<int>(recovered_freq.size())) <= expected_sparsity)) {
        NodePtr node = tree.GetLightestNode();
        if (node->level == info.Dimensions() * info.LogSignalWidth()) {
            //std::cout << "now we are peeling a freq" << '\n';
            auto filter = Filter(tree, node, info);
            complex_t recovered = 0;
            for (const auto& freq : total_freq) {
                recovered += freq.second * filter.FilterFrequency(freq.first);
            }
            complex_t filtered = 0;
            for (const auto& value : filter.FilterTime()) {
                filtered += value.second * x.ValueAtTime(-value.first);
            }
            auto value = static_cast<double>(info.SignalSize()) * filtered - recovered;
            auto key = Key(info, node->label);
            recovered_freq[key] += value;
            total_freq[key] += value;
            tree.RemoveNode(node);
        } else {
            node->AddChildren();
            int64_t zerotest_budget = expected_sparsity - (1-settings.assume_random_phase)*(tree.LeavesCount()-2 + static_cast<int>(recovered_freq.size()));
            
            //std::cout << "sparsity: " << expected_sparsity << ", zerotest budget: " << zerotest_budget << '\n';
            
            if (!ZeroTest(x, total_freq, tree, node->left, info, zerotest_budget, delta, settings)) {
                tree.RemoveNode(node->left);
            }
            if (!ZeroTest(x, total_freq, tree, node->right, info, zerotest_budget, delta, settings)) {
                tree.RemoveNode(node->right);
            }
        }
    }
    if (!settings.assume_random_phase) {
        if ((tree.LeavesCount() + static_cast<int>(recovered_freq.size())) > expected_sparsity) {
            return std::nullopt;
        }
    }
    return recovered_freq;
}

class Restorer {
public:
    Restorer(const SplittingTree* parent_tree, SplittingTree::NodePtr& parent_node, const FrequencyMap& known_freq)
        : recovered_freq_(), total_freq_(known_freq), tree_(parent_tree, parent_node) {}

    std::optional<FrequencyMap> TryRestore(const Signal& x, const SignalInfo& info, const std::vector<int>& sparsities,
                                           int rank, IndexGenerator& delta, const TransformSettings& settings) {
        
        
        int64_t expected_sparsity = sparsities[rank - 1];
        next_sparsity_ = sparsities[rank - 2];
        auto root = tree_.GetRoot();
        SparsityTest(x, info, expected_sparsity, sparsities, root, rank, delta, settings);
        while (!tree_.IsEmpty() && (tree_.LeavesCount() * next_sparsity_ + static_cast<int>(recovered_freq_.size())) <= expected_sparsity) {
            NodePtr v = tree_.GetLightestNode();
            // Check if we reached the bottom
            if (v->level == info.Dimensions() * info.LogSignalWidth()) {
                SparsityTest(x, info, expected_sparsity, sparsities, v, rank, delta, settings);
                // Sometimes errors are too big, and the node is not removed after SparsityTest
                if (!v->Removed()) {
                    tree_.RemoveNode(v);
                }
            } else {
                v->AddChildren();
                std::vector<NodePtr> children({v->left, v->right});
                bool skip_restore = false;
                if (settings.use_preemptive_tests) {
                    int64_t zerotest_budget = expected_sparsity - (1-settings.assume_random_phase)*((tree_.LeavesCount() -2) * next_sparsity_ + static_cast<int>(recovered_freq_.size()));
                    for (auto& node : children) {
                        if (!ZeroTest(x, total_freq_, tree_, node, info, zerotest_budget, delta, settings)) {
                            skip_restore = true;
                            tree_.RemoveNode(node);
                        }
                    }
                }
                if (!skip_restore) {
                    for (auto& node : children) {
                        SparsityTest(x, info, expected_sparsity, sparsities, node, rank, delta, settings);
                    }
                }
            }
        }
        if ((tree_.LeavesCount() * next_sparsity_ + static_cast<int>(recovered_freq_.size())) > expected_sparsity) {
            return std::nullopt;
        }
        return std::move(recovered_freq_);
    }

private:
    void SparsityTest(const Signal& x, const SignalInfo& info, int64_t expected_sparsity, const std::vector<int>& sparsities,
                                           SplittingTree::NodePtr& node, int rank, IndexGenerator& delta, const TransformSettings& settings) {
        std::optional<FrequencyMap> probable_freq(std::nullopt);
        if (rank == 2) {
            probable_freq = SparseFFT(x, info, next_sparsity_, &tree_, node, total_freq_, delta, settings);
        } else {
            Restorer restorer(&tree_, node, total_freq_);
            probable_freq = restorer.TryRestore(x, info, sparsities, rank - 1,
                                                delta, settings);
        }
        if (probable_freq) {
            auto new_total_freq = MapUnion(probable_freq.value(), total_freq_);
            int64_t zerotest_budget = expected_sparsity - (1-settings.assume_random_phase)*(std::max<int>(0, tree_.LeavesCount() -2) * next_sparsity_ + static_cast<int>(recovered_freq_.size()));
            if (!ZeroTest(x, new_total_freq, tree_, node, info, zerotest_budget, delta, settings)) {
                tree_.RemoveNode(node);
                recovered_freq_ = MapUnion(recovered_freq_, probable_freq.value());
                total_freq_.swap(new_total_freq);
            }
        }
    }

    FrequencyMap recovered_freq_;
    FrequencyMap total_freq_;
    int64_t next_sparsity_;
    SplittingTree tree_;
};

std::vector<int> PrepareCombSizes(int n, int d, int k) {
    const double Comb_cst = 10;
    std::vector<int> W_Combs(d);
    int total_log = static_cast<int>(floor(std::min(d * log2(n), log2(Comb_cst * k))));
    int base = total_log / d;
    for (int i = 0; i < total_log % d; ++i) {
        W_Combs[i] = 1 << (base + 1);
    }
    for (int i = total_log % d; i < d; ++i) {
        W_Combs[i] = 1 << base;
    }
    return W_Combs;
}

void ComputeCombBucketedSignal(int n, int d, int lvl, const Key& a, int* W_Combs, int* u_index, Key& in_index, const Signal& in, complex_t* u) {
    if (lvl < d) {
        int step = n / W_Combs[lvl];
        for (int h = 0; h < W_Combs[lvl]; ++h) {
            u_index[lvl] = h;
            in_index[lvl] = (h * step + a[lvl]) & (n - 1);
            ComputeCombBucketedSignal(n, d, lvl + 1, a, W_Combs, u_index, in_index, in, u);
        }
    } else {
        // dimension order for ffrw is row major (different from Key order)
        int flatten = 0;
        for (int i = 0; i < d; ++i) {
            flatten = flatten * W_Combs[i] + u_index[i];
        }
        u[flatten] += in.ValueAtTime(in_index);
    }
}

void HashToBinsComb(const SignalInfo& info, const Signal& in, const Key& a, int W_total, int* W_Combs, complex_t* u, const fftw_plan& p) {
    std::vector<int> u_index(info.Dimensions(), 0);
    Key in_index(info);
    for (int i = 0; i < W_total; ++i) {
        u[i] = 0;
    }
    ComputeCombBucketedSignal(info.SignalWidth(), info.Dimensions(), 0, a, W_Combs, u_index.data(), in_index, in, u);
    fftw_execute(p);
    for (int i = 0; i < W_total; ++i) {
        u[i] *= info.SignalSize() / W_total;
    }
}

int RestoreFrequencies(const SignalInfo& info, int Btotal, const Key& ai, const std::vector<complex_t*>& u, FrequencyMap& out) {
    Key i(info);
    int cnt = 0;
    for (int j = 0; j < Btotal; ++j) {
        if (NonZero(u[0][j])) {
            ++cnt;
            i.SetZero();
            for (int h = 0; h < info.Dimensions(); ++h) {
                complex_t alpha = u[0][j] / u[h + 1][j];
                i[h] = (((int64_t) ai[h]) * lround(std::arg(alpha) * info.SignalWidth() / (2 * M_PI))) & (info.SignalWidth() - 1);
                if (i[h] < 0) {
                    i[h] += info.SignalWidth();
                }
            }
            out[i] += u[0][j];
        }
    }
    return cnt;
}

int CombFiltration(const Signal& x, const SignalInfo& info, int64_t sparsity, FrequencyMap& out) {
    int n = info.SignalWidth();
    int d = info.Dimensions();
    int WCombTotal = 1;
    auto W_Combs = PrepareCombSizes(n, d, sparsity);
    for (int i = 0; i < d; ++i) {
        WCombTotal *= W_Combs[i];
    }
    std::vector<complex_t*> u(d + 1);
    for (int i = 0; i < d + 1; ++i) {
        u[i] = (complex_t*) fftw_malloc(sizeof(fftw_complex) * WCombTotal);
    }
    Key c(info);
    for (int i = 0; i < d + 1; ++i) {
        if (i > 0) {
            c[i - 1] = n - 1;
        }
        fftw_plan plan = fftw_plan_dft(d, W_Combs.data(), reinterpret_cast<fftw_complex*>(u[i]), reinterpret_cast<fftw_complex*>(u[i]), FFTW_FORWARD, FFTW_ESTIMATE);
        HashToBinsComb(info, x, c, WCombTotal, W_Combs.data(), u[i], plan);
        fftw_destroy_plan(plan);
        if (i > 0) {
            c[i - 1] = 0;
        }
    }
    Key ai(info);
    for (int i = 0; i < d; ++i) {
        ai[i] = 1;
    }
    auto cnt = RestoreFrequencies(info, WCombTotal, ai, u, out);
    for (int i = 0; i < d + 1; ++i) {
        fftw_free(u[i]);
    }
    return cnt;
}

FrequencyMap RecursiveSparseFFT(const Signal& x, const SignalInfo& info, int64_t sparsity, int rank, int64_t seed = 61, TransformSettings settings = {}) {
    assert(info.SignalSize() > 1);
    if (sparsity == 0) {
        return {};
    }
    if (info.IsSmallSignalWidth()) {
        PrepareCosSinTables(info.SignalWidth());
    }
    
    //std::cout << "rank: " << rank << ", signal size: " << info.SignalSize() << '\n';
    
    FrequencyMap prefiltered;
    if (settings.use_comb) {
        CombFiltration(x, info, sparsity, prefiltered);
        sparsity += prefiltered.size();
    }
    
    if (settings.use_projection_recovery && (info.SignalSize() / info.SignalWidth() < settings.zero_test_koef*pow(sparsity,3)) ) {
        // PFT is assumed to be always correct, therefore we can just overwrite frequencies and not change the sparsity
        auto res = ProjectionFT(x, info, 2);
//        printf("width: %i, recovered: %i\n", (int) info.SignalWidth(), (int) res.size());
        for (const auto& w: res) {
            prefiltered[w.first] = w.second;
        }
        sparsity = std::max<int64_t>(1, sparsity - res.size());
        //std::cout << "log n: " << info.LogSignalWidth() << ", rank: " << rank << ", sparsity after projection: " << sparsity << '\n';
    }
    
    NodePtr parent(nullptr);
    if (settings.assume_random_phase) {
        rank = 1;
        sparsity = std::min<int64_t>(sparsity, settings.random_phase_sparsity_koef);
    }
    
    IndexGenerator delta(info, sparsity, settings.zero_test_koef, seed);
    
    double step = pow(sparsity / 0.5, 1. / rank);
    std::vector<int> sparsities(rank);
    sparsities[rank - 1] = sparsity;
    double curspars = sparsity / step;
    for (int i = rank - 2; i >= 0; --i) {
        sparsities[i] = std::max<int>(1, int(curspars));
        curspars /= step;
    }
    std::optional<FrequencyMap> res;
        
    if (rank == 1) {
        res = SparseFFT(x, info, sparsity, nullptr, parent, prefiltered, delta, settings);
    } else {
        Restorer restorer(nullptr, parent, prefiltered);
        res = restorer.TryRestore(x, info, sparsities, rank, delta, settings);
    }
    if (!res) {
        return prefiltered;
    }

    return MapUnion(prefiltered, res.value());
}
