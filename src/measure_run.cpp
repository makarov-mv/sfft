#include "iostream"
#include "chrono"
#include "disfft.h"
#include "utility_test_fftw.h"

template <class Generator>
std::vector<complex_t> GenRandomSupport(const SignalInfo& info, int64_t sparsity, Generator& gen) {
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < sparsity; ++i) {
        out[dist(gen)] = 1;
    }
    return out;
}

std::vector<complex_t> GenDiracComb(const SignalInfo& info, int64_t sparsity) {
    assert((sparsity & (sparsity - 1)) == 0);
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = CalcKernel(i * sparsity / info.SignalWidth(), sparsity);
    }
    return out;
}

template <class Generator>
std::vector<complex_t> GenCombined(const SignalInfo& info, int64_t sparsity, Generator& gen) {
    assert((sparsity & (sparsity - 1)) == 0);
    std::vector<complex_t> out = GenRandomSupport(info, sparsity / 2, gen);
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / (sparsity / 2)) {
        out[i + 1] += CalcKernel(i * sparsity / 2 / info.SignalWidth(), sparsity / 2);
    }
    return out;
}

std::vector<complex_t> Gen4DComb(const SignalInfo& info, int64_t sparsity) {
    assert(info.Dimensions() == 4);
    int64_t root = int(sqrt(sparsity));
    assert(root * root == sparsity);
    assert(info.SignalWidth() % root == 0);
    int64_t step = info.SignalWidth() * info.SignalWidth() * info.SignalWidth();
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalWidth(); i += info.SignalWidth() / root) {
        for (int j = 0; j < info.SignalWidth(); j += info.SignalWidth() / root) {
            out[i * step + j + 1] = CalcKernel(j * root / info.SignalWidth(), root);
        }
    }
    return out;
}

void PrintDur(const std::chrono::nanoseconds& dur) {
    std::cout << std::chrono::duration<double, std::milli>(dur).count();
}

void PrintArr(const std::string& name, const std::vector<std::chrono::nanoseconds>& durs) {
    std::cout << name << " = [\n";
    for (auto p: durs) {
        std::cout << p.count() << ",\n";
    }
    std::cout << "]\n";
}

template <class Func>
auto RunBenchmark(const std::string& name, const std::vector<std::vector<complex_t>>& signals, Func function, std::vector<std::chrono::nanoseconds>& res) {
    using clock = std::chrono::system_clock;
    std::cout << name << ": ";
    auto start = clock::now();
    for (int i = 0; i < static_cast<int>(signals.size()); ++i) {
        (void) function(signals[i], i);
    }
    auto dur = clock::now() - start;
    dur /= signals.size();
    PrintDur(dur);
    res.push_back(dur);
    std::cout << ", ";
}

class Benchmark {
public:
    virtual ~Benchmark() = default;
    virtual int Sparsity() = 0;
    virtual std::vector<TransformSettings> PrepareSettings() = 0;
    virtual std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) = 0;
    virtual int Start() = 0;
    virtual int End() = 0;
    virtual SignalInfo GetInfo(int p) = 0;
};

class BenchmarkRandomSupport : public Benchmark {
public:
    BenchmarkRandomSupport(bool use_comb) : use_comb_(use_comb) {}

    int Sparsity() override {
        return 27;
    }

    std::vector<TransformSettings> PrepareSettings() override {
        std::vector<TransformSettings> res(5);
        TransformSettings settings;
        settings.use_comb = use_comb_;
        settings.zero_test_koef = 0.03;
        for (int i = 0; i < 4; ++i) {
            res[i] = settings;
        }
        settings.zero_test_koef = 1;
        settings.random_phase_sparsity_koef = 1;
        settings.assume_random_phase = true;
        res[4] = settings;
        return res;
    }

    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
        return GenRandomSupport(info, sparsity, gen);
    }

    int Start() override {
        return 3;
    }

    int End() override {
        return 7;
    }

    SignalInfo GetInfo(int p) override {
        return {3, 1 << p};
    }

private:
    bool use_comb_;
};

class BenchmarkCombinedSupport : public Benchmark {
public:
    BenchmarkCombinedSupport(bool use_comb) : use_comb_(use_comb) {}

    int Sparsity() override {
        return 32;
    }

    std::vector<TransformSettings> PrepareSettings() override {
        std::vector<TransformSettings> res(5);
        TransformSettings settings;
        settings.use_comb = use_comb_;
        settings.zero_test_koef = 0.5;
        for (int i = 0; i < 4; ++i) {
            res[i] = settings;
        }
        settings.zero_test_koef = 1;
        settings.random_phase_sparsity_koef = 4;
        settings.assume_random_phase = true;
        res[4] = settings;
        return res;
    }

    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
        return GenCombined(info, sparsity, gen);
    }

    int Start() override {
        return 3;
    }

    int End() override {
        return 7;
    }

    SignalInfo GetInfo(int p) override {
        return {3, 1 << p};
    }

private:
    bool use_comb_;
};

class Benchmark4DCombSupportWithComb : public Benchmark {
public:
    Benchmark4DCombSupportWithComb(bool use_comb) : use_comb_(use_comb) {}

    int Sparsity() override {
        return comb_sparsity + rand_supp;
    }

    std::vector<TransformSettings> PrepareSettings() override {
        std::vector<TransformSettings> res(5);
        TransformSettings settings;
        settings.use_comb = use_comb_;
        settings.zero_test_koef = 0.5;
        for (int i = 0; i < 4; ++i) {
            res[i] = settings;
        }
        settings.zero_test_koef = 1;
        settings.random_phase_sparsity_koef = 10;
        settings.assume_random_phase = true;
        res[4] = settings;
        return res;
    }

    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
        assert(sparsity == comb_sparsity + rand_supp);
        auto out = Gen4DComb(info, comb_sparsity);
        std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);
        for (int i = 0; i < rand_supp; ++i) {
            out[dist(gen)] = 1;
        }
        return out;
    }

    int Start() override {
        return 3;
    }

    int End() override {
        return 5;
    }

    SignalInfo GetInfo(int p) override {
        return {4, 1 << p};
    }

private:
    bool use_comb_;
    const int comb_sparsity = 64;
    const int rand_supp = 64;
};

int main() {
    std::mt19937_64 gen(123423);
    std::vector<int64_t> npow;
    std::vector<std::chrono::nanoseconds> dur_fftw;
    std::vector<std::chrono::nanoseconds> dur1;
    std::vector<std::chrono::nanoseconds> dur2;
    std::vector<std::chrono::nanoseconds> dur3;
    std::vector<std::chrono::nanoseconds> dur4;
    std::vector<std::chrono::nanoseconds> dur_phase;
    auto bench = BenchmarkRandomSupport(false);

    for (int64_t p = bench.Start(); p <= bench.End(); ++p) {
        SignalInfo info = bench.GetInfo(p);
        const int64_t sparsity = bench.Sparsity();
        const int samples = 100;

        auto settings = bench.PrepareSettings();

        std::vector<std::vector<complex_t>> signals;
        auto runner = FFTWRunner(info, FFTW_BACKWARD);
        for (int i = 0; i < samples; ++i) {
            auto out = bench.GenSignal(info, sparsity, gen);
            signals.emplace_back(runner.Run(out));
        }
        auto reverse = FFTWRunner(info, FFTW_FORWARD);
        npow.push_back(p);
        std::cout << "p = " << p << ", ";
        using Input = const std::vector<complex_t>&;
        RunBenchmark("fftw", signals, [&](Input in, int){return reverse.Run(in);}, dur_fftw);
        RunBenchmark("rank 1", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 1, i, settings.at(0));}, dur1);
        RunBenchmark("rank 2", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 2, i, settings.at(1));}, dur2);
        RunBenchmark("rank 3", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 3, i, settings.at(2));}, dur3);
        RunBenchmark("rank 4", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 4, i, settings.at(3));}, dur4);

        RunBenchmark("rand phase", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 1, i, settings.at(4));}, dur_phase);

        std::cout << std::endl;
    }
    std::cout << "p = [\n";
    for (auto p: npow) {
        std::cout << p << ",\n";
    }
    std::cout << "]\n";
    PrintArr("fftw", dur_fftw);
    PrintArr("rank1", dur1);
    PrintArr("rank2", dur2);
    PrintArr("rank3", dur3);
    PrintArr("rank4", dur4);
    PrintArr("rand_phase", dur_phase);
    return 0;
}

