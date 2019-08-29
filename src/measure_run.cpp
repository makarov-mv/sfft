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
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / sparsity) {
        out[i] = 1;
    }
    return out;
}

template <class Generator>
std::vector<complex_t> GenCombined(const SignalInfo& info, int64_t sparsity, Generator& gen) {
    std::vector<complex_t> out = GenRandomSupport(info, sparsity / 2, gen);
    for (int i = 0; i < info.SignalSize(); i += info.SignalSize() / (sparsity / 2)) {
        out[i] += 1;
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

int main() {
    std::mt19937_64 gen(123423);
    std::vector<int64_t> npow;
    std::vector<std::chrono::nanoseconds> dur_fftw;
    std::vector<std::chrono::nanoseconds> dur1;
    std::vector<std::chrono::nanoseconds> dur2;
    std::vector<std::chrono::nanoseconds> dur3;
    std::vector<std::chrono::nanoseconds> dur4;
    TransformSettings settings;
    settings.use_comb = true;
    settings.assume_random_phase = true;
    settings.random_phase_sparsity_koef = 1;
    for (int64_t p = 3; p <= 7; ++p) {
        SignalInfo info{3, 1 << p};
        const int64_t sparsity = 50;
        const int samples = 5;
        std::vector<std::vector<complex_t>> signals;
        auto runner = FFTWRunner(info, FFTW_BACKWARD);
        for (int i = 0; i < samples; ++i) {
            auto out = GenRandomSupport(info, sparsity, gen);
            signals.emplace_back(runner.Run(out));
        }
        auto reverse = FFTWRunner(info, FFTW_FORWARD);
        npow.push_back(p);
        std::cout << "p = " << p << ", ";
        using Input = const std::vector<complex_t>&;
        RunBenchmark("fftw", signals, [&](Input in, int){return reverse.Run(in);}, dur_fftw);
        RunBenchmark("rank 1", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 1, i, settings);}, dur1);
//        RunBenchmark("rank 2", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 2, i, settings);}, dur2);
//        RunBenchmark("rank 3", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 3, i, settings);}, dur3);
//        RunBenchmark("rank 4", signals, [&](Input in, int i){return RecursiveSparseFFT(DataSignal(info, in.data()), info, sparsity, 4, i, settings);}, dur4);
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
    return 0;
}

