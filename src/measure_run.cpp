#include "iostream"
#include "chrono"
#include "map"
#include "disfft.h"
#include "utility_test_fftw.h"
#include "utility_test.h"

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
        out[i] = 1;//CalcKernel(i, info.SignalWidth());
    }
    return out;
}

template <class Generator>
std::vector<complex_t> GenRandomSupportWithOvertones(const SignalInfo& info, int64_t sparsity, Generator& gen) {
    if (info.Dimensions() < 2) {
        return GenRandomSupport(info, sparsity, gen);
    }
    std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);
    std::vector<complex_t> out(info.SignalSize());
    int maxd = std::min(4, info.Dimensions() + 1);
    for (int i = 0; i < sparsity / maxd; ++i) {
        auto pos = dist(gen);
        out[pos] += 1;
        Key overtone(info, pos);
        for (int j = 0; j < maxd - 1; ++j) {
            out[overtone.IncreaseAt(j, info.SignalWidth() / 2).Flatten()] += 0.5;
        }
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

class SignalGenerator {
public:
    virtual ~SignalGenerator() = default;
    virtual std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) = 0;
};

class RandomSignalGenerator : public SignalGenerator {
public:
    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
        return GenRandomSupport(info, sparsity, gen);
    }
};

class RandomSignalGeneratorWithOvertones : public SignalGenerator {
public:
    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
        return GenRandomSupportWithOvertones(info, sparsity, gen);
    }
};

class DiracCombGenerator : public SignalGenerator {
public:
    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64&) override {
        return GenDiracComb(info, sparsity);
    }
};

class RandomCombGenerator : public SignalGenerator {
public:
    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
        return GenCombined(info, sparsity, gen);
    }
};

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

class Algorithm {
public:
    virtual ~Algorithm() = default;
    virtual void Prepare(const SignalInfo& info, int sparsity, TransformSettings settings) = 0;
    virtual void Run(const DataSignal& signal, int seed, std::chrono::nanoseconds& dur, const std::vector<complex_t>& FTsignal) = 0;
};

auto RunBenchmark(const std::string&, const DataSignal& signals, Algorithm& alg, std::vector<std::chrono::nanoseconds>& res, const std::vector<complex_t>& FTsignals, int samples) {
    //using clock = std::chrono::system_clock;
//    std::cout << name << ": ";
    
    for (int i = 0; i < samples; ++i) {
        std::chrono::nanoseconds dur;
        alg.Run(signals, i, dur, FTsignals);
        //dur /= signals.size();
        res.push_back(dur);
    }
//    PrintDur(dur);
//    std::cout << ", ";
}

//class Benchmark {
//public:
//    virtual ~Benchmark() = default;
//    virtual int Sparsity() = 0;
//    virtual std::vector<TransformSettings> PrepareSettings() = 0;
//    virtual std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) = 0;
//    virtual int Start() = 0;
//    virtual int End() = 0;
//    virtual SignalInfo GetInfo(int p) = 0;
//};
//
//class BenchmarkRandomSupport : public Benchmark {
//public:
//    BenchmarkRandomSupport(bool use_comb) : use_comb_(use_comb) {}
//
//    int Sparsity() override {
//        return 27;
//    }
//
//    std::vector<TransformSettings> PrepareSettings() override {
//        std::vector<TransformSettings> res(5);
//        TransformSettings settings;
//        settings.use_comb = use_comb_;
//        settings.zero_test_koef = 0.03;
//        for (int i = 0; i < 4; ++i) {
//            res[i] = settings;
//        }
//        settings.zero_test_koef = 1;
//        settings.random_phase_sparsity_koef = 1;
//        settings.assume_random_phase = true;
//        res[4] = settings;
//        return res;
//    }
//
//    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
//        return GenRandomSupport(info, sparsity, gen);
//    }
//
//    int Start() override {
//        return 3;
//    }
//
//    int End() override {
//        return 7;
//    }
//
//    SignalInfo GetInfo(int p) override {
//        return {3, 1 << p};
//    }
//
//private:
//    bool use_comb_;
//};
//
//class BenchmarkCombinedSupport : public Benchmark {
//public:
//    BenchmarkCombinedSupport(bool use_comb) : use_comb_(use_comb) {}
//
//    int Sparsity() override {
//        return 32;
//    }
//
//    std::vector<TransformSettings> PrepareSettings() override {
//        std::vector<TransformSettings> res(5);
//        TransformSettings settings;
//        settings.use_comb = use_comb_;
//        settings.zero_test_koef = 0.5;
//        for (int i = 0; i < 4; ++i) {
//            res[i] = settings;
//        }
//        settings.zero_test_koef = 1;
//        settings.random_phase_sparsity_koef = 4;
//        settings.assume_random_phase = true;
//        res[4] = settings;
//        return res;
//    }
//
//    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
//        return GenCombined(info, sparsity, gen);
//    }
//
//    int Start() override {
//        return 3;
//    }
//
//    int End() override {
//        return 7;
//    }
//
//    SignalInfo GetInfo(int p) override {
//        return {3, 1 << p};
//    }
//
//private:
//    bool use_comb_;
//};
//
//class Benchmark4DCombSupportWithComb : public Benchmark {
//public:
//    Benchmark4DCombSupportWithComb(bool use_comb) : use_comb_(use_comb) {}
//
//    int Sparsity() override {
//        return comb_sparsity + rand_supp;
//    }
//
//    std::vector<TransformSettings> PrepareSettings() override {
//        std::vector<TransformSettings> res(5);
//        TransformSettings settings;
//        settings.use_comb = use_comb_;
//        settings.zero_test_koef = 0.5;
//        for (int i = 0; i < 4; ++i) {
//            res[i] = settings;
//        }
//        settings.zero_test_koef = 1;
//        settings.random_phase_sparsity_koef = 10;
//        settings.assume_random_phase = true;
//        res[4] = settings;
//        return res;
//    }
//
//    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen) override {
//        assert(sparsity == comb_sparsity + rand_supp);
//        auto out = Gen4DComb(info, comb_sparsity);
//        std::uniform_int_distribution<int64_t> dist(0, info.SignalSize() - 1);
//        for (int i = 0; i < rand_supp; ++i) {
//            out[dist(gen)] = 1;
//        }
//        return out;
//    }
//
//    int Start() override {
//        return 3;
//    }
//
//    int End() override {
//        return 5;
//    }
//
//    SignalInfo GetInfo(int p) override {
//        return {4, 1 << p};
//    }
//
//private:
//    bool use_comb_;
//    const int comb_sparsity = 64;
//    const int rand_supp = 64;
//};

class FFTWAlgorithm : public Algorithm {
public:
    void Prepare(const SignalInfo& info, int, TransformSettings) override {
        runner_.emplace(info, FFTW_FORWARD);
    }

    void Run(const DataSignal& signal, int, std::chrono::nanoseconds& dur, const std::vector<complex_t>& FTsignal) override {
        auto start = std::chrono::system_clock::now();
        auto result = runner_.value().Run(signal.Data());
        dur = std::chrono::system_clock::now() - start;
        if (!std::equal(FTsignal.begin(), FTsignal.end(), result.begin(), CheckEqual)){
            dur = std::chrono::nanoseconds(-1);
        }
    }

private:
    std::optional<FFTWRunner> runner_;
};

class RecursiveAlgorithm : public Algorithm {
public:
    RecursiveAlgorithm(int rank) : rank_(rank) {}

    void Prepare(const SignalInfo& info, int sparsity, TransformSettings settings) override {
        if (settings.assume_random_phase) {
            throw std::runtime_error("incorrect settings for recursive");
        }
        info_ = info;
        sparsity_ = sparsity;
        settings_ = settings;
    }

    void Run(const DataSignal& signal, int seed, std::chrono::nanoseconds& dur, const std::vector<complex_t>& FTsignal) override {
        auto start = std::chrono::system_clock::now();
        FrequencyMap frequency = RecursiveSparseFFT(signal, info_.value(), sparsity_, rank_, seed, settings_);
        dur = std::chrono::system_clock::now() - start;
        
        auto result = GetSignalFromMap(frequency, info_.value());
        if (!std::equal(FTsignal.begin(), FTsignal.end(), result.begin(), CheckEqual)){
            dur = std::chrono::nanoseconds(-1);
        }
    }

private:
    int rank_;
    TransformSettings settings_;
    std::optional<SignalInfo> info_;
    int sparsity_;
};

class RandomPhaseAlgorithm : public Algorithm {
public:
    void Prepare(const SignalInfo& info, int sparsity, TransformSettings settings) override {
        if (!settings.assume_random_phase) {
            throw std::runtime_error("incorrect settings for random_phase");
        }
        info_ = info;
        sparsity_ = sparsity;
        settings_ = settings;
    }

    void Run(const DataSignal& signal, int seed, std::chrono::nanoseconds& dur, const std::vector<complex_t>& FTsignal) override {
        auto start = std::chrono::system_clock::now();
        FrequencyMap frequency = RecursiveSparseFFT(signal, info_.value(), sparsity_, 1, seed, settings_);
        dur = std::chrono::system_clock::now() - start;
        
        auto result = GetSignalFromMap(frequency, info_.value());
        if (!std::equal(FTsignal.begin(), FTsignal.end(), result.begin(), CheckEqual)){
            dur = std::chrono::nanoseconds(-1);
        }
    }

private:
    TransformSettings settings_;
    std::optional<SignalInfo> info_;
    int sparsity_;
};

class StreamBenchmark {
public:
    StreamBenchmark(std::istream& in) {
        CheckInput("samples:", in);
        in >> samples_;
        int alg_cnt;
        in >> alg_cnt;
        for (int i = 0; i < alg_cnt; ++i) {
            std::string name;
            in >> name;
            if (name == "recursive") {
                int rank;
                in >> rank;
                name.push_back(' ');
                name.append(std::to_string(rank));
                algs_[name] = std::make_unique<RecursiveAlgorithm>(rank);
            } else if (name == "fftw") {
                algs_[name] = std::make_unique<FFTWAlgorithm>();
            } else if (name == "random_phase") {
                algs_[name] = std::make_unique<RandomPhaseAlgorithm>();
            } else {
                throw std::runtime_error("no such alg");
            }
            settings_[name] = {};
        }
        int iters;
        in >> iters;
        for (int it_num = 0; it_num < iters; ++it_num) {
            CheckInput("id:", in);
            int p;
            in >> p;
            p_.push_back(p);

            CheckInput("info:", in);
            int64_t n;
            int d;
            in >> d >> n;
            info_.emplace_back(d, 1ll << n);

            CheckInput("sparsity:", in);
            int sparsity;
            in >> sparsity;
            sparsity_.push_back(sparsity);

            CheckInput("signal:", in);
            std::string signal_type;
            in >> signal_type;
            if (signal_type == "random") {
                generator_.push_back(std::make_unique<RandomSignalGeneratorWithOvertones>());
            } else if (signal_type == "comb") {
                generator_.push_back(std::make_unique<DiracCombGenerator>());
            } else if (signal_type == "combined") {
                generator_.push_back(std::make_unique<RandomCombGenerator>());
            } else {
                throw std::runtime_error("unknown signal type");
            }
            for (int i = 0; i < alg_cnt; ++i) {
                std::string name;
                in >> name;
                if (name == "recursive") {
                    int rank;
                    in >> rank;
                    name.push_back(' ');
                    name.append(std::to_string(rank));
                }
                TransformSettings settings;
                int field_cnt;
                in >> field_cnt;
                for (int j = 0; j < field_cnt; ++j) {
                    std::string field;
                    in >> field;
                    if (field == "use_comb") {
                        in >> settings.use_comb;
                    } else if (field == "zero_test_coef") {
                        in >> settings.zero_test_koef;
                    } else if (field == "use_projection_recovery") {
                        in >> settings.use_projection_recovery;
                    } else if (field == "assume_random_phase") {
                        in >> settings.assume_random_phase;
                    } else {
                        throw std::runtime_error("no such settings field: " + field);
                    }
                }
                settings_.at(name).push_back(settings);
            }
        }
    }

    int Sparsity(int iter) {
        return sparsity_[iter];
    }

    TransformSettings PrepareSettings(const std::string& alg, int iter) {
        return settings_.at(alg).at(iter);
    }

    std::vector<complex_t> GenSignal(const SignalInfo& info, int64_t sparsity, std::mt19937_64& gen, int iter) {
        return generator_.at(iter)->GenSignal(info, sparsity, gen);
    }

    int Start() {
        return 0;
    }

    int End() {
        return p_.size() - 1;
    }

    int GetP(int iter) {
        return p_.at(iter);
    }

    SignalInfo GetInfo(int iter) {
        return info_.at(iter);
    }

    std::map<std::string, std::unique_ptr<Algorithm>>& GetAlgs() {
        return algs_;
    }

    int GetSamples() {
        return samples_;
    }

private:
    void CheckInput(const std::string& string, std::istream& in) {
        std::string value;
        in >> value;
        if (value != string) {
            throw std::runtime_error(string + " != " + value);
        }
    }

    int samples_;
    std::vector<int> p_;
    std::vector<SignalInfo> info_;
    std::vector<int> sparsity_;
    std::vector<std::unique_ptr<SignalGenerator> > generator_;
    std::unordered_map<std::string, std::vector<TransformSettings> > settings_;
    std::map<std::string, std::unique_ptr<Algorithm>> algs_;
};

int main() {
    std::mt19937_64 gen(123423);
    std::vector<int64_t> npow;
    std::map<std::string, std::vector<std::chrono::nanoseconds>> dur;
    StreamBenchmark bench(std::cin);
    for (auto& w : bench.GetAlgs()) {
        dur[w.first] = {};
    }
    for (int64_t iter = bench.Start(); iter <= bench.End(); ++iter) {
        int p = bench.GetP(iter);
        SignalInfo info = bench.GetInfo(iter);
        const int64_t sparsity = bench.Sparsity(iter);
        const int samples = bench.GetSamples();

        
        auto runner = FFTWRunner(info, FFTW_BACKWARD);
        auto out = bench.GenSignal(info, sparsity, gen, iter);
        //for (int i = 0; i < samples; ++i) {
        std::vector<complex_t> FTsignals = out;
        std::vector<complex_t> signals = runner.Run(out);
        DataSignal datasignals(info, signals.data());
        //}
        
        npow.push_back(p);

        for (auto& alg : bench.GetAlgs()) {
            auto settings = bench.PrepareSettings(alg.first, iter);
//            std::cout << "p = " << p << ", ";
            alg.second->Prepare(info, sparsity, settings);
            RunBenchmark(alg.first, datasignals, *alg.second, dur.at(alg.first), FTsignals, samples);

//            std::cout << std::endl;
        }
    }
//    std::cout << "p = [\n";
    for (auto p: npow) {
        std::cout << p << " ";
    }
    std::cout << std::endl;
//    std::cout << "]\n";

    for (auto& alg : bench.GetAlgs()) {
        std::cout << alg.first << std::endl;
        for (auto p: dur.at(alg.first)) {
            std::cout << p.count() << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}

