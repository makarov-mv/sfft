#include "catch.hpp"
#include "projection_recovery.h"
#include "utility_test_fftw.h"
#include "signal.h"

bool RunPFTTest(const std::vector<complex_t>& out, const SignalInfo& info, int extra_measures = 6) {
    assert(info.SignalSize() == static_cast<int64_t>(out.size()));
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    auto result = GetSignalFromMap(ProjectionFT(x, info, extra_measures), info);
    return std::equal(out.begin(), out.end(), result.begin(), CheckEqual);
}

TEST_CASE("PFT 8") {
    SignalInfo info{2, 8};
    std::vector<complex_t> out(info.SignalSize());
    out[0] = 1;
    out[1] = 2;
    REQUIRE(RunPFTTest(out, info, 6));
}

TEST_CASE("PFT 8 shifted") {
    SignalInfo info{2, 8};
    std::vector<complex_t> out(info.SignalSize());
    out[info.SignalWidth()] = 1;
    out[info.SignalWidth()*2 + 1] = 2;
    REQUIRE(RunPFTTest(out, info, 6));
}

TEST_CASE("PFT 16") {
    SignalInfo info{2, 16};
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalWidth(); ++i) {
        out[i] = i;
    }
    REQUIRE(RunPFTTest(out, info));
}

TEST_CASE("PFT diag") {
    SignalInfo info{2, 16};
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalWidth(); ++i) {
        out[i + i * info.SignalWidth()] = i + 1;
    }
    REQUIRE(RunPFTTest(out, info));
}

TEST_CASE("PFT 256") {
    SignalInfo info{3, 128};
    std::vector<complex_t> out(info.SignalSize());
    for (int i = 0; i < info.SignalWidth(); ++i) {
        out[i] = i;
        out[info.SignalWidth() * i] = -i;
    }
    REQUIRE(RunPFTTest(out, info));
}

TEST_CASE("PFT zero recovery") {
    SignalInfo info{2, 16};
    std::vector<complex_t> out(info.SignalSize());
    out[0] = 1;
    out[1] = 1;
    out[info.SignalWidth()] = 1;
    out[info.SignalWidth() + 1] = 1;
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    auto res = ProjectionFT(x, info, 6);
    REQUIRE(res.empty());
}

TEST_CASE("PFT partial recovery") {
    SignalInfo info{2, 16};
    std::vector<complex_t> out(info.SignalSize());
    out[0] = 1;
    out[1] = 1;
    out[info.SignalWidth()] = 1;
    out[info.SignalWidth() + 1] = 1;
    out[info.SignalWidth()*3 + 1] = -1;
    out[info.SignalWidth() + 4] = -2;
    auto runner = FFTWRunner(info, FFTW_BACKWARD);
    auto in = runner.Run(out);
    auto x = DataSignal(info, in.data());
    auto res = ProjectionFT(x, info, 6);

    REQUIRE(res.size() == 2);
    Key k(info);
    k.SetFromFlatten(info.SignalWidth()*3 + 1);
    REQUIRE(res.find(k) != res.end());
    REQUIRE(CheckEqual(res.at(k), -1));
    k.SetFromFlatten(info.SignalWidth() + 4);
    REQUIRE(res.find(k) != res.end());
    REQUIRE(CheckEqual(res.at(k), -2));
}
