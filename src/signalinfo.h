#include "arithmetics.h"

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

    bool IsSmallSignalWidth() const {
        return signal_width_ <= (1 << 10);
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
    explicit Key(const SignalInfo& info) : info_(info), indices_(info.Dimensions(), 0) {
    }

    explicit Key(const SignalInfo& info, std::vector<int64_t> key) : info_(info), indices_(std::move(key)) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
    }

    explicit Key(const SignalInfo& info, const std::initializer_list<int64_t>& key) : info_(info), indices_(info.Dimensions(), 0) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
        std::copy(key.begin(), key.end(), indices_.begin());
    }

    explicit Key(const SignalInfo& info, int64_t flatten) : info_(info), indices_(info.Dimensions(), 0) {
        assert(info.Dimensions() == static_cast<int>(indices_.size()));
        SetFromFlatten(flatten);
    }

    void SetZero() {
        indices_.assign(indices_.size(), 0);
    }

    // leftmost dimension is highest in the tree
    const int64_t& operator[](int index) const {
        return indices_[index];
    }

    int64_t& operator[](int index) {
        return indices_[index];
    }

    // in flattened form, least significant dimension is highest
    int64_t Flatten() const {
        int64_t res = 0;
        for (int i = info_.Dimensions() - 1; i >= 0; --i) {
            res = (res << info_.LogSignalWidth()) | indices_[i];
        }
        return res;
    }

    bool operator==(const Key& other) const {
        assert(info_ == other.info_);
        return indices_ == other.indices_;
    }

    SignalInfo GetSignalInfo() const {
        return info_;
    }

    Key IncreaseAt(int index, int64_t value) const {
        std::vector<int64_t> new_indices(indices_);
        new_indices[index] += value + info_.SignalWidth();
        new_indices[index] %= info_.SignalWidth();
        return Key{info_, std::move(new_indices)};
    }

    void StoreDifference(const Key& a, const Key& b) {
        int64_t mod = info_.SignalWidth() - 1;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = (a.indices_[i] - b.indices_[i] + info_.SignalWidth()) & mod;
        }
    }

    Key operator-() const {
        std::vector<int64_t> new_indices(info_.Dimensions());
        for (int i = 0; i < info_.Dimensions(); ++i) {
            new_indices[i] = (-indices_[i] + info_.SignalWidth()) % info_.SignalWidth();
        }
        return Key{info_, std::move(new_indices)};
    }

    int64_t operator*(const Key& key) const {
        int64_t result = 0;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            result += key[i] * indices_[i];
        }
        return result;
    }

    void SetFromFlatten(int64_t flat) {
        int64_t mod = info_.SignalWidth() - 1;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = flat & mod;
            flat >>= info_.LogSignalWidth();
        }
    }

private:
    SignalInfo info_;
    std::vector<int64_t> indices_;
};