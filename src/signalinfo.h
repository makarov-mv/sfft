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
        return signal_width_ <= SMALL_SIGNAL_WIDTH;
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
    const static int MAX_DIM = 6;
    explicit Key(const SignalInfo& info) : info_(info) {
        SetZero();
    }

    Key(const Key& other) : info_(other.info_) {
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = other.indices_[i];
        }
    }

    Key& operator=(const Key& other) {
        info_ = other.info_;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = other.indices_[i];
        }
        return *this;
    }

    explicit Key(const SignalInfo& info, std::vector<int64_t> key) : info_(info) {
        assert(info.Dimensions() == static_cast<int>(key.size()));
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = key[i];
        }
    }

    explicit Key(const SignalInfo& info, const std::initializer_list<int64_t>& key) : info_(info) {
        assert(info.Dimensions() == static_cast<int>(key.size()));
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = data(key)[i];
        }
    }

    explicit Key(const SignalInfo& info, int64_t flatten) : info_(info) {
        SetZero();
        SetFromFlatten(flatten);
    }

    void SetZero() {
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = 0;
        }
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
        for (int i = 0; i < info_.Dimensions(); ++i) {
            if (indices_[i] != other.indices_[i]) {
                return false;
            }
        }
        return true;
    }

    SignalInfo GetSignalInfo() const {
        return info_;
    }

    Key IncreaseAt(int index, int64_t value) const {
        Key res(*this);
        res.indices_[index] += value + info_.SignalWidth();
        res.indices_[index] %= info_.SignalWidth();
        return res;
    }

    void StoreDifference(const Key& a, const Key& b) {
        int64_t mod = info_.SignalWidth() - 1;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = (a.indices_[i] - b.indices_[i] + info_.SignalWidth()) & mod;
        }
    }

    Key operator-() const {
        Key res(info_);
        for (int i = 0; i < info_.Dimensions(); ++i) {
            res.indices_[i] = (-indices_[i] + info_.SignalWidth()) % info_.SignalWidth();
        }
        return res;
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
    int64_t indices_[MAX_DIM];
};