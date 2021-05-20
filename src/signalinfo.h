#include "arithmetics.h"

class SignalInfo {
public:
    static const int MAX_DIM = 6;

    SignalInfo(std::vector<int> signal_width):
            dimensions_(signal_width.size()) {
        signal_size_ = 1;
        for (int i = 0; i < dimensions_; ++i) {
            signal_width_[i] = signal_width[i];
            signal_size_ *= signal_width[i];
            log_signal_width_[i] = CalcLog(signal_width[i]);
        }
    }

    explicit SignalInfo(int dimensions, int signal_width): SignalInfo(std::vector<int>(dimensions, signal_width)) {}

    int Dimensions() const {
        return dimensions_;
    }

    int SignalWidth(int pos) const {
        return signal_width_[pos];
    }

    std::vector<int> GetAllDimensions() const {
        std::vector<int> res(dimensions_);
        for (int i = 0; i < dimensions_; ++i) {
            res[i] = signal_width_[i];
        }
        return res;
    }

    int MaxLogWidth() const {
        int max = 0;
        for (int i = 0; i < dimensions_; ++i) {
            if (max < log_signal_width_[i]) {
                max = log_signal_width_[i];
            }
        }
        return max;
    }

    bool IsSmallSignalWidth() const {
        for (int i = 0; i < dimensions_; ++i) {
           if (signal_width_[i] > SMALL_SIGNAL_WIDTH) {
               return false;
           }
        }
        return true;
    }

    int LogSignalWidth(int pos) const {
        return log_signal_width_[pos];
    }

    int LogSignalSize() const {
        int sum = 0;
        for (int i = 0; i < dimensions_; ++i) {
            sum += log_signal_width_[i];
        }
        return sum;
    }

    int64_t SignalSize() const {
        return signal_size_;
    }

    bool operator==(const SignalInfo& other) const {
        if (dimensions_ != other.dimensions_) {
            return false;
        }
        for (int i = 0; i < dimensions_; ++i) {
            if (signal_width_[i] != other.signal_width_[i]) {
                return false;
            }
        }
        return true;
    }

private:
    int dimensions_;
    int signal_width_[MAX_DIM];
    int64_t signal_size_;
    int log_signal_width_[MAX_DIM];
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
            res = (res << info_.LogSignalWidth(i)) | indices_[i];
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
        res.indices_[index] += value + info_.SignalWidth(index);
        res.indices_[index] %= info_.SignalWidth(index);
        return res;
    }

    void StoreDifference(const Key& a, const Key& b) {
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = (a.indices_[i] - b.indices_[i] + info_.SignalWidth(i)) & (info_.SignalWidth(i) - 1);
        }
    }

    Key operator-() const {
        Key res(info_);
        for (int i = 0; i < info_.Dimensions(); ++i) {
            res.indices_[i] = (-indices_[i] + info_.SignalWidth(i)) % info_.SignalWidth(i);
        }
        return res;
    }

//    int64_t operator*(const Key& key) const {
//        int64_t result = 0;
//        for (int i = 0; i < info_.Dimensions(); ++i) {
//            result += key[i] * indices_[i];
//        }
//        return result;
//    }

    // use with MaxLogWidth to get correct exp powers
    int64_t GetProduct(const Key& key) const {
        int64_t result = 0;
        int max_log = info_.MaxLogWidth();
        for (int i = 0; i < info_.Dimensions(); ++i) {
            result += (key[i] * indices_[i]) * (1 << (max_log - info_.LogSignalWidth(i)));
        }
        return result;
    }

    void SetFromFlatten(int64_t flat) {
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = flat & (info_.SignalWidth(i) - 1);
            flat >>= info_.LogSignalWidth(i);
        }
    }

private:

    SignalInfo info_;
    int64_t indices_[MAX_DIM];
};