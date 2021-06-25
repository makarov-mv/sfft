#include "arithmetics.h"
#include <x86intrin.h>
#include <iostream>


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
    const static int MAX_DIM = 4;
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

    explicit Key(const SignalInfo& info, std::vector<int32_t> key) : info_(info) {
        assert(info.Dimensions() == static_cast<int>(key.size()));
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = key[i];
        }
    }

    explicit Key(const SignalInfo& info, const std::initializer_list<int32_t>& key) : info_(info) {
        assert(info.Dimensions() == static_cast<int>(key.size()));
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = data(key)[i];
        }
    }

    explicit Key(const SignalInfo& info, int32_t flatten) : info_(info) {
        SetZero();
        SetFromFlatten(flatten);
    }

    void SetZero() {
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = 0;
        }
    }

    // leftmost dimension is highest in the tree
    const int32_t& operator[](int index) const {
        return indices_[index];
    }

    int32_t& operator[](int index) {
        return indices_[index];
    }

    // in flattened form, least significant dimension is highest
    int32_t Flatten() const {
        int32_t res = 0;
        #pragma omp simd
        for (int i = info_.Dimensions() - 1; i >= 0; --i) {
            res = (res << info_.LogSignalWidth()) | indices_[i];
        }
        return res;
        /*
        if (info_.Dimensions() <= 1){
            return indices_[0];
        } else {
            __m128i _key1;
            __m128 _key1f, _key2f, _result;
            
            _key1 = _mm_load_si128 ((__m128i *)&indices_[0]);
            
            _key1f = _mm_cvtepi32_ps (_key1);
            _key2f = _mm_setr_ps (1, 1<< info_.LogSignalWidth(), 1<< 2*info_.LogSignalWidth(), 1<< 3*info_.LogSignalWidth());

            _result = _mm_dp_ps (_key1f, _key2f, (1<<8) - 1);
            
            return _mm_cvtt_ss2si (_result);
        }
        */
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

    Key IncreaseAt(int index, int32_t value) const {
        Key res(*this);
        res.indices_[index] += value + info_.SignalWidth();
        res.indices_[index] %= info_.SignalWidth();
        return res;
    }

    void StoreDifference(const Key& a, const Key& b) {
        int32_t mod_ = info_.SignalWidth() - 1;
        /*
        #pragma omp simd
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = (a.indices_[i] - b.indices_[i] + info_.SignalWidth()) & mod_;
        }
         */
        if (info_.Dimensions() <= 1){
            indices_[0] = (a.indices_[0] - b.indices_[0] + info_.SignalWidth()) & mod_;
        } else {
            __m128i _keya, _keyb, _keyc, _result, _mod;
        
            _keya = _mm_load_si128 ((__m128i *)&a[0]);
            _keyb = _mm_load_si128 ((__m128i *)&b[0]);
        
            _result = _mm_sub_epi32 (_keya, _keyb);
            
            _keyc = _mm_set1_epi32 ((int32_t) info_.SignalWidth());
            _mod = _mm_set1_epi32 ( mod_ );
        
            _result = _mm_add_epi32 (_result, _keyc);
            _result = _mm_and_si128 (_result,_mod);
        
            _mm_store_si128 ((__m128i *)&indices_[0], _result);
        }
        
    }

    Key operator-() const {
        Key res(info_);
        int32_t mod_ = info_.SignalWidth() - 1;
        /*
        #pragma omp simd
        for (int i = 0; i < info_.Dimensions(); ++i) {
            res.indices_[i] = (-indices_[i] + info_.SignalWidth()) & mod_;
        }
        */
        if (info_.Dimensions() <= 1){
            res.indices_[0] = (-indices_[0] + info_.SignalWidth()) & mod_;
        } else {
            __m128i _key, _width, _result, _mod;
            
            _width = _mm_set1_epi32 ((int32_t) info_.SignalWidth());
            
            _key = _mm_load_si128 ((__m128i *)&indices_[0]);
        
            _result = _mm_sub_epi32 (_width, _key);
            
            _mod = _mm_set1_epi32 ( mod_ );
            
            _result = _mm_and_si128 (_result,_mod);
                    
            _mm_store_si128 ((__m128i *)&res[0], _result);
        }
                
        return res;
    }

    int32_t operator*(const Key& key) const {
        /*
        int32_t result = 0;
        #pragma omp simd
        for (int i = 0; i < info_.Dimensions(); ++i) {
                result += key[i] * indices_[i];
            }
        return result;
        */
        if (info_.Dimensions() <= 1){
            return key[0] * indices_[0];
        } else {
            
            __m128i _key1, _key2;
            __m128 _key1f, _key2f, _result;
            
            _key1 = _mm_load_si128 ((__m128i *)&indices_[0]);
            _key2 = _mm_load_si128 ((__m128i *)&key[0]);
            
            _key1f = _mm_cvtepi32_ps (_key1);
            _key2f = _mm_cvtepi32_ps (_key2);

            _result = _mm_dp_ps (_key1f, _key2f, (1<<8) - 1);
            
            return _mm_cvtt_ss2si (_result);
            
            /*
            #pragma omp simd
            return key[0] * indices_[0] + key[1] * indices_[1] + key[2] * indices_[2] + key[3] * indices_[3];
            */
            
        }
        
    }

    void SetFromFlatten(int32_t flat) {
        int32_t mod_ = info_.SignalWidth() - 1;
        for (int i = 0; i < info_.Dimensions(); ++i) {
            indices_[i] = flat & mod_;
            flat >>= info_.LogSignalWidth();
        }
    }

private:
    SignalInfo info_;
    int32_t indices_[MAX_DIM] = {0};
};
