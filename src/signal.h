#pragma once
#include "arithmetics.h"
#include "signalinfo.h"

class Signal {
public:
    virtual ~Signal() = default;

    virtual complex_t ValueAtTime(const Key& key) const = 0;
};

class DataSignal: public Signal {
public:
    DataSignal(const SignalInfo& info, const complex_t* v):
            info_(info),
            values_(v) {
    }

    complex_t ValueAtTime(const Key& key) const override {
        return values_[key.Flatten()];
    }

    const complex_t* Data() const {
        return values_;
    }

private:
    const SignalInfo info_;
    const complex_t* values_;
};