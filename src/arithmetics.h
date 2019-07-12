#pragma once

#include "complex"
#include "math.h"


using complex_t  = std::complex<double>;

#define _USE_MATH_DEFINES
const double PI = M_PI;
const complex_t I = complex_t(0, 1);
const double EPS = 1e-9;

complex_t CalcKernel(double power, double base) {
    double phi = 2 * PI * power / base;
    return {cos(phi), sin(phi)};
}

bool NonZero(complex_t value) {
    return abs(value) > EPS;
}

int CalcLog(int64_t v) {
    int res = 0;
    while (v > 1) {
        v /= 2;
        ++res;
    }
    return res;
}