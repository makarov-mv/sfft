#pragma once

#include "complex"
#include "math.h"


using complex_t  = std::complex<double>;

#define _USE_MATH_DEFINES
const double PI = M_PI;
const complex_t I = complex_t(0, 1);
const double EPS = 1e-7;

inline complex_t CalcKernel(double power, double base) {
    double phi = 2 * PI * power / base;
    return {cos(phi), sin(phi)};
}

inline complex_t CalcKernelNormalized(double phi) {
    return {cos(phi), sin(phi)};
}

void PrepareCosSinTables(int width);

double GetTableCos(int64_t n, int width);
double GetTableSin(int64_t n, int width);

inline bool NonZero(complex_t value) {
    // two times faster than abs(value) > EPS
    return value.real() > EPS || value.real() < -EPS || value.imag() > EPS || value.imag() < -EPS;
}

inline int CalcLog(int64_t v) {
    int res = 0;
    while (v > 1) {
        v /= 2;
        ++res;
    }
    return res;
}