#pragma once

#include "complex"
#include "math.h"


using complex_t  = std::complex<double>;

#define _USE_MATH_DEFINES
const double PI = M_PI;
const complex_t I = complex_t(0, 1);
const double EPS = 1e-10;

complex_t CalcKernel(double power, double base) {
    return std::exp(2 * PI * I * power / base);
}

bool NonZero(complex_t value) {
    return abs(value) > EPS;
}