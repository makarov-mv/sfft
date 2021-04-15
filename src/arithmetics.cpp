#include "arithmetics.h"
#include "assert.h"

double cos_table[1 << 10];
double sin_table[1 << 10];

void PrepareCosSinTables(int width) {
    assert(width <= (1 << 10));
    for (int i = 0; i < width; ++i) {
        cos_table[i] = cos((2 * PI * i) / width);
        sin_table[i] = sin((2 * PI * i) / width);
    }
}

double GetTableCos(int64_t n, int width) {
    return cos_table[n & (width - 1)];
}
double GetTableSin(int64_t n, int width) {
    return sin_table[n & (width - 1)];
}