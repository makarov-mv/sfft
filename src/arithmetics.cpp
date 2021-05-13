#include "arithmetics.h"
#include "cassert"

double cos_table[SMALL_SIGNAL_WIDTH];
double sin_table[SMALL_SIGNAL_WIDTH];

void PrepareCosSinTables(int width) {
    assert(width <= SMALL_SIGNAL_WIDTH);
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