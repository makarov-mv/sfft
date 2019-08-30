import matplotlib.pyplot as plt
import numpy as np

milsec = 1e6

p = [
    3,
    4,
    5,
    6,
    7,
]
fftw = [
    2612,
    23903,
    244691,
    3574356,
    67041862,
]
rank1 = [
    1826934,
    1934921,
    2854730,
    4226215,
    7075255,
]
rank2 = [
    4541388,
    4375681,
    4895753,
    6393661,
    6523560,
]
rank3 = [
    6253416,
    6112755,
    6750564,
    6844653,
    7752159,
]
rank4 = [
    8235431,
    9018661,
    10772709,
    13318868,
    12748138,
]
rand_phase = [
    1729185,
    2187956,
    3371758,
    5021928,
    8959761,
]

graphs = {
    'fftw': fftw,
    'disFT, $n^3$': rank1,
    r'disFT, $n^{2\frac{1}{2}}$': rank2,
    r'disFT, $n^{2\frac{1}{3}}$': rank3,
    #r'disFT, $n^{2\frac{1}{4}}$': rank4,
    r'disFT, random phase' : rand_phase,
}

for gr in graphs.items():
    plt.plot(p, np.array(gr[1]) / milsec, label=gr[0])

plt.legend()
plt.title("Comparison of algorithms, 3 dimensions, prob=90%\n sparsity=32, using comb filter, random support + Dirac comb")
plt.xlabel(r"$\log_2 n$")
plt.ylabel("ms")
plt.grid()
plt.yscale('log')
#plt.show()
plt.savefig('with_comb_combined.png', dpi=300)
