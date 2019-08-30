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
    3092,
    24521,
    245324,
    3905806,
    74959133,
]
rank1 = [
    39262852,
    76248189,
    108173567,
    187097804,
    362859257,
]
rank2 = [
    21374912,
    25031167,
    28713843,
    38705956,
    46126198,
]
rank3 = [
    19939059,
    24825588,
    31005542,
    47313196,
    58460593,
]
rank4 = [
    30201766,
    41882062,
    53713228,
    108221484,
    116224935,
]
rand_phase = [
    4082713,
    7253853,
    9653051,
    27680461,
    26012096,
]

graphs = {
    'fftw': fftw,
    'disFT, $n^3$': rank1,
    r'disFT, $n^{2\frac{1}{2}}$': rank2,
    r'disFT, $n^{2\frac{1}{3}}$': rank3,
    r'disFT, $n^{2\frac{1}{4}}$': rank4,
    r'disFT, random phase' : rand_phase,
}

for gr in graphs.items():
    plt.plot(p, np.array(gr[1]) / milsec, label=gr[0])

plt.legend()
plt.title("Comparison of algorithms, 3 dimensions,\n sparsity=27, without comb filter, random support")
plt.xlabel(r"$\log_2 n$")
plt.ylabel("ms")
plt.grid()
plt.yscale('log')
#plt.show()
plt.savefig('without_comb.png', dpi=300)
