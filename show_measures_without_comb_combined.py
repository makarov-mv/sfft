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
    2514,
    22979,
    225030,
    3664961,
    53502254,
]
rank1 = [
    12001927,
    20316605,
    31955925,
    49560701,
    88525859,
]
rank2 = [
    13976127,
    15905738,
    19881683,
    23693780,
    28744511,
]
rank3 = [
    15203625,
    18066933,
    23497032,
    28070731,
    32884768,
]
rank4 = [
    22781848,
    28820350,
    38534212,
    46265741,
    55707130,
]
rand_phase = [
    4384830,
    8182715,
    17019835,
    18383747,
    29895710,
]

graphs = {
    'fftw': fftw,
    'disFT, $n^3$': rank1,
    r'disFT, $n^{2\frac{1}{2}}$': rank2,
    r'disFT, $n^{2\frac{1}{3}}$': rank3,
    #r'disFT, $n^{2\frac{1}{4}}$': rank4,
    r'disFT, random phase' : rand_phase,
}

ax = plt.axes()
ax.set_xticks(p)

for gr in graphs.items():
    plt.plot(p, np.array(gr[1]) / milsec, label=gr[0])

plt.legend()
plt.title("Comparison of algorithms, 3 dimensions,\n sparsity=32, without comb filter, random support + Dirac comb")
plt.xlabel(r"$\log_2 n$")
plt.ylabel("ms")
plt.grid()
plt.yscale('log')
#plt.show()
plt.savefig('without_comb_combined.png', dpi=300)
