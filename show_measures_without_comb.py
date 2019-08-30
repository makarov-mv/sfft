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
    2670,
    25942,
    263053,
    3702702,
    55129729,
]
rank1 = [
    3954893,
    7755332,
    11312490,
    15105190,
    19141955,
]
rank2 = [
    4360903,
    8921670,
    10592242,
    10726630,
    12558089,
]
rank3 = [
    9712045,
    18447490,
    21399027,
    26161412,
    28036076,
]
rank4 = [
    20814345,
    33410646,
    38975980,
    49523968,
    83686640,
]
rand_phase = [
    4526637,
    8362649,
    13632777,
    15887589,
    21287620,
]

ax = plt.axes()
ax.set_xticks(p)

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
