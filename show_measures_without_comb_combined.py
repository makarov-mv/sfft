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
    2817,
    24230,
    279066,
    3860319,
    57886542,
]
rank1 = [
    23047176,
    37546542,
    55125746,
    91882079,
    184553587,
]
rank2 = [
    25157936,
    28470570,
    32722411,
    54537677,
    58752167,
]
rank3 = [
    22894221,
    27141661,
    32155919,
    72939077,
    58776871,
]
rank4 = [
    29923491,
    40406856,
    45376116,
    122315533,
    87406394,
]
rand_phase = [
    3194544,
    7908017,
    9427138,
    18267691,
    68042355,
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
plt.title("Comparison of algorithms, 3 dimensions,\n sparsity=27, without comb filter, random support + Dirac comb")
plt.xlabel(r"$\log_2 n$")
plt.ylabel("ms")
plt.grid()
plt.yscale('log')
#plt.show()
plt.savefig('without_comb_combined.png', dpi=300)
