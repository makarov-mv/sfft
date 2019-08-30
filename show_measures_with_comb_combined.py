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
    2598,
    23235,
    260472,
    3800493,
    53884409,
]
rank1 = [
    6051085,
    10481633,
    17496975,
    22027618,
    32985893,
]
rank2 = [
    20241096,
    20262626,
    24204012,
    26450607,
    31337033,
]
rank3 = [
    17568775,
    18372300,
    22263247,
    25771086,
    29690744,
]
rank4 = [
    21155371,
    25178983,
    31189802,
    37479368,
    39469112,
]
rand_phase = [
    1807087,
    3220542,
    4158675,
    6111749,
    7367960,
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
plt.title("Comparison of algorithms, 3 dimensions, prob=90%\n sparsity=32, using comb filter, random support + Dirac comb")
plt.xlabel(r"$\log_2 n$")
plt.ylabel("ms")
plt.grid()
plt.yscale('log')
#plt.show()
plt.savefig('with_comb_combined.png', dpi=300)
