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
    2838,
    26736,
    284525,
    3875702,
    53453656,
]
rank1 = [
    747935,
    2003619,
    3665210,
    5638472,
    6206940,
]
rank2 = [
    825795,
    1724183,
    2042071,
    2887528,
    2741810,
]
rank3 = [
    936691,
    2599587,
    3119049,
    4643492,
    4020754,
]
rank4 = [
    1238314,
    5237258,
    4613761,
    9648533,
    6924271,
]
rand_phase = [
    582755,
    1979812,
    2621370,
    5593132,
    4582832,
]

graphs = {
    'fftw': fftw,
    'disFT, $n^3$': rank1,
    r'disFT, $n^{2\frac{1}{2}}$': rank2,
    r'disFT, $n^{2\frac{1}{3}}$': rank3,
    r'disFT, $n^{2\frac{1}{4}}$': rank4,
    r'disFT, random phase' : rand_phase,
}

ax = plt.axes()
ax.set_xticks(p)

for gr in graphs.items():
    plt.plot(p, np.array(gr[1]) / milsec, label=gr[0])

p = [
    6,
    7,
]
sfft = [
    185327698,
    234982192,
]

plt.plot(p, np.array(sfft) / milsec, label='sfft')

plt.legend()
plt.title("Comparison of algorithms, 3 dimensions,\n sparsity=27, using comb filter, random_support")
plt.xlabel(r"$\log_2 n$")
plt.ylabel("ms")
plt.grid()
plt.yscale('log')
#plt.show()
plt.savefig('with_comb.png', dpi=300)
