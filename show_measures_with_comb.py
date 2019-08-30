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
    2798,
    24306,
    239436,
    4455979,
    68843906,
]
rank1 = [
    10904300,
    28991625,
    55053279,
    85511819,
    136437876,
]
rank2 = [
    12977353,
    14904801,
    12672936,
    15947939,
    17667265,
]
rank3 = [
    4869344,
    5883737,
    7054355,
    10499698,
    13622528,
]
rank4 = [
    3526270,
    7618888,
    10094225,
    15546406,
    22984291,
]
rand_phase = [
    596548,
    1190308,
    2167206,
    3212119,
    5176930,
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
