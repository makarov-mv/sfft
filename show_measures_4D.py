import matplotlib.pyplot as plt
import numpy as np

milsec = 1e6

p = [
    3,
    4,
    5,
]
fftw = [
    53064,
    1094238,
    25550446,
]
rank1 = [
    1094949230,
    1287785893,
    2043846232,
]
rank2 = [
    1098628514,
    1629149667,
    1728558092,
]
rank3 = [
    734742690,
    1453657673,
    1742293185,
]
rank4 = [
    975534803,
    1281012828,
    1136224549,
]
rand_phase = [
    180301222,
    205272859,
    230259089,
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
    plt.plot(p, np.array(gr[1]) / milsec, '.', label=gr[0])

plt.legend()
plt.title("Comparison of algorithms, 4D, random support + 2D comb \n sparsity=64 + 64, using comb filter for nonrecursive")
plt.xlabel(r"$\log_2 n$")
plt.ylabel("ms")
plt.grid()
plt.yscale('log')
#plt.show()
plt.savefig('with_comb_4D.png', dpi=300)
