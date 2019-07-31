import matplotlib.pyplot as plt

graphs = {
    'fftw': fftw,
    'rank 1': rank1,
    'rank 2': rank2,
    'rank 3': rank3,
}

for gr in graphs.items():
    plt.plot(p, gr[1], label=gr[0])
plt.legend()
plt.xlabel("$\log_2 n$")
plt.ylabel("runtime, $\mu s$")
plt.yscale('log')
plt.show()
