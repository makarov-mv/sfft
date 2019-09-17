import argparse
import matplotlib.pyplot as plt
import numpy as np

def make_graph(name, title, p, algs, algs_names, data, extra_data_path=''):
    ax = plt.axes()
    ax.set_xticks(p)
    milsec = 1e6

    for alg in algs:
        plt.plot(p, np.array(data[alg]) / milsec, label=algs_names[alg])

    if extra_data_path:
        with open(extra_data_path, 'r') as f:
            label = f.readline()
            while len(label) > 1:
                label = label.strip()
                p1 = list(map(int, f.readline().split()))
                v1 = np.array(list(map(int, f.readline().split())))
                plt.plot(p1, v1 / milsec, label=label)
                label = f.readline()

    plt.legend()
    plt.title(title)
    plt.xlabel(r"$\log_2 n$")
    plt.ylabel("ms")
    plt.grid()
    plt.yscale('log')
    plt.savefig(name, dpi=300)


def get_args():
    parser = argparse.ArgumentParser(description="Run experiments. Note that not all of the algorithms work properly on signal of size greater than 2^22")
    parser.add_argument("--no_graph", action='store_true', help='whether to generate graph')
    parser.add_argument("--sparsity", action='store', default=32, type=int, help='signal sparsity')
    parser.add_argument("--samples", action='store', default=10, type=int, help='number of samples to average over')
    parser.add_argument("--dimensions", action='store', default=3, type=int, help='number of dimensions')
    parser.add_argument("--first_logn", action='store', default=3, type=int, help='logarithm of the first value of n to test')
    parser.add_argument("--last_logn", action='store', default=7, type=int, help='logarithm of the last value of n to test')
    parser.add_argument("--use_comb_filter", action='store', default=1, type=int, help='whether to use comb filter')
    parser.add_argument("--signal_type", choices=['random', 'comb', 'combined'], default='random', help='choose signal type to run experiments on')
    parser.add_argument("--zero_test_coef", action='store', default=1, type=float, help='ZeroTest coefficient')
    parser.add_argument("--path", action='store', default="./build/measure_run", type=str, help='location of measure_run executable')
    parser.add_argument("--recursive_algorithm_count", action='store', default=4, help="number of recursive algorithms to test")
    parser.add_argument("--output", action='store', default="graph.png", type=str, help='output graph name')
    parser.add_argument("--extra_data", action='store', default='', help='extra data to plot')

    args = parser.parse_args()
    return args
