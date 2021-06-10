import experiment_utils
import subprocess
import sys
import io
import argparse

parser = argparse.ArgumentParser(description="Plot experiment times")
parser.add_argument("--title", action='store', type=str, help='graph title')
parser.add_argument("--output", action='store', type=str, help='output graph name')

args = parser.parse_args()

algs = []

plot_names = {
    'fftw': 'FFTW',
    'recursive 1' : r'disFT, $n^{3}$',
    'recursive 2' : r'disFT, $n^{2\frac{1}{2}}$',
    'recursive 3' : r'disFT, $n^{2\frac{1}{3}}$',
    'recursive 4' : r'disFT, $n^{2\frac{1}{4}}$',
    'random_phase' : r'disFT, random phase'
}


p = list(map(int, input().split()))
results = dict()
while True:
    try:
        name = input().strip()
        assert(name not in algs)
        algs.append(name)
        values = list(map(int, input().split()))
        results[name] = values
    except EOFError:
        break

experiment_utils.make_graph(args.output, args.title, p, algs, plot_names, results, time_bounds=(0.1, 600))
