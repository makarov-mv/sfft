import subprocess
import sys
import io
import argparse
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.font_manager import FontProperties

font = FontProperties()
font.set_family('serif')
font.set_name('Times New Roman')
font.set_size('x-large')
font.set_weight('book')


linestyles = ['-', '--', '-.', ':']
color_ = ['mediumseagreen', 'dodgerblue', 'orangered', 'gold']
ecolor_ = ['limegreen', 'cornflowerblue', 'tomato', 'yellow']

def make_graph(name, title, p, algs, algs_names, data, samples, extra_data_path=''):
    ax = plt.axes()
    ax.set_xticks(p)
    milsec = 1e6
    
    runtimes_mean = np.zeros((len(data), len(p), len(data)))
    runtimes_std = np.zeros(np.shape(runtimes_mean))
    for j in range(len(data)):
        for i, alg in enumerate(algs):
            runtimes = np.reshape(np.array(data[j][alg])/milsec, (-1, samples))
            runtimes[runtimes<0] = np.nan
            fail_prob = np.mean(np.isnan(runtimes), axis=1)
            print(alg, fail_prob)
            runtimes_mean[j,:,i] = np.nanmean(runtimes, axis=1)
            runtimes_mean[j,fail_prob>0.1,i] = np.nan
            runtimes_std[j,:,i] = np.nanstd(runtimes, axis=1)
    
    for i, alg in enumerate(algs):
        best_idx = np.nanargmin(runtimes_mean[:,:,i], axis=0)
        print(best_idx)
        best_time = runtimes_mean[best_idx,np.arange(len(p)),i]
        best_std = runtimes_std[best_idx,np.arange(len(p)),i]
        plt.errorbar(2**p, best_time, yerr= best_std,linestyle=linestyles[i], color= color_[i], linewidth=2.2, capsize=6, ecolor=ecolor_[i], label=algs_names[alg])
        
    if extra_data_path:
        with open(extra_data_path, 'r') as f:
            label = f.readline()
            while len(label) > 1:
                label = label.strip()
                p1 = list(map(int, f.readline().split()))
                v1 = np.array(list(map(int, f.readline().split())))
                plt.plot(p1, v1 / milsec, label=label)
                label = f.readline()

    plt.legend(prop=font, loc='best', frameon=True,fancybox=True,framealpha=0.8,edgecolor='k')
    plt.title(title, fontproperties = font)
    plt.xlabel("Signal Size N", fontproperties = font)
    plt.ylabel("Runtime (ms)", fontproperties = font)
    plt.grid()
    plt.yscale('log')
    plt.xscale('log', base=2)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=13)
    plt.savefig(name, dpi=500, bbox_inches='tight')



def run_measure_run(algs, args):
    alg_settings = {
        'fftw' : {'use_comb' : args.use_comb_filter, 'zero_test_coef' : args.zero_test_coef},
        'recursive 1' : {'use_comb' : args.use_comb_filter, 'use_projection_recovery' : args.use_projection_recovery, 'zero_test_coef' : args.zero_test_coef},
        'recursive 2' : {'use_comb' : args.use_comb_filter, 'use_projection_recovery' : args.use_projection_recovery, 'zero_test_coef' : args.zero_test_coef}#,
        #'recursive 3' : {'use_comb' : args.use_comb_filter,  'use_projection_recovery' : args.use_projection_recovery, 'zero_test_coef' : args.zero_test_coef},
    }
        
    proc = subprocess.Popen(args.path, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='ascii')
    output = io.StringIO()
    
    print("samples:", args.samples, file=output)
    
    print(len(algs), file=output)
    print(*algs, file=output)
    
    start = args.first_logn
    finish = args.last_logn
    
    print(finish - start + 1, file=output)
    
    indent = ' ' * 3
    
    for p in range(start, finish + 1):
        print("id:", p, file=output)
        print('info:', args.dimensions, p, file=output)
        print('sparsity:', args.sparsity, file=output)
        print('signal:', args.signal_type, file=output)
        for alg in algs:
            print(alg, file=output)
            print(indent, len(alg_settings[alg]), file=output)
            for prop in alg_settings[alg].items():
                print(indent, *prop, file=output)
    
    output.flush()
    
    outs, errs = proc.communicate(input=output.getvalue())
    
    print(errs)
    
    inputs = io.StringIO(outs)
    
    #print(inputs.read())
    p = np.array(list(map(int, inputs.readline().split())))*args.dimensions
    results = dict()
    for _ in range(len(algs)):
        name = inputs.readline().strip()
        assert(name in algs)
        values = list(map(int, inputs.readline().split()))
        results[name] = values
    
    print(p)
    #print(results)
    return p, results
    



parser = argparse.ArgumentParser(description="Run experiments. Note that not all of the algorithms work properly on signal of size greater than 2^22")
parser.add_argument("--no_graph", action='store_true', help='whether to generate graph')
parser.add_argument("--sparsity", action='store', default=32, type=int, help='signal sparsity')
parser.add_argument("--samples", action='store', default=50, type=int, help='number of samples to average over')
parser.add_argument("--dimensions", action='store', default=2, type=int, help='number of dimensions')
parser.add_argument("--first_logn", action='store', default=5, type=int, help='logarithm of the first value of n to test')
parser.add_argument("--last_logn", action='store', default=12, type=int, help='logarithm of the last value of n to test')
parser.add_argument("--use_comb_filter", action='store', default=0, type=int, help='whether to use comb filter')
parser.add_argument("--use_projection_recovery", action='store', default=0, type=int, help='whether to use projection recovery')
parser.add_argument("--signal_type", choices=['random', 'comb', 'combined', 'twocomb'], default='twocomb', help='choose signal type to run experiments on')
parser.add_argument("--zero_test_coef", action='store', default=0.25, type=float, help='ZeroTest coefficient')
parser.add_argument("--path", action='store', default="./build-xcode/Release/measure_run", type=str, help='location of measure_run executable')
parser.add_argument("--recursive_algorithm_count", action='store', default=2, help="number of recursive algorithms to test")
parser.add_argument("--output", action='store', default="Graph.pdf", type=str, help='output graph name')
parser.add_argument("--extra_data", action='store', default='', help='extra data to plot')

args = parser.parse_args()

algs = [
    'fftw'
]

plot_names = {
    'fftw': 'FFTW',
    'recursive 1' : r'SparseFT, $k^{3}$',
    'recursive 2' : r'SparseFT, $k^{2\frac{1}{2}}$'#,
    #'recursive 3' : r'disFT, $k^{2\frac{1}{3}}$',
}

for i in range(1, args.recursive_algorithm_count + 1):
    algs.append('recursive ' + str(i))

title = '{}-dimensional {} signal, sparsity={}, {} projection'.format(
        args.dimensions, args.signal_type, args.sparsity, 'with' if args.use_projection_recovery else 'without'
    )

if args.signal_type == 'random':
    coeffs = 2**np.linspace(-np.log2(args.sparsity), -2, num=4)
else:
    coeffs = 2**np.linspace(-3.5, -0.5, num=4)
    
results = [0 for i in range(len(coeffs))]
for i in range(len(results)):
    args.zero_test_coef = coeffs[i]
    p, results[i] = run_measure_run(algs, args)

make_graph('{}-dimensional-{}{}'.format(
        args.dimensions, args.signal_type, '-Proj' if args.use_projection_recovery else '') + args.output, title, p, algs, plot_names, results, args.samples, args.extra_data)


