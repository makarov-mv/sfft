import experiment_utils
import subprocess
import sys
import io

args = experiment_utils.get_args()

algs = [
    'fftw',
    'random_phase'
]

for i in range(1, args.recursive_algorithm_count + 1):
    algs.append('recursive ' + str(i))

alg_settings = {
    'fftw' : {'use_comb' : args.use_comb_filter, 'zero_test_coef' : args.zero_test_coef},
    'recursive 1' : {'use_comb' : args.use_comb_filter, 'zero_test_coef' : args.zero_test_coef},
    'recursive 2' : {'use_comb' : args.use_comb_filter, 'zero_test_coef' : args.zero_test_coef},
    'recursive 3' : {'use_comb' : args.use_comb_filter, 'zero_test_coef' : args.zero_test_coef},
    'recursive 4' : {'use_comb' : args.use_comb_filter, 'zero_test_coef' : args.zero_test_coef},
    'random_phase' : {'use_comb' : args.use_comb_filter, 'assume_random_phase' : 1}
}

plot_names = {
    'fftw': 'FFTW',
    'recursive 1' : r'disFT, $n^{3}$',
    'recursive 2' : r'disFT, $n^{2\frac{1}{2}}$',
    'recursive 3' : r'disFT, $n^{2\frac{1}{3}}$',
    'recursive 4' : r'disFT, $n^{2\frac{1}{4}}$',
    'random_phase' : r'disFT, random phase'
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

p = list(map(int, inputs.readline().split()))
results = dict()
for _ in range(len(algs)):
    name = inputs.readline().strip()
    assert(name in algs)
    values = list(map(int, inputs.readline().split()))
    results[name] = values

print(p)
print(results)
if not args.no_graph:
    title = 'Comparison of algorithms, {} dimensions, {} signal, sparsity={}\n{} comb filter'.format(
        args.dimensions, args.signal_type, args.sparsity, 'with' if args.use_comb_filter else 'without'
    )
    experiment_utils.make_graph(args.output, title, p, algs, plot_names, results, args.extra_data)
