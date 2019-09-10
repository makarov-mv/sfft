#!/bin/bash

#python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=32 --use_comb_filter=1 --zero_test_coef=0.03 --output=random_with_comb.png
python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=27 --use_comb_filter=1 --extra_data=./measures.txt --output=random_with_comb.png
python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=27 --use_comb_filter=0 --output=random_without_comb.png

python3 run_experiment.py --path=./build/measure_run --signal_type=comb --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=32 --use_comb_filter=1 --output=comb_with_comb.png
python3 run_experiment.py --path=./build/measure_run --signal_type=comb --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=32 --use_comb_filter=0 --output=comb_without_comb.png

python3 run_experiment.py --path=./build/measure_run --signal_type=combined --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=32 --use_comb_filter=1 --output=combined_with_comb.png
python3 run_experiment.py --path=./build/measure_run --signal_type=combined --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=32 --use_comb_filter=0 --output=combined_without_comb.png

python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=1 --first_logn=6 --last_logn=21 --sparsity=32 --use_comb_filter=1 --output=1d_random.png
python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=2 --first_logn=4 --last_logn=10 --sparsity=32 --use_comb_filter=1 --output=2d_random.png
python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=3 --first_logn=3 --last_logn=7 --sparsity=32 --use_comb_filter=1 --output=3d_random.png
python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=4 --first_logn=2 --last_logn=5 --sparsity=32 --use_comb_filter=1 --output=4d_random.png
python3 run_experiment.py --path=./build/measure_run --signal_type=random --dimensions=5 --first_logn=2 --last_logn=4 --sparsity=32 --use_comb_filter=1 --output=5d_random.png