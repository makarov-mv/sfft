#!/bin/bash
# run run_sfft.sh first

python3 run_and_plot_experiments.py --signal_type=comb --dimensions=1 --first_logn=14 --last_logn=24 --sparsity=32 --use_comb_filter=0 --use_projection_recovery=0 --extra_data=./sfft-v1v2/comb_data.txt

python3 run_and_plot_experiments.py --signal_type=random --dimensions=1 --first_logn=14 --last_logn=24 --sparsity=32 --use_comb_filter=0 --use_projection_recovery=0 --extra_data=./sfft-v1v2/rand_data.txt
