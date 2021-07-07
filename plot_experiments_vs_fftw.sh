#!/bin/bash

python3 run_and_plot_experiments.py --signal_type=twocomb --dimensions=3 --first_logn=3 --last_logn=8 --sparsity=32 --use_comb_filter=0 --use_projection_recovery=0

python3 run_and_plot_experiments.py --signal_type=comb --dimensions=3 --first_logn=3 --last_logn=8 --sparsity=32 --use_comb_filter=0 --use_projection_recovery=0

python3 run_and_plot_experiments.py --signal_type=random --dimensions=3 --first_logn=3 --last_logn=8 --sparsity=32 --use_comb_filter=0 --use_projection_recovery=1

python3 run_and_plot_experiments.py --signal_type=combined --dimensions=3 --first_logn=3 --last_logn=8 --sparsity=32 --use_comb_filter=0 --use_projection_recovery=1

