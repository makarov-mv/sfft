#!/bin/bash

python3 run_experiments_vs_sparsity.py --signal_type=twocomb --logn=7 --dimensions=3 --first_logk=5 --last_logk=8 --use_comb_filter=0 --use_projection_recovery=0

python3 run_experiments_vs_sparsity.py --signal_type=comb --logn=7 --dimensions=3 --first_logk=5 --last_logk=8 --use_comb_filter=0 --use_projection_recovery=0

python3 run_experiments_vs_sparsity.py --signal_type=combined --logn=7 --dimensions=3 --first_logk=5 --last_logk=8 --use_comb_filter=0 --use_projection_recovery=0
