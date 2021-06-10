#!/bin/bash

source experiment_lists.sh
MEASURE_PATH=./build/measure_run

python3 gen_uneven_info.py --dim=1 --type=random > ./test_results/rand_1_data.txt
python3 gen_uneven_info.py --dim=2 --type=random > ./test_results/rand_2_data.txt
python3 gen_uneven_info.py --dim=3 --type=random > ./test_results/rand_3_data.txt
python3 gen_uneven_info.py --dim=2 --type=cubes > ./test_results/cubes_2_data.txt
python3 gen_uneven_info.py --dim=3 --type=cubes > ./test_results/cubes_3_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=16 --type=random > ./test_results/rand_3_sp16_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=64 --type=random > ./test_results/rand_3_sp64_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=128 --type=random > ./test_results/rand_3_sp128_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=256 --type=random > ./test_results/rand_3_sp256_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=16 --type=cubes > ./test_results/cubes_3_sp16_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=64 --type=cubes > ./test_results/cubes_3_sp64_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=128 --type=cubes > ./test_results/rand_3_sp128_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=256 --type=cubes > ./test_results/rand_3_sp256_data.txt
python3 gen_uneven_info.py --dim=4 --type=cubes > ./test_results/cubes_4_data.txt
python3 gen_uneven_info.py --dim=4 --type=random > ./test_results/random_4_data.txt

python3 gen_uneven_info.py --dim=3 --type=random --use_projection_recovery=0  > ./test_results/rand_3_npr_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=16 --type=random --use_projection_recovery=0  > ./test_results/rand_3_sp16_npr_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=64 --type=random --use_projection_recovery=0  > ./test_results/rand_3_sp64_npr_data.txt
python3 gen_uneven_info.py --dim=3 --type=cubes --use_projection_recovery=0  > ./test_results/cubes_3_npr_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=16 --type=cubes --use_projection_recovery=0  > ./test_results/cubes_3_sp16_npr_data.txt
python3 gen_uneven_info.py --dim=3 --sparsity=64 --type=cubes --use_projection_recovery=0  > ./test_results/cubes_3_sp64_npr_data.txt


for i in ${!names[@]}; do
  date +"%T"
  echo "Running ${names[$i]}"
  $MEASURE_PATH < ./test_results/${names[$i]}_data.txt > ./test_results/${names[$i]}_res.txt
done