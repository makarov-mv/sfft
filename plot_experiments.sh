#!/bin/bash

source experiment_lists.sh

for i in ${!names[@]}; do
  echo "PLotting ${names[$i]}"
  python3 plot_time.py --title="${titles[$i]}" --output=test_results/${names[$i]}.png < test_results/${names[$i]}_res.txt
done