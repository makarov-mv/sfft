#!/bin/bash

./build/generate_graphs -N -R 100 -W > rand_data.txt
./build/generate_graphs -N -R 100 -C > comb_data.txt
