#!/bin/bash

names=(
#"rand_1" "rand_2" "rand_3" "cubes_2"
"cubes_3" "rand_3_sp64" "rand_3_sp16" "rand_3_sp128" "rand_3_sp256" "cubes_3_sp64" "cubes_3_sp16" "cubes_3_sp128" "cubes_3_sp256" "cubes_4" "random_4")
titles=(
#"dim:1 sparsity:32 type:overtones"
#"dim:2 sparsity:32 type:overtones"
#"dim:3 sparsity:32 type:overtones"
#"dim:2 sparsity:32 type:cubes"
"dim:3 sparsity:32 type:cubes"
"dim:3 sparsity:64 type:overtones"
"dim:3 sparsity:16 type:overtones"
"dim:3 sparsity:128 type:overtones"
"dim:3 sparsity:256 type:overtones"
"dim:3 sparsity:64 type:cubes"
"dim:3 sparsity:16 type:cubes"
"dim:3 sparsity:128 type:cubes"
"dim:3 sparsity:256 type:cubes"
"dim:4 sparsity:32 type:cubes"
"dim:4 sparsity:32 type:random"
)

#names=("rand_3_npr" "rand_3_sp16_npr" "rand_3_sp64_npr" "cubes_3_npr" "cubes_3_sp16_npr" "cubes_3_sp64_npr")
#titles=(
#"dim:3 sparsity:32 type:random PR:disabled"
#"dim:3 sparsity:16 type:random PR:disabled"
#"dim:3 sparsity:64 type:random PR:disabled"
#"dim:3 sparsity:32 type:cubes PR:disabled"
#"dim:3 sparsity:16 type:cubes PR:disabled"
#"dim:3 sparsity:64 type:cubes PR:disabled"
#)