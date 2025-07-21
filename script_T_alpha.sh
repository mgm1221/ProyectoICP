#!/bin/bash

input_dir="FFMS_all_instances"
Ts=("1000" "100" "20" "5" "2" "1")
alphas=("0.95" "0.9" "0.8")
times=10
files=("100-600")
output_file="resultados_T_alpha"

for sufix in "${files[@]}"; do
    > "${output_file}_${sufix}.txt"
    for file in $input_dir/$sufix-*00[1-5].txt; do
        for T in "${Ts[@]}"; do
            for alpha in "${alphas[@]}"; do
                for run in $(seq 1 $times); do
                    ./a.out -i "$file" -th 0.775 -T $T -alpha $alpha -num_threads 1 -threads_per_block 1 \
                    | awk -v T="$T" -v alpha="$alpha" '{print T, alpha, $4, $7}' >> "${output_file}_${sufix}.txt"
                done
            done
        done
    done
done

