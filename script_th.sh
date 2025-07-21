#!/bin/bash

input_dir="FFMS_all_instances"
thresholds=("0.75" "0.775" "0.8")
times=5
files=("100-300" "100-600")
output_file="resultados_th_secuencial"

for sufix in "${files[@]}"; do
    > "${output_file}_${sufix}.txt"

    for threshold in "${thresholds[@]}"; do
        for file in "$input_dir"/${sufix}-*00[1-9].txt; do
            for run in $(seq 1 $times); do
                ./a.out -i "$file" -th "$threshold" -T 20 -alpha 0.9 -num_threads 128 -threads_per_block 8 -iterations 10 \
                | awk -v th="$threshold" '{print th, $4, $7}'  >> "${output_file}_${sufix}.txt"
            done
        done
    done
done
