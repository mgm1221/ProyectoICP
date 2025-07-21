#!/bin/bash

input_dir="FFMS_all_instances"
times=5
files=("100-300" "100-600")
num_threads_list=(128 64 32)
threads_per_block_list=(32 16 8)
output_file="resultados_blocks_threads_v2"

for sufix in "${files[@]}"; do
    > "${output_file}_${sufix}.txt"

    for file in $input_dir/${sufix}-*00[1-5].txt; do
        for num_threads in "${num_threads_list[@]}"; do
            for threads_per_block in "${threads_per_block_list[@]}"; do

                for run in $(seq 1 $times); do
                    ./a.out -i "$file" -th 0.775 -T 20 -alpha 0.9 -num_threads "$num_threads" -threads_per_block "$threads_per_block" -iterations 10\
                    | awk -v th="0.775" -v nt="$num_threads" -v tpb="$threads_per_block" -v blk="$blocks" \
                        '{print th, nt, tpb, blk, $4, $7}' >> "${output_file}_${sufix}.txt"
                done
            done
        done
    done
done
