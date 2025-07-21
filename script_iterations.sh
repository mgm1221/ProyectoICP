#!/bin/bash

input_dir="FFMS_all_instances"
times=10
files=("100-300" "100-600")
iterations_list=(20 10 5 1)

output_file="resultados_iter"

for sufix in "${files[@]}"; do
    > "${output_file}_${sufix}.txt"

    for file in $input_dir/${sufix}-*00[1-5].txt; do
        for iterations in "${iterations_list[@]}"; do
                for run in $(seq 1 $times); do
                    ./a.out -i "$file" -th 0.775 -T 5000 -alpha 0.95 -num_threads 10 -threads_per_block 10  -iterations "$iterations"\
                    | awk -v th="0.775" -v it="$iterations" \
                        '{print th, it, $4, $7}' >> "${output_file}_${sufix}.txt"
               	done
        done
    done
done

