#!/bin/bash

input_dir="../FFMS_all_instances"
threshold=0.8
time=60
files=("200-600" "200-800")
output_file="60s_0.8"

# Crea un archivo temporal para guardar los trabajos
job_file=$(mktemp)

# Genera los comandos para GNU Parallel
for sufix in "${files[@]}"; do
    # Limpia o crea el archivo de salida correspondiente
    > "${output_file}_${sufix}.txt"
    
    # Encuentra todos los archivos que coincidan con el sufijo
    for file in "$input_dir"/"$sufix"-*.txt; do
        # Genera el comando para cada archivo
        echo "./a.out -i \"$file\" -c $time -t $threshold -n 190 -r 0.01 -e 0.1 -x 0.3 --tn | awk '{print \$1}' >> \"${output_file}_${sufix}.txt\"" >> "$job_file"
    done
done

# Ejecuta los trabajos en paralelo utilizando 8 núcleos
cat "$job_file" | parallel -j 6

# Limpia el archivo temporal
rm "$job_file"