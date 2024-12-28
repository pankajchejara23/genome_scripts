#!/usr/bin/env bash
# author: Pankaj chejara
# Script to iterate over paired read sequences to apply trimmomatic

samples=()

for filename in ./raw_data/*.fastq.gz; do
    base=$(basename "$filename" .fastq.gz)
    nf=$(echo $base | sed -e 's/.......$//');
    if ! [[ ${samples[@]} =~ $nf ]]
    then
      samples+=("$nf");
    fi
  done

for sample in "${samples[@]}"; do
  forward="./raw_data/${sample}_R1_001.fastq.gz"
  backward="./raw_data/${sample}_R2_001.fastq.gz"
  baseout="./trimm_outputs/${sample}.fastq.gz"
  trimmomatic PE -threads 20 $forward $backward -baseout $baseout CROP:200"
done;
