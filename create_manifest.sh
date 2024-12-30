#!/usr/bin/bash
# author: Pankaj chejara
# Script to create manifest file for qiime2

samples=()

for filename in ./raw_files/*.fastq.gz; do
    base=$(basename "$filename" .fastq.gz)
    nf=$(echo $base | sed -e 's/.......$//');
    if ! [[ ${samples[@]} =~ $nf ]]
    then
      samples+=("$nf");
    fi
  done

echo "sample-id    forward-absolute-filepath    reverse-absolute-filepath" > manifest.csv

for sample in "${samples[@]}"; do
  forward="$PWD/trimm_outputs/${sample}_1P.fastq.gz"
  backward="$PWD/trimm_outputs/${sample}_2P.fastq.gz"
  sampleid=$(echo $sample | sed -e 's/_.*//')
  echo "$sampleid    $forward    $backward" >> manifest.csv
  #echo "$sampleid,$backward,backward" >> manifest.csv
done;
