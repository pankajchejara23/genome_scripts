#!/bin/sh
# Script to create manifest file for QIIME import

## author: pankaj chejara
## email: pankaj.chejara@metrosert.ee

usage() {
  echo "Usage: $0 -i <input_directory> -o <output_directory>"
  echo "  -i    Directory containing paired read sequences"
  echo "  -o    Direcotry to store results"
  exit 1
}

# Parse command-line options using getopts
while getopts ":i:o:" opt; do
  case $opt in
    i)
      input_dir=$OPTARG
      ;;
    o)
      output_dir=$OPTARG
      ;;
    *)
      usage
      ;;
  esac
done

# Check if the input directory was provided
if [ -z "$input_dir" ]; then
  echo "Error: Input directory is required."
  usage
fi

# Array to store unique sample ids
samples=()

# Iterate over each file in raw_data directory; 
# File names are like "Ademona1-2065_S1_L001_R1_001.fastq.gz"; Change as per your needs
for filename in ./$input_dir/*.fastq.gz; do
    # Get the file name without extension
    base=$(basename "$filename" .fastq.gz)

    # Remove last 7 characters from the filename to get sampleid
    nf=$(echo $base | sed -e 's/.......$//');

    # Add the sample id in the array only if it does not exist
    if ! [[ ${samples[@]} =~ $nf ]]
    then
      samples+=("$nf");
    fi
  done

# Add header to the TSV file
echo "sample-id    forward-absolute-filepath    reverse-absolute-filepath" > manifest.tsv

# Iterate over each unique sample id
for sample in "${samples[@]}"; do
  # Build forward read filename (Note: this is the file name produced by trimmomatic)
  forward="$PWD/${output_dir}/${sample}_1P.fastq.gz"

  # Build backward read filename
  backward="$PWD/${output_dir}/${sample}_2P.fastq.gz"

  # Extract sample id (e.g., Ademona1-2065)
  sampleid=$(echo $sample | sed -e 's/_.*//')

  # Adding row to the TSV file
  echo "$sampleid    $forward    $backward" >> manifest.csv

done;
