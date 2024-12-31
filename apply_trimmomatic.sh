#!/bin/sh
# Script to iteratively apply trimmomatic on paired read sequences

## author: pankaj chejara
## email: pankaj.chejara@metrosert.ee

# Function to display usage
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

# Check if the output directory was provided
if [ -z "$output_dir" ]; then
  echo "Error: Output directory is required."
  usage
fi

# array to store sample ids
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

# Iterate over each unique sample id
for sample in "${samples[@]}"; do
  # Build forward read file name
  forward="./${input_dir}/${sample}_R1_001.fastq.gz"

  # Build backward read file name
  backward="./${input_dir}/${sample}_R2_001.fastq.gz"

  # Filename suffix to be used for resultant file
  baseout="./${output_dir}/${sample}.fastq.gz"

  # Apply trimmomatic command
  trimmomatic PE -threads 20 $forward $backward -baseout $baseout CROP:200"
done;
