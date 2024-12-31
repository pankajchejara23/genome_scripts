#!/bin/sh
# This script automates the process of downloading dataset from NCBI along with metadata infromation.
## Usage: ./ncbi_downloader.sh

## author: pankaj chejara
## email: pankaj.chejara@metrosert.ee

# Function to display usage
usage() {
  echo "Usage: $0 -p <project_identifier>"
  echo "  -p    Project identifier (e.g., PRJXXXX)"
  exit 1
}

# Parse command-line options using getopts
while getopts ":p:" opt; do
  case $opt in
    p)
      project=$OPTARG
      ;;
    *)
      usage
      ;;
  esac
done

# Check if the project identifier was provided
if [ -z "$project" ]; then
  echo "Error: Project identifier is required."
  usage
fi

echo "Preparing for the download"

# check if project directory exists
if [ -d "$project" ]; then
  echo "  $project already exists. cd to $project"
else
  echo "  creating a new directory $project"
  mkdir $project
fi

# change to project directory
cd $project
echo "  changing to $project directory"


# check if raw_data directory exists
if [ -d "raw_data" ]; then
  echo "  raw_data already exists. cd to $raw_data"
else
  echo "  creating a new directory raw_data"
  mkdir raw_data
fi

# run search query and store metadata
esearch -db sra -query $project  | efetch -format runinfo > ${project}_metadata.txt


# fetch accessions numbers from metada file into a text file
cat ${project}_metadata.txt | cut -f 1 -d ',' |grep "SRR" > ${project}_sra_id.txt


# change directory to raw_data
cd raw_data
echo "  changing to raw_data directory"

 # download sra file, fastq file and zip fastq file into gz
cat ../${project}_sra_id.txt | while read sra_id; do prefetch $sra_id; fasterq-dump $sra_id; gzip ${sra_id}*.fastq;done

echo "Sequence data has been downloaded successfully"
