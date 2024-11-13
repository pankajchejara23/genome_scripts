#!/bin/sh
# This script automates the process of downloading dataset from NCBI along with metadata infromation.
## Usage: ncbi_downloader

## author: pankaj chejara
## email: pankaj.chejara@metrosert.ee

echo "Enter the project identifier (hint: it will be starting from PR...)"
read project

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
esearch -db sra -query $project  | efetch -format runinfo > ${sra}_metadata.txt


# fetch accessions numbers from metada file into a text file
cat ${sra}_metadata.txt | cut -f 1 -d ',' |grep "SRR" > ${sra}_sra_id.txt


# change directory to raw_data
cd raw_data
echo "  changing to raw_data directory"

 # download sra file, fastq file and zip fastq file into gz
cat ../${sra}_sra_id.txt | while read sra_id; do prefetch $sra_id; fasterq-dump $sra_id; gzip ${sra_id}*.fastq;done

echo "Sequence data has been downloaded successfully"
