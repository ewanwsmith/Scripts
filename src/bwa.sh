#!/bin/bash

# Set the BWA directory
BWA_DIR="/Users/e.smith.5/Documents/PhD/Software/bwa"

# Function to display usage
usage() {
  echo "Usage: $0 -r <reference.fasta> -s <read.fq.gz> | -p <read1.fq> <read2.fq>"
  exit 1
}

# Parse input arguments
while getopts "r:s:p:" opt; do
  case $opt in
    r) REFERENCE=$OPTARG ;;
    s) READ=$OPTARG ; MODE="single" ;;
    p) READ1=$OPTARG ; READ2=${!OPTIND} ; MODE="paired" ;;
    *) usage ;;
  esac
done

# Check that reference is provided
if [ -z "$REFERENCE" ]; then
  usage
fi

# Change directory to the BWA directory
cd "$BWA_DIR" || { echo "Failed to change directory to $BWA_DIR"; exit 1; }

# Index the reference fasta file
./bwa index "$REFERENCE"

# Run alignment based on the mode
if [ "$MODE" == "single" ]; then
  # Extract the directory and filename without extension for single-end read
  OUTPUT_DIR=$(dirname "$READ")
  BASENAME=$(basename "$READ" .fq.gz)
  OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_aln-se.sam.gz"
  
  # Run bwa mem for single-end reads
  ./bwa mem "$REFERENCE" "$READ" | gzip -3 > "$OUTPUT_FILE"
  echo "Single-end alignment completed: $OUTPUT_FILE"

elif [ "$MODE" == "paired" ]; then
  # Extract the directory and filenames without extension for paired-end reads
  OUTPUT_DIR=$(dirname "$READ1")
  BASENAME1=$(basename "$READ1" .fq)
  BASENAME2=$(basename "$READ2" .fq)
  OUTPUT_FILE="$OUTPUT_DIR/${BASENAME1}_${BASENAME2}_aln-pe.sam.gz"
  
  # Run bwa mem for paired-end reads
  ./bwa mem "$REFERENCE" "$READ1" "$READ2" | gzip -3 > "$OUTPUT_FILE"
  echo "Paired-end alignment completed: $OUTPUT_FILE"

else
  usage
fi