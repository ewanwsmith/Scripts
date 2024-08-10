#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i input.fasta -o output.fasta -s start -e end [-z y|n]"
    echo
    echo "Arguments:"
    echo "  -i input.fasta    Input FASTA file."
    echo "  -o output.fasta   Output FASTA file."
    echo "  -s start          Start position (1-based indexing by default)."
    echo "  -e end            End position (1-based indexing by default)."
    echo "  -z y|n            Adjust for zero-based indexing. Use 'y' for yes, 'n' for no (default: n)."
    exit 1
}

# Set default value for zero-based indexing adjustment
zero_based="n"

# Parse command line arguments
while getopts ":i:o:s:e:z:" opt; do
    case $opt in
        i) input_file="$OPTARG" ;;
        o) output_file="$OPTARG" ;;
        s) start="$OPTARG" ;;
        e) end="$OPTARG" ;;
        z) zero_based="$OPTARG" ;;
        *) usage ;;
    esac
done

# If no arguments are provided, display usage information
if [[ $# -eq 0 ]]; then
    usage
fi

# Check if all mandatory arguments are provided
if [[ -z "$input_file" || -z "$output_file" || -z "$start" || -z "$end" ]]; then
    usage
fi

# Ensure start and end are integers and start is less than end
if ! [[ "$start" =~ ^[0-9]+$ ]] || ! [[ "$end" =~ ^[0-9]+$ ]] || [[ "$start" -ge "$end" ]]; then
    echo "Error: Start and end must be integers, and start must be less than end."
    exit 1
fi

# Adjust for zero-based indexing if needed
if [[ "$zero_based" == "y" ]]; then
    start=$((start - 1))
    end=$((end - 1))
fi

# Process the input FASTA file
{
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            # Header line, print it to the output file
            echo "$line"
        else
            # Sequence line, trim and print to the output file
            echo "${line:$start:$(($end - $start + 1))}"
        fi
    done
} < "$input_file" > "$output_file"

echo "Processing complete. Trimmed sequences saved to $output_file."