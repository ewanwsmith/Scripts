#!/bin/bash

# Set the directory to search (current directory by default)
search_dir="${1:-.}"

# Use find to get all files, sort them by size, and print the largest one
largest_file=$(find "$search_dir" -type f -exec du -h {} + | sort -rh | head -n 1)

# Echo the largest file's path and size
echo "Largest file: $largest_file"