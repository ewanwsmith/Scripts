#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <path-to-folder>"
  exit 1
fi

# Get the path to the folder from the argument
folder_path=$1

# Extract the folder name for the plot title
folder_name=$(basename "$folder_path")

# Output CSV file
output_file="$folder_path/extracted_bic.csv"

# Write the header to the CSV file (file is now the first column)
echo "file,n_haps,BIC" > "$output_file"

# Loop through all .out files in the specified directory
for file in "$folder_path"/*.out; do
  # Check if the file exists
  if [ -f "$file" ]; then
    # Read the first line of the file
    first_line=$(head -n 1 "$file")
    # Extract the number after "Final BIC" using awk
    bic_number=$(echo "$first_line" | awk '{for (i=1; i<=NF; i++) if ($i=="BIC") print $(i+1)}')
    # Get the base name of the file (without path)
    file_name=$(basename "$file")
    # Extract the run number from the filename
    n_haps=$(echo "$file_name" | awk -F'_' '{print $2}')
    # Check if the number was found and is numeric
    if [[ -n "$bic_number" && "$n_haps" =~ ^[0-9]+$ ]]; then
      # Write the file name, n_haps number, and BIC number to the CSV file
      echo "$file_name,$n_haps,$bic_number" >> "$output_file"
    fi
  fi
done

echo "Extraction complete. Results saved in $output_file."

# Path for bic_plot.R
r_script_path="$folder_path/bic_plot.R"

# Always create or overwrite bic_plot.R
echo "Creating or overwriting bic_plot.R."

# Write the R script with dynamic x-axis label and title
cat <<EOL > "$r_script_path"
library(csv)
library(ggplot2)

# Read the CSV file
bic <- read.csv("$output_file")

# Print the bic object to check its contents
print(bic)

# Remove rows where BIC is NA
bic <- bic[!is.na(bic$BIC), ]

# Ensure that n_haps is treated as numeric and sort by n_haps
bic\$n_haps <- as.numeric(as.character(bic\$n_haps))
bic <- bic[order(bic\$n_haps), ]

# Create the plot with points only (no lines) and dynamic axis and title
bic_plot <- ggplot(bic, aes(x = n_haps, y = BIC)) +
  geom_point() +
  labs(title="$folder_name", x="# Haplotypes", y="BIC") +
  theme_minimal()

# Define the path for saving the plot
output_file <- "$folder_path/bic_plot.jpeg"

# Save the plot as a .jpeg file
ggsave(output_file, plot = bic_plot, width = 8, height = 6)

# No need to plot the graph in the R script
EOL

echo "bic_plot.R created or overwritten."

# Run the R script
echo "Running bic_plot.R"
Rscript "$r_script_path" "$output_file"

# Find the row with the lowest BIC and print the file name
lowest_bic_file=$(awk -F',' 'NR > 1 && $3 != "N/A" { if (min == "" || $3 < min) { min=$3; file=$1 } } END { print file }' "$output_file")

# Print the result after everything is complete
echo "File with the lowest BIC: $lowest_bic_file"