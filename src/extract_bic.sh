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

# Write the header to the CSV file
echo "file,n_haps,BIC" > "$output_file"

# Loop through all .out files in the specified directory
for file in "$folder_path"/*.out; do
  if [ -f "$file" ]; then
    first_line=$(head -n 1 "$file")
    bic_number=$(echo "$first_line" | awk '{for (i=1; i<=NF; i++) if ($i=="BIC") print $(i+1)}')
    file_name=$(basename "$file")
    n_haps=$(echo "$file_name" | awk -F'_' '{print $2}')
    if [[ -n "$bic_number" && "$n_haps" =~ ^[0-9]+$ ]]; then
      echo "$file_name,$n_haps,$bic_number" >> "$output_file"
    fi
  fi
done

echo "Extraction complete. Results saved in $output_file."

# Path for bic_plot.R
r_script_path="$folder_path/bic_plot.R"

echo "Creating or overwriting bic_plot.R."

# Write the R script
cat <<EOL > "$r_script_path"
library(ggplot2)

# Read the CSV file
bic <- read.csv("$output_file")

# Remove rows where BIC is NA
bic <- bic[!is.na(bic$BIC), ]

# Ensure n_haps is numeric and sort by n_haps
bic\$n_haps <- as.numeric(as.character(bic\$n_haps))
bic <- bic[order(bic\$n_haps), ]

# Create the plot
bic_plot <- ggplot(bic, aes(x = n_haps, y = BIC)) +
  geom_point() +
  labs(title="$folder_name", x="# Haplotypes", y="BIC") +
  theme_minimal()

# Save the plot
ggsave("$folder_path/bic_plot.jpeg", plot = bic_plot, width = 8, height = 6)

# Find the row with the lowest BIC
min_bic <- bic[which.min(bic\$BIC), ]

# Print the file with the lowest BIC
cat("File with the lowest BIC:", min_bic\$file, "\n")
EOL

echo "bic_plot.R created or overwritten."

# Run the R script
Rscript "$r_script_path"

echo "R script executed."