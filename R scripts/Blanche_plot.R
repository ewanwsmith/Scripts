#!/usr/bin/env Rscript

# Function to install and load packages
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, repos = "http://cran.us.r-project.org")
    library(package, character.only = TRUE)
  }
}

# Install and load required packages
install_and_load("optparse")
install_and_load("ggplot2")
install_and_load("viridis")

# Define command line arguments
option_list <- list(
  make_option(c("--in"), type = "character", default = NULL,
              help = "Path to input .dat file", metavar = "character"),
  make_option(c("--out"), type = "character", default = NULL,
              help = "Path to output .jpeg file", metavar = "character"),
  make_option(c("--names"), type = "character", default = NULL,
              help = "Path to sample names .txt file", metavar = "character")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if no arguments are provided and display help message
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  quit(status = 0)
}

# Check if input file is provided
input_file <- opt$`in`
if (is.null(input_file)) {
  stop("Error: --in argument is required.")
}

# Read the .dat file
data <- read.table(input_file, header = FALSE, col.names = c("x", "y"))

# Check if sample names file is provided
names_file <- opt$names
if (!is.null(names_file)) {
  sample_names <- readLines(names_file)
  if (length(sample_names) != nrow(data)) {
    stop("Error: The number of sample names does 
    not match the number of rows in the .dat file.")
  }
  data$sample <- sample_names
}

# Create a dot plot
plot <- ggplot(data, aes(x = x, y = y, color = sample)) +
  geom_point() +
  labs(title = "Blanche output",
       x = "Dimension 1", y = "Dimension 2", color = "") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE)

# Define output file path
output_file <- if (!is.null(opt$out)) {
  opt$out
} else {
  paste0(gsub("\\.dat$", "", input_file), ".jpeg")
}

# Save the plot as a .jpeg file
ggsave(output_file, plot = plot, device = "jpeg")

# Inform the user
cat("Plot saved as", output_file, "\n")