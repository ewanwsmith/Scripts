#!/usr/bin/env julia

using Pkg

# Activate the environment located in the current directory (the folder the script is in)
Pkg.activate(@__DIR__)

# Optionally, instantiate the environment to install all required packages
Pkg.instantiate()

# Install required Julia packages if they are not already installed
function install_packages()
    packages = ["ArgParse", "RCall"]
    installed_pkgs = Pkg.dependencies()

    # Check if each package is installed and install if it is not
    for pkg in packages
        if !any(p -> p.name == pkg, values(installed_pkgs))
            Pkg.add(pkg)
        end
    end
end

install_packages()

using ArgParse
using RCall

# Function to install required R packages
function install_R_packages()
    R"""
    required_packages <- c("Biostrings", "stringr")
    installed_packages <- rownames(installed.packages())
    for (pkg in required_packages) {
        if (!pkg %in% installed_packages) {
            install.packages(pkg, dependencies = TRUE)
        }
    }
    """
end

install_R_packages()

# Function to handle the .dat file processing
function pull_frames(dat_file_path::String)
    # Define the paths for the .fa files in Julia before calling R
    nucleotide_fasta_path = replace(dat_file_path, ".dat" => "_nucleotide.fa")
    aa_fasta_path = replace(dat_file_path, ".dat" => "_amino_acid.fa")

    # Pass the file paths to R
    @rput dat_file_path
    @rput nucleotide_fasta_path
    @rput aa_fasta_path

    R"""
    library(Biostrings)
    library(stringr)
    process_file <- function(dat_file_path, fasta_nucleotide_path, fasta_amino_acid_path) {
        lines <- tryCatch({
            readLines(dat_file_path)
        }, error = function(e) {
            cat("Error reading file:", e$message, "\n")
            return(NULL)  # Return NULL to handle error gracefully
        })
        if (is.null(lines)) return()  # Exit if file could not be read

        con_nuc <- file(fasta_nucleotide_path, open = "w")
        con_aa <- file(fasta_amino_acid_path, open = "w")
        
        for (i in seq(1, length(lines), by = 3)) {
            if (i + 2 > length(lines)) {
                break
            }
            header <- strsplit(lines[i], " ")[[1]]
            if (length(header) < 3) {
                next
            }
            start_pos <- header[1]
            end_pos <- header[2]
            annotation <- paste(header[3:length(header)], collapse=" ")
            full_annotation <- paste0('"accession":', start_pos, '..', end_pos, ' | ', annotation)
            dna_seq <- DNAString(lines[i+1])
            aa_seq <- AAString(lines[i+2])
            
            writeLines(paste0(">", full_annotation), con_nuc)
            writeLines(as.character(dna_seq), con_nuc)
            writeLines(paste0(">", full_annotation), con_aa)
            writeLines(as.character(aa_seq), con_aa)
        }
        
        close(con_nuc)
        close(con_aa)
    }
    """
    # Call the R function from Julia with paths
    R"process_file($(dat_file_path), $(nucleotide_fasta_path), $(aa_fasta_path))"
    
    # Julia prints and retains the paths in variables
    println("nucleotide .fa : ", nucleotide_fasta_path)
    println("amino acid .fa : ", aa_fasta_path)

    return nucleotide_fasta_path, aa_fasta_path
end

# Main function to handle command-line arguments
function main()
    # Setup argument parser
    s = ArgParseSettings()
    @add_arg_table s begin
        "--in", "-i"
        help = "Path to the .dat file"
        required = true
    end

    parsed_args = parse_args(s)
    dat_file_path = parsed_args["in"]
    
    # Call the pull_frames function with the provided file path
    nucleotide_fasta_path, aa_fasta_path = pull_frames(dat_file_path)
end

# Execute the main function if the script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end