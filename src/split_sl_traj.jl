using Pkg

# Add necessary packages
Pkg.add("DataFrames")
Pkg.add("BioSequences")
Pkg.add("CSV")
Pkg.add("CategoricalArrays")
Pkg.add("ArgParse")

using DataFrames
using BioSequences
using CSV
using CategoricalArrays
using ArgParse
using RCall

# Ensure the required R libraries are installed
R"""
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("Biostrings")
}
"""

# Function to read a .fasta file and return a DataFrame
function read_fasta_to_dataframe(fasta_path::String)
    # Initialize an empty DataFrame
    df = DataFrame(Accession = String[], Start = Int[], End = Int[], Protein = String[], Sequence = String[])
    
    # Open the fasta file and parse the sequences
    open(fasta_path, "r") do file
        accession = ""
        Start = 0
        End = 0
        protein = ""
        seq = ""
        for line in eachline(file)
            if startswith(line, ">")
                if !isempty(accession)
                    push!(df, (accession, Start, End, protein, seq))
                end
                # Parse the title line
                parts = strip(line[2:end])
                acc_part, rest = split(parts, ":")
                range_part, prot_part = split(rest, "|")
                accession = acc_part
                Start, End = parse.(Int, split(range_part, ".."))
                protein = replace(strip(prot_part), r"[()]" => "")
                seq = ""
            else
                seq *= strip(line)
            end
        end
        if !isempty(accession)
            push!(df, (accession, Start, End, protein, seq))
        end
    end
    
    return df
end

# Function to read .out file and split it based on protein information
function split_out_file(out_file::String, protein_info::DataFrame, out_dir::String)
    println("Reading the .out file...")
    # Read .out file line by line
    lines = readlines(out_file)
    
    # Initialize DataFrame to store locus and remainder
    out_data = DataFrame(locus = Int[], remainder = String[])
    
    for line in lines
        # Split the line into components
        parts = split(line)
        # Extract the first integer as locus
        locus = parse(Int, parts[1])
        # Join the rest as the remainder
        remainder = join(parts[2:end], " ")
        # Push to DataFrame
        push!(out_data, (locus, remainder))
    end

    println(".out file read successfully. Processing proteins...")
    println("Data from .out file:")
    println(out_data)
    
    for row in eachrow(protein_info)
        protein_name = replace(replace(row[:Protein], " " => "_"), r"[()]" => "")
        start_pos = row[:Start]
        end_pos = row[:End]

        # Filter rows where locus is between start and stop positions
        filtered_data = filter(row -> row[:locus] >= start_pos && row[:locus] <= end_pos, out_data)

        if nrow(filtered_data) > 0
            # Create new .out file for the protein in the specified output directory
            new_out_file = joinpath(out_dir, replace(basename(out_file), ".out" => "_$protein_name.out"))
            open(new_out_file, "w") do io
                for row in eachrow(filtered_data)
                    println(io, "$(row.locus) $(row.remainder)")
                end
            end
            println("Written file: $new_out_file")
        else
            println("No data for protein $protein_name in range $start_pos - $end_pos")
        end
    end
    println("All proteins processed successfully.")
end

# Define the function to read the .out file
function read_out_file(filepath::String)
    # Initialize an empty DataFrame
    df = DataFrame(pos = Int[], original_base = String[], variant_base = String[], rest = String[])
    
    # Open the file and process each line
    open(filepath, "r") do file
        for line in eachline(file)
            # Remove any leading/trailing whitespace
            stripped_line = strip(line)
            
            # Continue only if the line is not empty
            if !isempty(stripped_line)
                # Split the line into parts, typically space-delimited
                parts = split(stripped_line)
                
                # Ensure there are enough parts to avoid index errors
                if length(parts) < 3
                    continue  # Skip this line if it doesn't have enough parts
                end

                # Extract the position, original base, and variant base
                pos = parse(Int, parts[1])
                original_base = parts[2]  # Should be a single character string
                variant_base = parts[3]   # Should be a single character string
                
                # Concatenate the remaining parts to form the 'rest' string
                rest = join(parts[4:end], " ")

                # Append to the DataFrame
                push!(df, (pos, original_base, variant_base, rest))
            end
        end
    end
    
    return df
end

# match proteins by position 
function call_proteins(variants_df::DataFrame, fasta_df::DataFrame)
    # Initialize an empty DataFrame for results, pre-defining the columns by combining both DataFrame's columns
    cols = names(variants_df)
    append!(cols, names(fasta_df))
    result_df = DataFrame(; [Symbol(col) => Any[] for col in cols]...)

    # Iterate over each row in variants_df
    for variant in eachrow(variants_df)
        pos = variant[:pos]  # Accessing the 'pos' column in variants_df

        # Filter fasta_df to find rows where 'pos' is between 'Start' and 'End'
        matches = filter(row -> row[:Start] <= pos <= row[:End], fasta_df)

        # For each match, concatenate the row from variants_df with the row from fasta_df
        for match in eachrow(matches)
            # Create a new row by extracting the data from both variant and match
            new_row = [variant[coln] for coln in names(variants_df)]
            append!(new_row, [match[coln] for coln in names(fasta_df)])
            push!(result_df, new_row)
        end
    end

    return result_df
end

# retype columns in the dataframe
function convert_variant_columns!(variants::DataFrame)
    # Convert 'pos' to Int64 if it's not already an integer
    variants.pos = variants.pos isa AbstractArray{Int} ? variants.pos : parse.(Int, string.(variants.pos))
    
    # Convert 'original_base', 'variant_base', and 'Sequence' to String
    variants.original_base = string.(variants.original_base)
    variants.variant_base = string.(variants.variant_base)
    variants.Sequence = string.(variants.Sequence)
    
    # Convert 'Protein' to categorical
    variants.Protein = categorical(variants.Protein)
    
    # Convert 'Start' and 'End' to Int64 if they're not already integers
    variants.Start = variants.Start isa AbstractArray{Int} ? variants.Start : parse.(Int, string.(variants.Start))
    variants.End = variants.End isa AbstractArray{Int} ? variants.End : parse.(Int, string.(variants.End))
end

# create variant sequence
function substitute_variants(dataframe::DataFrame)
    # Create a new column for Base_Position
    dataframe.adj_pos = dataframe.pos .- dataframe.Start

    # Create a new column for variant_sequence
    dataframe.variant_sequence = Vector{String}(undef, nrow(dataframe))

    # Iterate over each row in the DataFrame
    for i in 1:nrow(dataframe)
        # Extract the original sequence
        sequence = string(dataframe[i, :Sequence])

        # Extract the variant position
        adjusted_variant_position = dataframe[i, :adj_pos]

        # Check if the variant position is within the sequence length
        if 1 <= adjusted_variant_position <= length(sequence)
            # Extract the variant base
            variant_base = string(dataframe[i, :variant_base])

            # Check if the character at Base_Position is equal to Original_Base
            original_base = string(sequence[adjusted_variant_position])
            if original_base != dataframe[i, :original_base]
                println("Warning: Original_Base in row $i does not match base at position $adjusted_variant_position in the Sequence.")
            end

            # Create the new sequence with the substitution
            new_sequence = string(sequence[1:adjusted_variant_position-1], variant_base, sequence[adjusted_variant_position+1:end])

            # Update the variant_sequence column
            dataframe[i, :variant_sequence] = new_sequence
        end
    end

    return dataframe
end

# Function to split a DNA sequence into codons
function split_into_codons(sequence::String)
    return [sequence[i:i+2] for i in 1:3:length(sequence)-2]
end

# Function to find original codons
function find_original_codons(df::DataFrame)
    codons = []
    codon_positions = []
    for row in eachindex(df[:, 1])
        sequence = string(df[row, :Sequence])
        base_position = df[row, :adj_pos]

        if base_position < 1 || base_position > length(sequence)
            throw(ArgumentError("Invalid Base_Position for row $row"))
        end

        codon_index = (base_position - 1) ÷ 3 + 1
        codon_position = (base_position - 1) % 3 + 1
        push!(codons, split_into_codons(sequence)[codon_index])
        push!(codon_positions, codon_position)
    end

    df[!, :original_codon] = codons
    df[!, :codon_position] = codon_positions
    df.codon_position = categorical(df.codon_position)
    
    return df
end

# Function to find variant codons
function find_variant_codons(df::DataFrame)
    codons = []
    for row in eachindex(df[:, 1])
        sequence = string(df[row, :variant_sequence])
        base_position = df[row, :adj_pos]

        if base_position < 1 || base_position > length(sequence)
            throw(ArgumentError("Invalid Base_Position for row $row"))
        end

        codon_index = (base_position - 1) ÷ 3 + 1
        push!(codons, split_into_codons(sequence)[codon_index])
    end

    df[!, :variant_codon] = codons
    
    return df
end

# Function to map codons to amino acids
function codon_to_aa(codon::String)
    codon_dict = Dict("TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
                      "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
                      "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
                      "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
                      "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
                      "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
                      "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
                      "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
                      "TAT" => "Y", "TAC" => "Y", "TAA" => "*", "TAG" => "*",
                      "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
                      "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
                      "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
                      "TGT" => "C", "TGC" => "C", "TGA" => "*", "TGG" => "W",
                      "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
                      "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
                      "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G")

    return get(codon_dict, codon, "Unknown")
end

# Function to translate original & variant codons and determine synonymity
function translate_codons(df::DataFrame)
    df.original_aa = map(codon_to_aa, df.original_codon)
    df.variant_aa = map(codon_to_aa, df.variant_codon)
    df.is_synonymous = ifelse.(df.original_aa .== df.variant_aa, "Yes", "No")

    return df
end

function filter_outfile(df::DataFrame, infile_path::String)
    # Read the .out file into a DataFrame
    out_df = CSV.read(infile_path, DataFrame, delim=' ', header=false)
    
    # Get the "pos" column from the input DataFrame
    pos_column = df.pos
    
    # Filter the out_df to only include rows whose first integer is in the "pos" column
    filtered_out_df = out_df[in.(out_df.Column1, Ref(pos_column)), :]
    
    # Construct the output file path
    outfile_path = infile_path[1:end-4] * "_ns.out"
    
    # Write the filtered DataFrame to the new .out file
    CSV.write(outfile_path, filtered_out_df, delim=' ', writeheader=false)
    
    return outfile_path  # Return the path to the filtered output file
end

# Function to extract data from a .fasta file using Biostrings in R
function pull_fasta(fasta_path::String)
    println("Reading and processing the fasta file...")
    # Read the fasta file using Biostrings in R
    R"""
    library(Biostrings)
    library(stringr)
    dna_seqs <- readDNAStringSet($fasta_path)
    headers <- names(dna_seqs)
    
    # Parse headers to extract Accession, Start, End, and Protein
    accession <- str_extract(headers, "^[^:]+")
    range <- str_extract(headers, "\\d+\\.\\.\\d+")
    start <- as.integer(sub("\\.\\..*", "", range))
    end <- as.integer(sub("^.*\\.\\.", "", range))
    protein <- str_extract(headers, "\\|\\s*(.*)$")
    protein <- str_trim(str_sub(protein, str_locate(protein, "\\|")[[1]] + 1)) # Trim pipe and whitespace
    """

    # Convert R data to Julia DataFrame
    accessions = try rcopy(R"accession") catch _ missing end
    starts = try rcopy(R"start") catch _ missing end
    ends = try rcopy(R"end") catch _ missing end
    proteins = try rcopy(R"protein") catch _ missing end
    sequences = try rcopy(R"sapply(dna_seqs, as.character)") catch _ missing end

    # Handle missing values by replacing with empty arrays if necessary
    accessions = ismissing(accessions) ? String[] : accessions
    starts = ismissing(starts) ? Int[] : starts
    ends = ismissing(ends) ? Int[] : ends
    proteins = ismissing(proteins) ? String[] : proteins
    sequences = ismissing(sequences) ? String[] : sequences

    # Create a DataFrame from the extracted data
    fasta_df = DataFrame(
        Protein = replace.(proteins, r"[()]" => ""),
        Start = starts,
        End = ends,
        Sequence = sequences
    )

    println("Fasta file processed successfully.")
    println("Extracted data from fasta file:")
    println(fasta_df)
    return fasta_df
end

# Parse command-line arguments
function parse_command_line_args()
    s = ArgParseSettings(description = "A script for processing variant data.")
    @add_arg_table s begin
        "--ref"
        help = "Path to the reference fasta file"
        required = true
        arg_type = String
        "--in"
        help = "Path to the input variant file"
        required = true
        arg_type = String
        "--by"
        help = "Mode for filtering: 's' for non-synonymous variants only, 'p' for a separate file for each protein in the reference .fasta, 'b' for running both modes sequentially"
        required = true
        arg_type = String
        "--out"
        help = "Path to the output directory"
        arg_type = String
        default = ""
    end
    return parse_args(Base.ARGS, s)
end

function run_mode_s(fasta_path::String, variant_path::String)
    fasta_df = read_fasta_to_dataframe(fasta_path)
    variants = read_out_file(variant_path)
    variants = call_proteins(variants, fasta_df)
    convert_variant_columns!(variants)
    variants = substitute_variants(variants)
    variants = find_original_codons(variants)
    variants = find_variant_codons(variants)
    variants = translate_codons(variants)

    # filter by synonymity 
    ns_variants = filter(row -> row[:is_synonymous] == "No", variants)

    return filter_outfile(ns_variants, variant_path)
end

function main()
    args = parse_command_line_args()
    fasta_path = args["ref"]
    variant_path = args["in"]
    mode = args["by"]
    out_dir = args["out"] == "" ? dirname(variant_path) : args["out"]

    if mode == "s"
        run_mode_s(fasta_path, variant_path)
    elseif mode == "p"
        protein_info = pull_fasta(fasta_path)
        split_out_file(variant_path, protein_info, out_dir)
    elseif mode == "b"
        ns_variant_path = run_mode_s(fasta_path, variant_path)
        protein_info = pull_fasta(fasta_path)
        split_out_file(ns_variant_path, protein_info, out_dir)
    else
        println("Invalid mode. Use 's' for non-synonymous variants only, 'p' for a separate file for each protein in the reference .fasta, 'b' for running both modes sequentially")
    end
end

main()
