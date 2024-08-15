#!/usr/bin/env julia

using DataFrames, CSV, FASTX
using ArgParse

function process_haplotype_data(haps_file::String, out_file::String)
    haplo_df = CSV.read(haps_file, DataFrame)
    loci = []
    open(out_file, "r") do f
        for line in eachline(f)
            push!(loci, parse(Int, split(line)[1]))
        end
    end
    unique_loci = unique(loci)
    result_df = DataFrame(locus = unique_loci)
    for haplotype in haplo_df.Haplotype
        result_df[!, haplotype] = [haplotype[j] for j in 1:length(unique_loci)]
    end
    return result_df
end

function adjust_locus(dataframe::DataFrame, start_locus::Int = 0)
    adj_locus_values = dataframe.locus .- start_locus
    insertcols!(dataframe, 2, :adj_locus => adj_locus_values)
    return dataframe
end

function create_modified_fasta(haplotypes::DataFrame, fasta_path::String; output_fasta_path::String="", suffix::String="")
    if output_fasta_path == ""
        fasta_dir = dirname(fasta_path)
        fasta_filename = basename(fasta_path)
        fasta_base, fasta_ext = splitext(fasta_filename)
        output_fasta_path = joinpath(fasta_dir, "$(fasta_base)_$(suffix)_haplotypes$fasta_ext")
    end

    reader = FASTA.Reader(open(fasta_path, "r"))
    original_sequence = ""
    for record in reader
        original_sequence = String(FASTA.sequence(record))
    end
    close(reader)

    sequence_length = length(original_sequence)
    writer = open(output_fasta_path, "w")

    for col in 3:ncol(haplotypes)
        haplotype_name = names(haplotypes)[col]
        if all(haplotypes[:, col] .== 'X')
            continue
        end
        
        modified_sequence_chars = collect(original_sequence)

        for row in 1:nrow(haplotypes)
            position = haplotypes[row, :adj_locus]
            base = haplotypes[row, col]
            
            if position > sequence_length || position < 1
                println("Skipping out-of-bounds position $position for haplotype $haplotype_name")
                continue
            end

            modified_sequence_chars[position] = base
        end

        modified_sequence = join(modified_sequence_chars)
        println(writer, ">$haplotype_name")
        println(writer, modified_sequence)
    end

    close(writer)
    println("FASTA file successfully written to: $output_fasta_path")
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--haps", "-i"
        help = "Path to the input haplotypes CSV file"
        arg_type = String
        required = true
    end

    @add_arg_table! s begin
        "--out", "-o"
        help = "Path to the .out file"
        arg_type = String
        required = true
    end
    
    @add_arg_table! s begin
        "--start", "-s"
        help = "Starting locus to adjust"
        arg_type = Int
        default = 0
    end

    @add_arg_table! s begin
        "--fa", "-f"
        help = "Path to the input FASTA file"
        arg_type = String
        required = true
    end

    @add_arg_table! s begin
        "--out-fasta", "-O"
        help = "Path to the output FASTA file (optional, defaults to the same directory as the input FASTA)"
        arg_type = String
        default = ""
    end
    
    @add_arg_table! s begin
        "--suffix", "-x"
        help = "Suffix to add to the output FASTA file (optional)"
        arg_type = String
        default = ""
    end

    parsed_args = parse_args(s)

    # Correct access to the "out-fasta" key
    output_fasta_path = parsed_args["out-fasta"]
    suffix = parsed_args["suffix"]

    haps = process_haplotype_data(parsed_args["haps"], parsed_args["out"])
    haps = adjust_locus(haps, parsed_args["start"])
    
    create_modified_fasta(haps, parsed_args["fa"], output_fasta_path=output_fasta_path, suffix=suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0
        println("No arguments provided. Use --help or -h to see the usage instructions.")
        s = ArgParseSettings()
        @add_arg_table! s begin
            "--haps", "-i"
            help = "Path to the input haplotypes CSV file"
            arg_type = String
        end

        @add_arg_table! s begin
            "--out", "-o"
            help = "Path to the .out file"
            arg_type = String
        end

        @add_arg_table! s begin
            "--start", "-s"
            help = "Starting locus to adjust"
            arg_type = Int
        end

        @add_arg_table! s begin
            "--fa", "-f"
            help = "Path to the input FASTA file"
            arg_type = String
        end

        @add_arg_table! s begin
            "--out-fasta", "-O"
            help = "Path to the output FASTA file (optional, defaults to the same directory as the input FASTA)"
            arg_type = String
        end

        @add_arg_table! s begin
            "--suffix", "-x"
            help = "Suffix to add to the output FASTA file (optional)"
            arg_type = String
        end

        parse_args(s)
    else
        main()
    end
end