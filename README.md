# scripts

somewhere to stick scripts

***Bash scripts:***

**trim_fasta.sh**
cuts fasta sequences to between a start and end point. 

**Usage**
```
./trim_fasta.sh -i -o -s -e -z`
```

**Options**
- -i a path to the input fasta file to be cut
- -o a path to the location in which to save the trimmed fasta 
- -s the start point
- -e the end point
- -z (y/n) whether to adjust for zero indexing -  defaults to n


**extract_bic.sh**
pulls BICs from VeTrans output and creates a simple ggplot

This script takes one input:
- a path to the folder containing the VeTrans output

Usage
```
./extract_bic/sh <path>
```

**find_largest_file.sh**
finds the largest file in a folder. Useful for checking files before uploading to Github

This script takes one input:
- a path to the folder to check

Usage
```
.find_largest_file.sh <path>
```

**bwa.sh**
performs filter and mem steps of bwa in a single step 

This script takes two inputs:
- -r a path to the reference genome to be indexed 
- -s a path to the single-read file to be aligned
OR
- -p two paths to the paired-end files to be aligned

**Usage**
```
$0 -r <reference.fasta> -s <read.fq.gz> | -p <read1.fq> <read2.fq>
```


***Julia scripts:***

**haps_to_csv.jl**
converts Inferred_haplotypes.out output from find_haplotypes_multi or Inference_n_0.out output from reconstruct_haplotypes_multi to a .csv for easy plotting

**Usage**
```
julia haps_to_csv.jl --in <path> --times <path>
```

**Options**
- --in a path to the single_locus_trajectories.out file to be split
- --times a path to a Times.in file, also required by sl_traj

**pull_reading_frame_fastas** 
creates two .fasta files from reading_frames.dat files generated by reading_frames - one containing nucleotide reading frames and the other containing translated protein sequences. 

**Usage**
```
julia pull_reading_frame_fastas.jl --in <path>
```

**Options**
- --in a path to the reading_frames.dat file to be split

**split_sl_traj**
splits single_locus_trajectory.out files from [[sl_traj]] by protein, synonymity, or both

**Usage**
```
julia split_sl_traj.jl --in <path> --ref <path> --by <s, p, or b>
```

**Options**
- --in a path to the single_locus_trajectories.out file to be split
- --ref a path to a .fasta file of reference coding regions
- --by 's' to split by synonymity, 'p' to split by protein, or 'b' to split by both
- --out (optional) a path to an output folder, defaults to the folder containing the input file

***R scripts:***

**Blanche_plot.R**
creates a simple plot for Blanche output

This script takes three inputs:
- --in a path to the output_points.dat file output by Blanche
- --out a path to save the output file as a .jpeg
- --names a path to a .txt file containing sample labels. This must be the length of the number of samples. The plot will colour points by this variable. 

Usage
```
./Blanche_plot.R --in --out --names
```