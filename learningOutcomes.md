## Week 1
 - pros, cons and rationale for using HPC
 - logging into `cedar` or `graham` using `ssh`
 - filesystem and policies
 - UNIX commands
  - `pwd`
  - `ls` and options `-l -a -t -r -h`
  - `cd` with `..`, `-`, `/`
  - `mkdir`, `rmdir`, `rm -rf` (careful with the last one)
  - `wget`
 - how to locate genome fasta and Gencode .gtf files from the Internet

## Week 2 
 - `gunzip`ping genome fasta and Gencode files
 - editing ~/.bashrc with `nano`, saving and exiting
  - adding custom path as an environment variable
 - `module` command with `load, list, spider` to load pre-installed software
 - submitting a batch job with slurm with `sbatch`
  - writing a job script with necessary `#!/bin/bash` and `#SBATCH` lines
 - understanding `genomeGenerate` script with `STAR` aligner, and submitting to slurm
 - checking status of your job submission with `sq`, `sacct` and `seff`

## Week 3
 - using FileZilla client to transfer files from your local computer to HPC
 - concatenating fastq.gz files together with `cat`
 - learning how to use `fastqc` by Google searching
 - understanding and submitting `STAR` align script to slurm

## Week 4
 - checking alignment rate output from `STAR`
 - loading `HTSeq` with python `virtualenv`
 - understanding and submitting `htseq-count` script to slurm
