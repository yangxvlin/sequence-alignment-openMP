#!/bin/bash
#
#SBATCH --partition=physical
#SBATCH --output=parallel-out.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=xuliny@student.unimelb.edu.au
#SBATCH --mail-type=End
#SBATCH --time=0-0:10:00

# Load required modules
module load GCC/8.2.0

# Launch code
g++ -fopenmp -std=c++14 seqalign_parallel.cpp -o parallel
parallel.exe < seq.dat