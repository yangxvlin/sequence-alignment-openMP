#!/bin/bash
#
#SBATCH --partition=physical
#SBATCH --output=sequential-out.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=xuliny@student.unimelb.edu.au
#SBATCH --mail-type=End
#SBATCH --time=0-0:10:00

# Load required modules
module load gcc/10.1.0

# Launch code
g++ -std=c++14 seqalign.cpp -o sequential -O3
./sequential < seq.dat