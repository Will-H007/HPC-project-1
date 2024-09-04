#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=work
#SBATCH --account=courses0101
#SBATCH --mem=200G
#SBATCH --time=00:15:00

module load gcc

# Compile the C++ program with OpenMP support
g++ -std=c++17 -fopenmp -o hpc_project_1 ./hpc_project_1.cpp

if [ -f ./hpc_project_1 ]; then
    # Loop through thread counts in reverse order
    for threads in 64 32 16 8 4; do
        echo "Running with $threads threads"
        export NUM_THREADS=$threads
        srun --ntasks=1 --cpus-per-task=$threads ./hpc_project_1
    done
else
    echo "Compilation failed: Executable not found."
    exit 1
fi
