#!/bin/bash
#SBATCH --partition=Centaurus
#SBATCH --job-name=nbody_bench
#SBATCH --output=benchmark_%j.out
#SBATCH --error=benchmark_%j.err
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G

module load g++

# Build the executable
make nbody_sim

# Loop over each input file and run sim
for input in input1.txt input2.txt input3.txt; do
    base_name=${input%.txt}  # strip .txt
    echo "Running benchmark with $input..."

    # Run sim and write to each tsv
    ./nbody_sim < $input > simulation_output_${base_name}.tsv

    # plot to generate pdf
    python3 plot.py simulation_output_${base_name}.tsv trajectories_${base_name}.pdf
done
