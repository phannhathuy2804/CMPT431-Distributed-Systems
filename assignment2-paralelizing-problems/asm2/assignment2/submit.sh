#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun /home/nhathuyp/sfuhome/Study/CMPT431/asm2/assignment2/curve_area --nPoints 10000000 --coeffA 2.0 --coeffB 4.0  --rSeed 15