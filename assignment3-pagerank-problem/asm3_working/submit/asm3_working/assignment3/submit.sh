#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

# srun ./page_rank_pull_parallel --nThreads 4 --nIterations 10 --inputFile "./input_graphs/roadNet-CA"
# srun ./page_rank_push_parallel --nThreads 4 --nIterations 10 --inputFile "./input_graphs/roadNet-CA"
srun python scripts/page_rank_tester.pyc --execPath=./ --scriptPath=scripts/page_rank_evaluator.pyc --inputPath=/scratch/input_graphs/
