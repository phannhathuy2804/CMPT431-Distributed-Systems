#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00
#SBATCH --mem=1G
#SBATCH --partition=slow

# srun ./page_rank_pull_parallel --nThreads 4 --nIterations 100 --strategy 1 --granularity 10000  --inputFile "/scratch/input_graphs/roadNet-CA"
# srun ./page_rank_pull_parallel --nThreads 4 --nIterations 100 --strategy 2 --granularity 10000  --inputFile "/scratch/input_graphs/roadNet-CA"
# srun ./page_rank_pull_parallel --nThreads 4 --nIterations 100 --strategy 3 --granularity 10000  --inputFile "/scratch/input_graphs/roadNet-CA"
# srun ./page_rank_pull_parallel --nThreads 4 --nIterations 100 --strategy 4 --granularity 10000  --inputFile "/scratch/input_graphs/roadNet-CA"


srun ./page_rank_push_parallel_atomic --nThreads 4 --nIterations 20 --strategy 1 --granularity 10000  --inputFile "/scratch/input_graphs/web-Google"
# srun ./page_rank_push_parallel_atomic --nThreads 4 --nIterations 100 --strategy 2 --granularity 10000  --inputFile "/scratch/input_graphs/roadNet-CA"
# srun ./page_rank_push_parallel_atomic --nThreads 4 --nIterations 100 --strategy 3 --granularity 10000  --inputFile "/scratch/input_graphs/roadNet-CA"
srun ./page_rank_push_parallel_atomic --nThreads 4 --nIterations 20 --strategy 4 --granularity 1  --inputFile "/scratch/input_graphs/web-Google"
#srun python scripts/page_rank_tester.pyc --execPath=./ --scriptPath=scripts/page_rank_evaluator.pyc --inputPath=/scratch/input_graphs/
