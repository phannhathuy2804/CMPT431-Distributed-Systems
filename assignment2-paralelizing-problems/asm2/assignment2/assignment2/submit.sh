#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun python ./test_scripts/curve_area_tester.pyc --execPath=./curve_area_parallel --scriptPath=./test_scripts/curve_area_evaluator.pyc
srun python ./test_scripts/heat_transfer_tester.pyc --execPath=./heat_transfer_parallel --scriptPath=./test_scripts/heat_transfer_evaluator.pyc
