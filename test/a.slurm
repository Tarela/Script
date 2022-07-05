#!/bin/bash
#SBATCH -n 2
#SBATCH --mem=1000
#SBATCH -t 15:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o a.out
#SBATCH -e a.err
#Run program
python a.py > a.txt
