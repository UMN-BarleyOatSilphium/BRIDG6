#!/bin/bash
#PBS -l walltime=2:00:00,mem=62gb,nodes=1:ppn=24
#PBS -m abe -M tnaik@umn.edu
#PBS -q mesabi
#PBS -N ibd
cd /home/pepsico/tnaik/Pepsi/Analysis/Identity/input
module load python3
module load pandas
module load numpy
module load multiprocessing
module load csv
python3 ./SNPCalls_Results_0.py 

