#!/bin/bash
#BSUB -P "espaloma"
#BSUB -J "hydration-openeye[1-642]"
#BSUB -n 1
#BSUB -R rusage[mem=8]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 12 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
#BSUB -W 24:00
#BSUB -o stdout/out_%J_%I.stdout
#BSUB -eo stdout/out_%J_%I.stderr
#BSUB -L /bin/bash

source ~/.bashrc
OPENMM_CPU_THREADS=1

echo "changing directory to ${LS_SUBCWD}"
cd "$LS_SUBCWD"
conda activate hydration

# export environment information to file
conda list

# Report node in use
hostname

# Open eye license activation/env
export OE_LICENSE=~/.openeye/oe_license.txt

# Report CUDA info
env | sort | grep 'CUDA'

# Report GPU info
nvidia-smi -L
nvidia-smi --query-gpu=name --format=csv

# Run
python hydration.py freesolv --index $LSB_JOBINDEX --toolkit OpenEye --method am1bccelf10 --forcefield "gaff-2.11" --filepath openeye --niterations 1000
