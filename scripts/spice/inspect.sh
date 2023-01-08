#BSUB -o %J.stdout
#BSUB -R "rusage[mem=50]"
#BSUB -W 7:59
#BSUB -n 1

python _inspect.py

