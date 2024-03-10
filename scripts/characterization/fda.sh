#BSUB -J "data[1-1038]"
#BSUB -o %J.stdout
#BSUB -W 0:10

python run.py --in_path fda.smi --idx $((LSB_JOBINDEX - 1)) --out fda.csv 

