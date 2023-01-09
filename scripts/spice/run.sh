for learning_rate in 1e-3; do # 1e-4 1e-5; do
for weight_decay in 1e-4; do # ;1e-6 1e-8 1e-10; do
for batch_size in 512; do # ; 128 512; do
for n_epochs in 5000; do

    bsub -q gpuqueue -o %J.stdout -gpu "num=1:j_exclusive=yes" -R "rusage[mem=50] span[ptile=1]" -W 47:59 -R V100 -n 1\
    python train.py \
        --learning_rate $learning_rate \
        --weight_decay $weight_decay \
        --batch_size $batch_size \
        --n_epochs $n_epochs 

done; done; done; done
