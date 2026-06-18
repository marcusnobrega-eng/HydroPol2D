#!/bin/bash
#SBATCH --job-name=stan100
#SBATCH --output=stan100_%j.out
#SBATCH --error=stan100_%j.err
#SBATCH --partition=gpu
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gpus=1
#SBATCH --mail-type=END,FAIL

module load matlab

cd "$SLURM_SUBMIT_DIR"

echo "==== GPU CHECK ===="
nvidia-smi

matlab -batch "run('HP2D_Wrapper.m')"