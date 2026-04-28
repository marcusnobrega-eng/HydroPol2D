#!/bin/bash
#SBATCH --job-name=stanfbig
#SBATCH --output=stanfbig_%j.out
#SBATCH --error=stanfbig_%j.err
#SBATCH --partition=gpu
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gpus=1
#SBATCH --mail-type=END,FAIL

module load matlab

cd "$SLURM_SUBMIT_DIR"

echo "==== GPU CHECK ===="
nvidia-smi

matlab -batch "run('HP2D_Wrapper.m')"