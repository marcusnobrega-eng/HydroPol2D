#!/bin/bash
#SBATCH --job-name=hydropol2d
#SBATCH --output=hydropol2d_%j.out
#SBATCH --error=hydropol2d_%j.err
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

matlab -batch "run(fullfile(getenv('SLURM_SUBMIT_DIR'),'Applications','HPC','HP2D_Wrapper.m'))"
