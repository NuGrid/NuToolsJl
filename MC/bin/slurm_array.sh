#!/bin/bash
# SLURM job array for MC/bin/run_mc.jl: each array task runs its own block of seeds into the shared OUT.
# Edit the #SBATCH lines (account/partition/time) and the variables below, then
#     sbatch slurm_array.sh
# Example: 10000 samples in 20 tasks of 500 seeds, 32 cores each.
#SBATCH --job-name=ppn_mc
#SBATCH --array=0-19
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_%A_%a.log

NUTOOLS=$HOME/NuToolsJl                          # this repository
TEMPLATE=$HOME/runs/mc_template/ppn              # finished-setup ppn run dir; its parent holds NPDATA
OUT=$HOME/runs/mc_10k                            # shared by all tasks
PER_TASK=500                                     # seeds per array task

first=$(( SLURM_ARRAY_TASK_ID * PER_TASK + 1 ))
last=$(( first + PER_TASK - 1 ))
julia "$NUTOOLS/MC/bin/run_mc.jl" --template "$TEMPLATE" --out "$OUT" --seeds "$first-$last" \
      --jobs "$SLURM_CPUS_PER_TASK"
