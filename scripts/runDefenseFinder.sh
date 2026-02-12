#!/bin/bash
# Slurm array job to run defensefinder on Bristol BluePebble cluster
# Author: Liam Shaw

#SBATCH --job-name=defensefinder
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:10:00 
#SBATCH --mem=200M
#SBATCH --account=panm038524
#SBATCH --array=1-1751%50 # limit to 50 jobs simultaneously 

set -euo pipefail
# To reduce chance of launch failure 
sleep $((RANDOM % 15))

# Set thread values (unclear if it matters)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Get name of fasta file from samples.txt list 
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

hostname
# activate defensefinder (assumes we have already activated mamba with initMamba.sh on login node
source activate defensefinder-2.0.1

# Make temporary output folder
OUT_TMP=/tmp/defense_finder_tmp_${SLURM_ARRAY_TASK_ID}
OUT_FINAL=/user/home/xr24099/trieste/defense-finder-outputs
mkdir -p $OUT_TMP
# Check that file exists
if [ -z "$name" ]; then
    echo "No input file for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [ ! -f "$name" ]; then
    echo "File $name does not exist"
    exit 1
fi
# Then copy across if it does
cp "$name" $OUT_TMP

# Change to /tmp on node to run defense-finder
cd $OUT_TMP
defense-finder run -a "$name" --workers ${SLURM_CPUS_PER_TASK} -o $OUT_TMP
# Cp back tsv output
cp ${OUT_TMP}/*.tsv "${OUT_FINAL}/"
# Clean up
rm -r "${OUT_TMP}"
