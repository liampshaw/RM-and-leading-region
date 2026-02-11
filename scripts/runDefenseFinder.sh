#!/bin/bash
# Slurm array job to run defensefinder on Bristol BluePebble cluster
# Author: Liam Shaw

#SBATCH --job-name=defensefinder
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00 
#SBATCH --mem=200M
#SBATCH --account=panm038524
#SBATCH --array=1-1751 # if running on n=1,751 plasmids in the original dataset

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

hostname
. ~/initMamba.sh
conda activate defensefinder-2.0.1

# Use /tmp on node to run job
OUT_TMP=/tmp/defense_finder_tmp_${SLURM_ARRAY_TASK_ID}
OUT_FINAL=/user/home/xr24099/trieste/defense-finder-outputs
# Copy input across
mkdir -p $OUT_TMP
cp "$name" $OUT_TMP
cd $OUT_TMP
# run defensefinder
defense-finder run "$name" --workers ${SLURM_CPUS_PER_TASK} -o $OUT_TMP
# Copy files back
rsync -a --remove-source-files "${OUT_TMP}/" "${OUT_FINAL}/"
rmdir "${OUT_TMP}"
