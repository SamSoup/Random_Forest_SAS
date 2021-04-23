#!/bin/bash
#SBATCH -J genInput
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -p normal
#SBATCH -o genSamples.o%j
#SBATCH -e genSamples.e%j
#SBATCH -t 48:00:00
#SBATCH -A Recombining-sex-chro
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Sam.Su@utexas.edu

module load launcher boost Rstats





export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=commands

$LAUNCHER_DIR/paramrun
