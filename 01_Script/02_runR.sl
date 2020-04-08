#!/bin/bash

#SBATCH --time=120:00:00         # Walltime
#SBATCH --nodes=1               # Use 1 Node (Unless code is multi-node parallelized)
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=24
#SBATCH -o slurm-%j.out
#SBATCH --job-name=FEISMCsp

# Export file and working directory
export FILENAME=$HOME/FEIS_replication/01_Script/02_FEIS_Simulation_Supplementary-Analyses.R

export WORK_DIR=$HOME/FEIS_replication/02_Data

# Load the default version of R
module load R

# Take advantage of all the threads (linear algebra)
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# Create scratch & copy everything over to scratch
# mkdir -p $SCR_DIR
# cd $SCR_DIR
# cp -p $WORK_DIR/* . 

cd $WORK_DIR

# Run the R script in batch
Rscript $FILENAME > $FILENAME.out

# Copy results over + clean up
# cd $WORK_DIR
# cp -pR $SCR_DIR/* .
# rm -rf $SCR_DIR

echo "End of program at `date`"
