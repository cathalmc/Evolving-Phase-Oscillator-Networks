#!/bin/bash -l

# Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cores per node required
#SBATCH -n 1

# Set the walltime of the job to 10 mins (format is hh:mm:ss)
#SBATCH -t 24:00:00

# E-mail on begin (b), abort (a), and end (end) of job
#SBATCH --mail-type=ALL

# E-mail address of recipient
#SBATCH --mail-user=14369856@ucdconnect.ie

# Specify the jobname
#SBATCH --job-name=geneticAlg

# Specify the error and output file names
#SBATCH --error="ErrorSW.out"
#SBATCH --output="OutputSW.out"

# Setup the environment
module load anaconda
module load gcc openmpi/3.1.4
conda activate --stack mpynn4

# Change to model working directory
model_dir="${HOME}/storage/newGeneticAlg"

cd ${model_dir}

# Setup the model script
model_script=runSim.py
model="${model_dir}/${model_script}"

python ${model} 
