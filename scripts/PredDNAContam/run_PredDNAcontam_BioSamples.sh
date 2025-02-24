#!/bin/bash

#SBATCH --account=development
#SBATCH --ntasks=36                
#SBATCH --mem=180G
#SBATCH --time=00:15:00                   
#SBATCH --qos=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com
#SBATCH --job-name=ML_Contamination_Model_test
#SBATCH --output=run_ml_model_RF_test_%j.out
#SBATCH --error=run_ml_model_RF_test_%j.err

# Load the conda environment
source /path/to/bin/miniconda3/etc/profile.d/conda.sh
conda activate /path/to/conda_env

# Define the path to the Python script
python_path="/path/to/python/script/run_PredDNAcontam_BioSamples.py"

# Run the Python script with the configuration file
python $python_path 