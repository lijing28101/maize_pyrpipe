#!/bin/bash
#SBATCH -J jobname
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 24:00:00
#SBATCH -C EGRESS
#SBATCH --mail-user=usingh@iastate.edu
#SBATCH --mail-type=ALL

# move to your appropriate pylon5 directory
cd /pylon5/mc5pl7p/usingh/urmi/pyrpipeTest/maizeorph/maize_pyrpipe 

#load required modules
source /pylon5/mc5pl7p/usingh/lib/myAnacondaInstallation/etc/profile.d/conda.sh
source ~/.bashrc
conda activate maizeorphan 

# run commands
python maize_orphan_prediction.py


