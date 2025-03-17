#!/bin/bash 

#SBATCH --job-name=run_apoe
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=100:45:00
#SBATCH --account=p274
#SBATCH --output logs/slurm-%j.txt


 
######################
# setting environment
######################
echo "SETTING UP COLOSSUS ENVIRONMENT"
echo "LOADING SINGULARITY MODULE"
module purge
module load R/4.2.1-foss-2022a
echo `which R`
#module load matlab
#echo `which matlab`

#echo "SOURCING FREESURFER"
#export FREESURFER_HOME=/cluster/projects/p274/tools/mri/freesurfer/current
#source $FREESURFER_HOME/SetUpFreeSurfer.sh
#echo "SOURCING FSL"
#FSLDIR=/cluster/projects/p274/tools/mri/fsl/current
#. ${FSLDIR}/etc/fslconf/fsl.sh
#PATH=${FSLDIR}/bin:${PATH}
#export FSLDIR PATH
export LANG=en_US.utf8


wd=$1

echo "$wd"


if [ -e logs/slurm.apoe.mem.txt ]; then
	exit 1
fi

mv logs/slurm-${SLURM_JOBID}.txt logs/slurm.apoe.mem.txt


Rscript /cluster/projects/p274/projects/p039_image_brain/scripts/memory_long_scripts/sbatch_scripts/run_apoe_memory.R $wd

