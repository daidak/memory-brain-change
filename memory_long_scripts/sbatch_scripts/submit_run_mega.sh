#!/bin/bash 

#SBATCH --job-name=run_mega
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=50:00:00
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


feature=$1
wd=$2

echo "$feature $wd"


if [ -e logs/slurm.mega.$feature.txt ]; then
	exit 1
fi

mv logs/slurm-${SLURM_JOBID}.txt logs/slurm.mega.$feature.txt


Rscript /cluster/projects/p274/projects/p039_image_brain/scripts/memory_long_scripts/sbatch_scripts/run_mega.R $feature $wd

