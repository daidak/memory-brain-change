#!/bin/bash 

#SBATCH --job-name=compute_Zchange
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=00:15:00
#SBATCH --account=p274
#SBATCH --output scripts/normative_scripts/logs/slurm-%j.txt
##SBATCH --partition=hugemem



######################
# setting environment
######################
echo "SETTING UP COLOSSUS ENVIRONMENT"
echo "LOADING SINGULARITY MODULE"
module purge
module load R/4.1.0-foss-2021a
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


mmodels=$1
model_name=$2
ds=$3
phenotypes=$4

echo "$mmodels $model_name $ds $phenotypes"


if [ -e scripts/normative_scripts/logs/slurm.$ds.$model_name.$phenotypes.txt ]; then
	exit 1
fi

mv scripts/normative_scripts/logs/slurm-${SLURM_JOBID}.txt scripts/normative_scripts/logs/slurm.$ds.$model_name.$phenotypes.txt


Rscript scripts/normative_scripts/compute_Zchange.r $mmodels $model_name $phenotypes
