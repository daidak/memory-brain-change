#!/bin/bash 

#SBATCH --job-name=sum_norm_rep
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=05:30:00
#SBATCH --account=p274
#SBATCH --output scripts/normative_scripts/logs/slurm-%j.txt
##SBATCH --partition=bigmem



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


ds=$1
model_name=$2
nreps=$3
phenotypes=$4
wd=$5

echo "$ds $model_name $nreps $phenotypes $wd"


if [ -e scripts/normative_scripts/logs/slurm.summarize_normative_reps_$ds.txt ]; then
	exit 1
fi

mv scripts/normative_scripts/logs/slurm-${SLURM_JOBID}.txt scripts/normative_scripts/logs/slurm.summarize_normative_reps_$ds.txt


Rscript scripts/normative_scripts/summarize_normative_reps.r $ds $model_name $nreps $phenotypes $wd
