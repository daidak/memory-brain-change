#!/bin/bash 

#SBATCH --job-name=impute_outlier
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=02:00:00
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


outdir=$1
value=$2
m=$3
dataset=$4
model_name=$5

echo "$outdir $value $m $dataset $model_name"


if [ -e scripts/normative_scripts/logs/slurm.impute.outlier.$dataset.$value.$model_name.txt ]; then
	exit 1
fi

mv scripts/normative_scripts/logs/slurm-${SLURM_JOBID}.txt scripts/normative_scripts/logs/slurm.impute.outlier.$dataset.$value.$model_name.txt


Rscript scripts/normative_scripts/impute_outlier_data.r $outdir $value $m $model_name

