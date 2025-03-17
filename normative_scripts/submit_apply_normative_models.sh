#!/bin/bash 

#SBATCH --job-name=normative_model
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=04:00:00
#SBATCH --account=p274
#SBATCH --output scripts/normative_scripts/logs/slurm-%j.txt
##SBATCH --partition=hugemem



######################
# setting environment
######################
module purge
module load SciPy-bundle/2021.05-foss-2021a
module load Tkinter/3.9.5-GCCcore-10.3.0
source /cluster/projects/p274/projects/p039_image_brain_change/py3_env/bin/activate
echo `which python`

echo $PATH
echo "SETTING UP COLOSSUS ENVIRONMENT"
echo "LOADING SINGULARITY MODULE"

pip list
#module load R/3.6.3-foss-2020a
#echo `which R`
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

model_name=${1}
site_file=${2}
modeldir=${3}
phenotypes=${4}
dataset=${5}
mm=$(basename $modeldir)

echo "$model_name $site_file $modeldir $phenotypes"
mv scripts/normative_scripts/logs/slurm-${SLURM_JOBID}.txt scripts/normative_scripts/logs/slurm-normative-model-${dataset}-${model_name}-${mm}.txt

python scripts/normative_scripts/apply_normative_models_new_samples.py $model_name $site_file $modeldir $phenotypes
   
