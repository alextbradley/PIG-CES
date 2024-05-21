#!/bin/bash

JOBNO=JOB_NUMBER

# move to run directory
cd ../run

export JULIADEPOT="/work/n02/n02/aleey/wavi/WAVIhpc/.julia_AlexEKI"

# copy the julia script from utilities
cp $W_ROOT/utilities/zip_mat.jl .

#execute the zipping script
export SINGULARITYENV_JULIA_DEPOT_PATH="/opt/julia"
singularity exec -B ${JULIADEPOT}:/opt/julia,$(pwd) $IMGNAME julia zip_mat.jl

# make target directory
EXPT_NAME=$(basename $(dirname $(dirname $(dirname $(pwd)))))
W_HOMEDIR=$W_HOMEROOT/ARCHER2_EKI/${EXPT_NAME}/EKI_${JOBNO}/run
ssh $W_HOMEHOST "mkdir -p $W_HOMEDIR"

rsync -avzL *.nc $W_HOMEHOST:$W_HOMEDIR

rm outfile.nc
