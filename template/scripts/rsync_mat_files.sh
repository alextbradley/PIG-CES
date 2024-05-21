#!/bin/bash

JOBNO=JOB_NUMBER

# move to run directory
cd ../run

# make target directory
EXPT_NAME=$(basename $(dirname $(dirname $(dirname $(pwd)))))
W_HOMEDIR=$W_HOMEROOT/ARCHER2_EKI/${EXPT_NAME}/EKI_${JOBNO}/run
ssh $W_HOMEHOST "mkdir -p $W_HOMEDIR"

rsync -avzL *.mat $W_HOMEHOST:$W_HOMEDIR
