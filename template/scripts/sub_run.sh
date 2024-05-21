#!/bin/bash

# clean run directory and link all required files
./prep_run.sh

# set job id
JOBNO=JOB_NUMBER

# Set the source code directory
JULIADEPOT="/work/n02/n02/aleey/wavi/WAVIhpc/.julia_AlexEKI"

# record start times
TIMEQSTART="$(date +%s)"
echo Start-time `date` >> times


# submit the job
sbatch -J EKI_$JOBNO \
       -A $HECACC \
       --export HECACC=$HECACC,TIMEQSTART=$TIMEQSTART,IMGNAME=$IMGNAME,JDEPOT=$JULIADEPOT,JOBNO=$JOBNO,TIMEQSTART=$TIMEQSTART \
       ./run_repeat.sh

