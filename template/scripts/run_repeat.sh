#!/bin/bash 
################################################################################
# Run the model for as long as we can, then prepare for a restart and submit the next job.
################################################################################

#SBATCH --partition=standard
##SBATCH --qos=standard
#SBATCH --qos=taskfarm
##SBATCH --qos=short
##SBATCH --reservation=shortqos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
##SBATCH --time=00:20:00
##SBATCH --account=$HECACC
#SBATCH --chdir=../run
#SBATCH --no-requeue


# pass the julia depot path through singularity
export SINGULARITYENV_JULIA_DEPOT_PATH="/opt/julia"

# function to return number of seconds left in this job
# the squeue command returns either hh:mm:ss or mm:ss
# so handle both cases.
# We should add in 1-00:00:00 for a day

function hmsleft()
{
        local lhms
        lhms=$(squeue  -j $SLURM_JOB_ID -O TimeLeft | tail -1)
        echo $lhms
}
function secsleft() {
    if [[ ${#hms} < 6 ]]
    then
        echo secs=$(echo $hms|awk -F: '{print ($1 * 60) + $2 }')
    else
        echo secs=$(echo $hms|awk -F: '{print ($1 * 3600) + ($2 * 60) + $3 }')
    fi
}

# start timer
timeqend="$(date +%s)"
elapsedqueue="$(expr $timeqend - $TIMEQSTART)"
timestart="$(date +%s)"
echo >> times
echo Queue-time seconds $elapsedqueue >> times
echo Run start `date` >> times
hms=$(hmsleft)
echo Walltime left is $hms>>walltime
rem_secs=$(secsleft)  # function above
echo Walltime left in seconds is $rem_secs >> walltime
# Subtract 3 minutes
RUNTIME="$(($rem_secs-180))"
echo Will run for $RUNTIME sec >> walltime


echo "received from SLURM"  HECACC=$HECACC,IMGNAME=$IMGNAME,JDEPOT=$JDEPOT,JOBNO=$JOBNO,TIMEQSTART=$TIMEQSTART

# Run the job but leave 3 minutes at the end
timeout $RUNTIME singularity exec -B ${JDEPOT}:/opt/julia,${MNT_CASES},${W_CASES},$(pwd) ${IMGNAME} julia driver.jl


# Get the exit code
OUT=$?
echo 'job chain: leave_time activated, exit code' $OUT

# end timer
timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo >> times
echo Run end `date` >> times
echo Run-time seconds $elapsedtotal >> times

if [ $OUT == 0 ]; then

  # Simulation completed

  echo 'job chain: finished'

elif [ $OUT == 124 ]; then
  echo 'job chain did not finish'
  # Ran out of time, resubmit job
  TIMEQSTART="$(date +%s)"
    # Find the most recent pickup file
  unset -v PICKUP_FILE
  if ! ls PChkpt_*.jld2 1> /dev/null 2>&1 ; then
    echo 'job chain: fail, no pickup files'
    exit 1
  fi
  for file in PChkpt_*.jld2; do
    [[ $file -nt $PICKUP_FILE ]] && PICKUP_FILE=$file
  done
  # Extract the middle bit of this filename
  PICKUP=${PICKUP_FILE#PChkpt_}
  PICKUP=${PICKUP%.jld2}

  re='^[0-9]+$'
  if [[ $PICKUP =~ $re ]]; then

    echo 'job chain: pickup from permanent checkpoint'

    # Save the timestep, with any leading zeros removed
    NITER0=$(echo $PICKUP | sed 's/^0*//')

  else
  echo 'job chain: fail, problem w pickup (WAVI only currently supports pCkpt pickups)' $PICKUP
    exit 1
  fi

  #edit the driver namelist: replace the first instance of niter0 = * with appropriate checkpoint number
  NITER0_LINE="niter0 = $NITER0"
  echo $NITER0_LINE
  sed -i '0,/.*niter0.*/s//'"$NITER0_LINE"'/' driver.jl

  #submit the next job
  sbatch --job-name=EKI_$JOBNO --account=$HECACC --export=HECACC=$HECACC,IMGNAME=$IMGNAME,JDEPOT=$JDEPOT,JOBNO=$JOBNO ../scripts/run_repeat.sh

else
  echo 'job chain: fail, simulation died, exit code' $OUT
fi







