#!/bin/bash
#SBATCH --job-name=test-cpu
#SBATCH --mem=8G               # max memory per node
#SBATCH --cpus-per-task=2      # max CPU cores per MPI process
#SBATCH --time=0-20:00         # time limit (DD-HH:MM)
#SBATCH --partition=rng-short  # rng-short is default, but use rng-long if time exceeds 7h
#SBATCH --nice=2000            # allow other priority jobs to go first (note, this is different from the linux nice command below)

# Worker script that runs one instance of ExploreASL on a (remote) terminal

if 0; then
    echo $XASLDIR
    echo $DATAFOLDER
fi 


# run Script
nice -n 10 matlab -nodesktop -nosplash -r "cd('$XASLDIR');ExploreASL('$DATAFOLDER', 0, 1);exit;"
echo "xASL has ran over subdirectory ${DATAFOLDER}" 

exit 0