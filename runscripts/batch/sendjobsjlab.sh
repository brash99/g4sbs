#!/bin/bash

if [[ $# -ne 2 ]]; then
	echo "illegal number of parameters"
	exit 
fi

export runlist=$1
export nevents=$2
export homedir=$HOME
export OUTPUT_DIR=$homedir/g4sbs/data

echo "Reading runs from $runlist"
echo "Number of events per run: $nevents"

#submit a job for each file in filelist
for run in `cat $runlist`
do
    echo "The next run is $run"
    export RUN_NUMBER=$run
    sbatch g4sbsJLab.sh $nevents
#    qsub g4sbs.sh 
#    ./g4sbs.sh $nevents
    sleep 1
done

echo "All finished ... :)"
