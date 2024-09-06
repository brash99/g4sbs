#!/bin/bash

#Initialize values
start=${1:-7001}
end=${2:-7001}
nevents=${3:-100}

echo $nevents > batch/nevents_test

#Generate the runlist
echo "Generating run list ..."
echo $start > batch/runlist_test
for (( i=$(($start + 1)) ; i<=$end ; i++ ))
do
echo "$i" >> batch/runlist_test
done
echo "Done"

#Simulation
echo "Submitting simulation jobs ... nevents = $nevents"
cd ~/g4sbs/runscripts/batch
./sendjobsjlab.sh runlist_test nevents_test
echo "Done"





