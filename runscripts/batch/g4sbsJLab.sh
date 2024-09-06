#!/bin/bash
#SBATCH --job-name=g4sbs
#SBATCH --export=ALL
#SBATCH --partition=priority
#SBATCH --time=2:00:00

export homedir=$HOME

export nevents_file=$1
nevents=`cat $nevents_file`
export tempdir=$homedir/g4sbs/runscripts/batch

export MACRO_PATH=$homedir/g4sbs/scripts/
export MCMACRO=$tempdir/g4sbs_$RUN_NUMBER.mac

echo "/g4sbs/filename $tempdir/g4sbs_$RUN_NUMBER.root" > $MCMACRO
cat $tempdir/gep_12GeV2_2023_base.mac >> $MCMACRO
echo "/g4sbs/run $nevents" >> $MCMACRO

cd $tempdir
source ~/.bashrc

export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH
echo "****************** Execution Host = "
hostname
echo "****************** g4sbs Starting ..."
g4sbs $MCMACRO #>& /dev/null
echo "****************** g4sbs Finished!"

cp    ${tempdir}/"g4sbs_$RUN_NUMBER.root"   ${OUTPUT_DIR}/
rm -f ${tempdir}/"g4sbs_$RUN_NUMBER.root"
rm -f $MCMACRO
