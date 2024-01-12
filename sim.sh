#!/bin/bash
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Jan22/setUp.sh
source /eos/user/u/ursovsnd/start_ali.sh

set -ux
echo "Starting script."
DIR=$1
ProcId=$2
LSB_JOBINDEX=$((ProcId+1))
NTOTAL=$4
NJOBS=$3


#echo $MUONS
#echo $DIR
#echo $SUB
N=$(( NTOTAL/NJOBS + ( LSB_JOBINDEX == NJOBS ? NTOTAL % NJOBS : 0 ) ))
FIRST=$(((NTOTAL/NJOBS)*(LSB_JOBINDEX-1)))

python $SNDSW_ROOT/shipLHC/run_simSND.py --PG --pID 211 --Estart 300 --Eend 301 -n $N --HX 1 --EVx -38. --EVy 44. --EVz 200.


xrdcp sndLHC* root://eosuser.cern.ch//eos/user/u/ursovsnd/private/SND_Data/"$DIR"/"$LSB_JOBINDEX"/sndLHC.Pythia8-TGeant4.root 
xrdcp geofile_full*.root root://eosuser.cern.ch//eos/user/u/ursovsnd/private/SND_Data/"$DIR"/"$LSB_JOBINDEX"/geofile_full.Pythia8-TGeant4.root
