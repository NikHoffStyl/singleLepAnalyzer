#!/bin/bash

condorDir=$PWD
theDir=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh
sleep 5
source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-centos7-gcc8-opt/setup.sh
sleep 5
source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-centos7-gcc8-opt/setup.sh

cd $theDir
eval `scramv1 runtime -sh`
pwd
python -u doHists.py $condorDir \
					--iPlot=${2} \
					--region=${3} \
					--isCategorized=${4} \
					--isEM=${5} \
					--nhott=${6} \
					--nttag=${7} \
					--nWtag=${8} \
					--nbtag=${9} \
					--njets=${10} \
					