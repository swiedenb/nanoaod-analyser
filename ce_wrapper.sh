#!/bin/env sh

#~ eval "$(cat ./set_env.sh)"
echo 'sourcing LCG'
source ./set_env.sh
echo 'Creating directories'
mkdir output
mkdir output/2018
mkdir output/2017
mkdir output/2016
mkdir output/2018/snapshots
mkdir output/2017/snapshots
mkdir output/2016/snapshots
echo 'starting the analysis'
./guent.her "$@" || { echo 'running music failed' ; exit 1; }

tar cjf AnalysisOutput.tar.bz2 output/
#tar czf MusicOutDir.tar.gz AnalysisOutput
exit 0
