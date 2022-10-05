#!/bin/env sh

#~ eval "$(cat ./set_env.sh)"
source ./set_env.sh
mkdir output
mkdir output/2018
mkdir output/2017
mkdir output/2016
mkdir output/2018/snapshots
mkdir output/2017/snapshots
mkdir output/2016/snapshots
./guent.her "$@" || { echo 'running music failed' ; exit 1; }

tar cjf AnalysisOutput.tar.bz2 output/
#tar czf MusicOutDir.tar.gz AnalysisOutput
exit 0
