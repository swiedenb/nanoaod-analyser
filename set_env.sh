# this is the aachen 3a nanoaod analyser 
source /cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-centos7-gcc7-opt/setup.sh
export MY_ANALYSIS_PATH=$PWD

echo "Set analysis path to "
echo $MY_ANALYSIS_PATH

export LHAPDF=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-omkpbe2/
export LHAPATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:$LHAPDF/share/LHAPDF/
export LD_LIBRARY_PATH=$LHAPDF/lib:$LD_LIBRARY_PATH

echo "Set LHAPATH path to "
echo $LHAPATH
