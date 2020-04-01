#!/bin/bash

myhome=/afs/cern.ch/user/k/kparasch/work/public/ecloud_sixtracklib

source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4cuda9/latest/x86_64-centos7-gcc62-opt/setup.sh
unset PYTHONHOME
unset PYTHONPATH
source $myhome/miniconda3/bin/activate ""
export PATH=$myhome/miniconda3/bin:$PATH

export PYTHONPATH=$myhome/packages/PyFRIENDS:$PYTHONPATH
export PYTHONPATH=$myhome/incoherent-ecloud-machinery/Tools:$PYTHONPATH

which pip

pinch_name=Pinch$1
workers=24
nPinches=3000
final_destination=$myhome/incoherent-ecloud-machinery/Average_pinch


if [[ -d "$pinch_name" ]]
then
    echo "ERROR: $pinch_name exists, exiting.."
    exit 1
fi    

mkdir $pinch_name

cp -r AverPinch_source/* Pinch$1/
cp Simulation_parameters.py Pinch$1/


echo "============== start submit file ============="
tee temp_submit_file.sub << EOF
executable  = AverPinch_source/aver_exec.sh
arguments = $pinch_name $workers $nPinches $final_destination
output = $pinch_name/$pinch_name.out
error = $pinch_name/$pinch_name.err
log = $pinch_name/$pinch_name.log
transfer_input_files = $pinch_name
RequestCpus = 24
stream_output = True
stream_error = True
+BigMemJob = True
+JobFlavour = "espresso"
queue
EOF
echo "=============== end submit file =============="

condor_submit temp_submit_file.sub
rm temp_submit_file.sub
#+JobFlavour = "nextweek"

