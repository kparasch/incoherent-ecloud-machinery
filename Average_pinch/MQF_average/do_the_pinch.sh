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

pinch_name="wrong_pinch_name"
workers=24
nPinches=4000
final_destination=/eos/project/e/ecloud-simulations/kparasch/Pinches/MQF/


#### Argument Parser ####
for i in "$@"
do
case $i in
    --pinch_name=*)
    pinch_name="${i#*=}"
    shift
    ;;
    --states=*)
    states="${i#*=}"
    shift
    ;;
    --workers=*)
    workers="${i#*=}"
    shift
    ;;
    --nPinches=*)
    nPinches="${i#*=}"
    shift
    ;;
    --intensity=*)
    intensity="${i#*=}"
    shift
    ;;
    *)

    ;;
esac
done
#########################


if [[ -d "simulations/$pinch_name" ]]
then
    echo "ERROR: simulations/$pinch_name exists, exiting.."
    exit 1
fi    

mkdir simulations/$pinch_name

cp LHC_chm_ver.mat simulations/$pinch_name/
cp prepare_pinch.py simulations/$pinch_name/
cp -r pyecloud_config simulations/$pinch_name/
cp Simulation_parameters.py simulations/$pinch_name/
cp replaceline.py simulations/$pinch_name/
python3 replace_intensity.py simulations/$pinch_name/Simulation_parameters.py $intensity

echo "============== start submit file ============="
tee temp_submit_file.sub << EOF
executable  = aver_exec.sh
arguments = $pinch_name $workers $nPinches $final_destination $states
output = simulations/$pinch_name/htcondor.out
error = simulations/$pinch_name/htcondor.err
log = simulations/$pinch_name/htcondor.log
transfer_input_files = simulations/$pinch_name
RequestCpus = 24
+BigMemJob = True
+JobFlavour = "nextweek"
queue
EOF
echo "=============== end submit file =============="

condor_submit temp_submit_file.sub
rm temp_submit_file.sub

#+MaxRunTime = 259200
#+JobFlavour = "nextweek"
#stream_output = True
#stream_error = True

