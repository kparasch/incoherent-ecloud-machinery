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

pinch_folder=/eos/user/k/kparasch/Pinches/

MTI=1.0
MLI=1.0
DTO=1.0
DLO=1.0
do_symm=0

#### Argument Parser ####
for i in "$@"
do
case $i in
    --pinch_name=*)
    pinch_name="${i#*=}"
    shift
    ;;
    --MTI=*)
    MTI="${i#*=}"
    shift
    ;;
    --MLI=*)
    MLI="${i#*=}"
    shift
    ;;
    --DTO=*)
    DTO="${i#*=}"
    shift
    ;;
    --DLO=*)
    DLO="${i#*=}"
    shift
    ;;
    --do_symm)
    do_symm=1
    shift
    ;;
    *)

    ;;
esac
done
#########################

out_name=refined_${pinch_name}_MTI${MTI}_MLI${MLI}_DTO${DTO}_DLO${DLO}_do_symm${do_symm}

echo output name will be $out_name

if [[ -d "$out_name" ]]
then
    echo "ERROR: $out_name exists, exiting.."
    exit 1
fi    

mkdir $out_name

cp -r RefinePinch_source/* $out_name/


echo "============== start submit file ============="
tee temp_submit_file.sub << EOF
executable  = RefinePinch_source/refine_exec.sh
arguments = $pinch_name $MTI $MLI $DTO $DLO $do_symm $pinch_folder
output = $out_name/$pinch_name.out
error = $out_name/$pinch_name.err
log = $out_name/$pinch_name.log
transfer_input_files = $out_name
RequestCpus = 24
stream_output = True
stream_error = True
+BigMemJob = True
+MaxRunTime = 259200
queue
EOF
echo "=============== end submit file =============="

condor_submit temp_submit_file.sub
rm temp_submit_file.sub

