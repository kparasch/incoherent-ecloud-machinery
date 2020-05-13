#!/bin/bash

myhome=/afs/cern.ch/user/k/kparasch/work/public/ecloud_sixtracklib

source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4cuda9/latest/x86_64-centos7-gcc62-opt/setup.sh
unset PYTHONHOME
unset PYTHONPATH
source $myhome/miniconda3/bin/activate ""
export PATH=$myhome/miniconda3/bin:$PATH

export PYTHONPATH=$myhome/packages/PyFRIENDS:$PYTHONPATH
export PYTHONPATH=$myhome/incoherent-ecloud-machinery/Tools:$PYTHONPATH

repository=/eos/user/e/ecincohe

which pip

job_name="wrong_job_name"
IMO=0
skip_turns=10000
turns_per_checkpoint=1000000
last_checkpoint=10
max_tau=0.4
ptau_max=0.

intensity=0
do_ecloud=
#### Argument Parser ####
for i in "$@"
do
case $i in
    --job_name=*)
    job_name="${i#*=}"
    shift
    ;;
    --turns_per_checkpoint=*)
    turns_per_checkpoint="${i#*=}"
    shift
    ;;
    --skip_turns=*)
    skip_turns="${i#*=}"
    shift
    ;;
    --IMO=*)
    IMO="${i#*=}"
    shift
    ;;
    --last_checkpoint=*)
    last_checkpoint="${i#*=}"
    shift
    ;;
    --ptau_max=*)
    ptau_max="${i#*=}"
    shift
    ;;
    --intensity=*)
    intensity="${i#*=}"
    shift
    ;;
    --ecloud)
    ecloud=--ecloud
    shift
    ;;
    *)

    ;;
esac
done
#########################

case "$intensity" in
0.7) pinch=$repository/Pinches/refined_LHC_ArcDip_1.35sey_0.7e11ppb_symm2D_MTI4.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
1.2) pinch=$repository/Pinches/refined_LHC_ArcDip_1.35sey_1.2e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
0) pinch=nopinch ;;
*) echo "intensity not supported, exiting..." ; exit 1 ;;
esac

line_folder=$repository/Lines/Line_IMO_${IMO}
output_file=$job_name
copy_destination=$repository/Tracking_Data
tracking_arguments="${ecloud} --copy_destination ${copy_destination} --line_folder ${line_folder} \
--pinch ${pinch} --ptau_max ${ptau_max} --max_tau ${max_tau} --output ${job_name}.h5 \
--skip_turns ${skip_turns} --turns_per_checkpoint ${turns_per_checkpoint} \
--last_checkpoint ${last_checkpoint}"

if [[ -f "${copy_destination}/${job_name}.h5" ]]
then
    echo "ERROR: ${copy_destination}/${job_name}.h5 exists, exiting.."
    exit 1
fi    
if [[ -d "$job_name" ]]
then
    echo "ERROR: $job_name exists, exiting.."
    exit 1
fi    
if [[ "$job_name" == "wrong_job_name" ]]
then
    echo "ERROR: $job_name is wrong, exiting.."
    exit 1
fi    
if [[ -d "$job_name" ]]
then
    echo "ERROR: $job_name exists, exiting.."
    exit 1
fi    

python initialize_file.py ${copy_destination}/${job_name}.h5

mkdir $job_name

cp ../InjectionLine/009_track_for_long.py $job_name/
cp ../Tools/kostas_filemanager.py         $job_name/
cp ../Tools/distribution.py               $job_name/
cp ../Tools/RF_bucket.py                  $job_name/
cp ../Tools/normalization.py              $job_name/
cp ../Tools/ecloud_sixtracklib_helpers.py $job_name/


echo "============== start executable file ============="
tee $job_name/$job_name.sh << EOF
#!/bin/bash

myhome=/afs/cern.ch/user/k/kparasch/work/public/ecloud_sixtracklib
source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4cuda9/latest/x86_64-centos7-gcc62-opt/setup.sh
unset PYTHONHOME
unset PYTHONPATH
source \$myhome/miniconda3/bin/activate ""
export PATH=\$myhome/miniconda3/bin:\$PATH

cd $job_name
echo \$1
time python 009_track_for_long.py ${tracking_arguments} --device 'opencl:0.0' --seed \$1

EOF
echo "=============== end executable file =============="

chmod +x $job_name/$job_name.sh

#$job_name/job.sh
echo "============== start submit file ============="
tee $job_name/submit_file.sub << EOF
executable  = $job_name/$job_name.sh
arguments = \$(ClusterId)\$(ProcId)
output = $job_name/htcondor.out
error = $job_name/htcondor.err
log = $job_name/htcondor.log
transfer_input_files = $job_name
request_GPUs = 1
request_CPUs = 1
+MaxRunTime = 518400
queue
EOF
echo "=============== end submit file =============="

condor_submit $job_name/submit_file.sub

#+JobFlavour = "espresso"
#stream_output = True
#stream_error = True
#+MaxRunTime = 259200

