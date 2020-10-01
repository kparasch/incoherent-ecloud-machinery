#!/bin/bash

myhome=/afs/cern.ch/user/k/kparasch/work/public/ecloud_sixtracklib

source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4cuda9/latest/x86_64-centos7-gcc62-opt/setup.sh
unset PYTHONHOME
unset PYTHONPATH
source $myhome/miniconda3/bin/activate ""
export PATH=$myhome/miniconda3/bin:$PATH

export PYTHONPATH=$myhome/packages/PyFRIENDS:$PYTHONPATH
export PYTHONPATH=$myhome/incoherent-ecloud-machinery/Tools:$PYTHONPATH

repository=/eos/project/e/ecloud-simulations/kparasch

which pip

job_name="wrong_job_name"
max_tau=0.4
ptau_max=0.

qx=62.270
qy=60.295
qprime=20
IMO=40
VRF=8

intensity=0
seyMB=1.30
seyMQ=1.30
MB=false
MQ=false
#### Argument Parser ####
for i in "$@"
do
case $i in
    --job_name=*)
    job_name="${i#*=}"
    shift
    ;;
    --qx=*)
    qx="${i#*=}"
    shift
    ;;
    --qy=*)
    qy="${i#*=}"
    shift
    ;;
    --qprime=*)
    qprime="${i#*=}"
    shift
    ;;
    --IMO=*)
    IMO="${i#*=}"
    shift
    ;;
    --VRF=*)
    VRF="${i#*=}"
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
    --MB)
    MB=true
    shift
    ;;
    --MQ)
    MQ=true
    shift
    ;;
    --seyMB=*)
    seyMB="${i#*=}"
    shift
    ;;
    --seyMQ=*)
    seyMQ="${i#*=}"
    shift
    ;;
    *)

    ;;
esac
done
#########################

if  [ "MB" == "true" ]
then
    case "$intensity" in
        0.30) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.30e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.35) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.35e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.40) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.40e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.45) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.45e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.50) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.50e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.55) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.55e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.60) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.60e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.65) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.65e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.70) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.70e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.75) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.75e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.80) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.80e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.85) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.85e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.90) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.90e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.95) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_0.95e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.00) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_1.00e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.05) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_1.05e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.10) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_1.10e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.15) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_1.15e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.20) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_1.20e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.25) MBpinch=--mb ${repository}/refined_Pinches/MB/refined_LHC_MB_450GeV_sey${seyMB}_1.25e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0) MBpinch=
        ;;
        *) echo "intensity not supported, exiting..." ; exit 1 ;;
    esac
else
    MBpinch=
fi

if  [ "MQ" == "true" ]
then
    case "$intensity" in
        0.30) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.30e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.35) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.35e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.40) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.40e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.45) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.45e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.50) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.50e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.55) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.55e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.60) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.60e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.65) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.65e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.70) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.70e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.75) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.75e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.80) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.80e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.85) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.85e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.90) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.90e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.95) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_0.95e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.00) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_1.00e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.05) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_1.05e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.10) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_1.10e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.15) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_1.15e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.20) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_1.20e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.25) MQFpinch=--mqf ${repository}/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey${seyMQ}_1.25e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0) MQFpinch=
        ;;
        *) echo "intensity not supported, exiting..." ; exit 1 ;;
    esac

    case "$intensity" in
        0.30) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.30e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.35) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.35e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.40) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.40e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.45) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.45e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.50) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.50e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.55) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.55e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.60) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.60e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.65) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.65e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.70) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.70e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.75) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.75e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.80) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.80e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.85) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.85e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.90) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.90e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0.95) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_0.95e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.00) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_1.00e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.05) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_1.05e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.10) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_1.10e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.15) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_1.15e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.20) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_1.20e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        1.25) MQDpinch=--mqd ${repository}/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey${seyMQ}_1.25e11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 ;;
        0) MQDpinch=
        ;;
        *) echo "intensity not supported, exiting..." ; exit 1 ;;
    esac

else
    MQFpinch=
    MQDpinch=
fi

line=$repository/EC_line_v1.0/EC_line_v1.0_qx${qx}_qy${qy}_qprime${qprime}_IMO${IMO}A_VRF${VRF}MV.pkl
output_file=$job_name
copy_destination=$repository/Tracking_Data/DA_450GeV
tracking_arguments="--copy_destination ${copy_destination} --simulation_input ${line} \
--ptau_max ${ptau_max} --max_tau ${max_tau} --output ${job_name}.h5 \
${MBpinch}$ ${MQFpinch} ${MQDpinch}"

if [[ -f "${copy_destination}/${job_name}.h5" ]]
then
    echo "ERROR: ${copy_destination}/${job_name}.h5 exists, exiting.."
    exit 1
fi    
if [[ -d "simulations/$job_name" ]]
then
    echo "ERROR: simulations/$job_name exists, exiting.."
    exit 1
fi    
if [[ "$job_name" == "wrong_job_name" ]]
then
    echo "ERROR: $job_name is wrong, exiting.."
    exit 1
fi    

mkdir simulations/$job_name

cp ../Trackers/ecloud_track_v1.0.py       simulations/$job_name/
cp ../Tools/kostas_filemanager.py         simulations/$job_name/
cp ../Tools/distribution.py               simulations/$job_name/
cp ../Tools/RF_bucket.py                  simulations/$job_name/
cp ../Tools/normalization.py              simulations/$job_name/
cp ../Tools/ecloud_sixtracklib_helpers.py simulations/$job_name/

echo "============== start executable file ============="
tee simulations/$job_name/$job_name.sh << EOF
#!/bin/bash

myhome=/afs/cern.ch/user/k/kparasch/work/public/ecloud_sixtracklib
source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4cuda9/latest/x86_64-centos7-gcc62-opt/setup.sh
unset PYTHONHOME
unset PYTHONPATH
source \$myhome/miniconda3/bin/activate ""
export PATH=\$myhome/miniconda3/bin:\$PATH

rm -r ~/.nv/ComputeCache

cd $job_name
echo \$1
time python ecloud_track_v1.0.py ${tracking_arguments} --device 'opencl:0.0' --seed \$1

EOF
echo "=============== end executable file =============="

chmod +x simulations/$job_name/$job_name.sh

#$job_name/job.sh
echo "============== start submit file ============="
tee simulations/$job_name/submit_file.sub << EOF
executable  = simulations/$job_name/$job_name.sh
arguments = \$(ClusterId)\$(ProcId)
output = simulations/$job_name/htcondor.out
error = simulations/$job_name/htcondor.err
log = simulations/$job_name/htcondor.log
transfer_input_files = simulations/$job_name
request_GPUs = 1
request_CPUs = 1
+MaxRunTime = 86400
queue
EOF
echo "=============== end submit file =============="

condor_submit simulations/$job_name/submit_file.sub
