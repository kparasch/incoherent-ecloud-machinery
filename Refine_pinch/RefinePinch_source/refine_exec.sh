#!/bin/bash

myhome=/afs/cern.ch/user/k/kparasch/work/public/ecloud_sixtracklib

source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4cuda9/latest/x86_64-centos7-gcc62-opt/setup.sh
unset PYTHONHOME
unset PYTHONPATH
source $myhome/miniconda3/bin/activate ""
export PATH=$myhome/miniconda3/bin:$PATH

export PYTHONPATH=$myhome/packages/PyFRIENDS:$PYTHONPATH
export PYTHONPATH=$myhome/incoherent-ecloud-machinery/Tools:$PYTHONPATH

pinch_name=$1
MTI=$2
MLI=$3
DTO=$4
DLO=$5
do_symm=$6
pinch_folder=$7
out_name=refined_${pinch_name}_MTI${MTI}_MLI${MLI}_DTO${DTO}_DLO${DLO}_do_symm${do_symm}

cd $out_name

time python refine_pinch.py $pinch_name $MTI $MLI $DTO $DLO $do_symm $pinch_folder

cd ..

