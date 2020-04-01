#!/bin/bash

source /home/kparasch/workspace/ecloud_sixtracklib/miniconda3/bin/activate ""

export PYTHONPATH=/home/kparasch/workspace/ecloud_sixtracklib/packages/PyFRIENDS:$PYTHONPATH
export PYTHONPATH=/home/kparasch/workspace/ecloud_sixtracklib/incoherent-ecloud-machinery/Tools:$PYTHONPATH

which pip

pinch_name=Pinch$1
workers=2
nPinches=3
final_destination=/home/kparasch/workspace/ecloud_sixtracklib/incoherent-ecloud-machinery/Average_pinch


if [[ -d "$pinch_name" ]]
then
    echo "ERROR: $pinch_name exists, exiting.."
    exit 1
fi    

mkdir $pinch_name

cp -r AverPinch_source/* Pinch$1/
cp Simulation_parameters.py Pinch$1/

cd Pinch$1

python prepare_pinch.py Pinch$1 $workers $nPinches $final_destination

cd ..
