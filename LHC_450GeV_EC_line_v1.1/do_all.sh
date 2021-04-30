#!/bin/bash

#qx=62.270
#qy=60.295
#qprime=20.0
#IMO=40.0
#VRF=8.0
qx=$1
qy=$2
qprime=$3
IMO=$4
VRF=$5

python 000_run_mask.py --noblock --qx0 ${qx} --qy0 ${qy} --qprime ${qprime} --I_MO ${IMO} --VRF400 ${VRF}
python 001_make_line.py --noblock 
python 002_setup_eclouds.py --noblock 
python 003_collimators.py  
python 004_final_input.py 
cp simulation_input.pkl Lines/EC_line_v1.1_qx${qx}_qy${qy}_qprime${qprime}_IMO${IMO}A_VRF${VRF}MV.pkl
