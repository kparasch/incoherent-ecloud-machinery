#!/bin/bash

intensity=1.00
for ii in 0 1 3 4 5 6 7 8 9
#for ii in 0 1 2 3 4 5 6 7 8 9
do
    readlink -f /eos/user/e/ecincohe/Quadrupole_initial_states/injection_ArcQuad_MQF_intensity_${intensity}e11ppb_sey_1.30_VRF_8MV_${ii}/MP_state* > /eos/user/e/ecincohe/Quadrupole_initial_states/state_lists/injection_ArcQuad_MQF_intensity_${intensity}e11ppb_sey_1.30_VRF_8MV_${ii}
done
