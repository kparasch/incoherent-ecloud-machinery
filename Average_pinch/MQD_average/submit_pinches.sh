#!/bin/bash

states_folder=/eos/user/e/ecincohe/Quadrupole_initial_states/state_lists

intensity=0.70
for ii in 7
#for ii in 0 1 2 3 4 5 6 7 8 9
do
    ./do_the_pinch.sh --pinch_name=LHC_MQD_450GeV_sey1.30_${intensity}e11ppb_400npinches_${ii} --states=${states_folder}/injection_ArcQuad_MQF_intensity_${intensity}e11ppb_sey_1.30_VRF_8MV_${ii} --nPinches=400 --intensity=${intensity}
done

#./do_the_pinch.sh --pinch_name=LHC_ArcDipReal_450GeV_sey1.30_1.20e11ppb --edensity=/eos/user/e/ecincohe/e-cloud_buildup_profiles/edensity_LHC_ArcDipReal_450GeV_sey1.30_1.20e11ppb_bl_1.20ns.pkl

