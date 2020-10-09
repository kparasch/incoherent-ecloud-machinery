#!/bin/bash

ptau_max=0.e-4

qx=62.270
qy=60.295
qprime=20
IMO=40
VRF=8
intensity=0.00

seyMB=1.30
seyMQ=1.30
#MB=--MB ; [[ "$MB" == "--MB" ]] && MB_name=_${seyMB}seyMB || MB_name=;
#MQ=--MQ ; [[ "$MQ" == "--MQ" ]] && MQ_name=_${seyMQ}seyMQ || MQ_name=;
#
#job_name=FMA_LHC_450GeV${MB_name}${MQ_name}_${intensity}e11ppb_${qx}qx_${qy}qy_${qprime}qprime_${IMO}IMO_${VRF}VRF_${ptau_max}ptau
#./submit_fma_ecloud_v1.0.sh --job_name=${job_name} --ptau_max=${ptau_max} \
#                     --qx=${qx} --qy=${qy} --qprime=${qprime} --IMO=${IMO} --VRF=${VRF} \
#                     --seyMB=${seyMB} --seyMQ=${seyMQ} ${MB} ${MQ} --intensity=${intensity}

for intensity in 1.20 1.10 1.00 0.90 0.80 0.70 0.60 0.50 0.40 0.30 
do

    MB=--MB ; [[ "$MB" == "--MB" ]] && MB_name=_${seyMB}seyMB || MB_name=;
    MQ=--MQ ; [[ "$MQ" == "--MQ" ]] && MQ_name=_${seyMQ}seyMQ || MQ_name=;
    job_name=FMA_LHC_450GeV${MB_name}${MQ_name}_${intensity}e11ppb_${qx}qx_${qy}qy_${qprime}qprime_${IMO}IMO_${VRF}VRF_${ptau_max}ptau
    ./submit_fma_ecloud_v1.0.sh --job_name=${job_name} --ptau_max=${ptau_max} \
                         --qx=${qx} --qy=${qy} --qprime=${qprime} --IMO=${IMO} --VRF=${VRF} \
                         --seyMB=${seyMB} --seyMQ=${seyMQ} ${MB} ${MQ} --intensity=${intensity}

    MB=--MB ; [[ "$MB" == "--MB" ]] && MB_name=_${seyMB}seyMB || MB_name=;
    MQ=     ; [[ "$MQ" == "--MQ" ]] && MQ_name=_${seyMQ}seyMQ || MQ_name=;
    job_name=FMA_LHC_450GeV${MB_name}${MQ_name}_${intensity}e11ppb_${qx}qx_${qy}qy_${qprime}qprime_${IMO}IMO_${VRF}VRF_${ptau_max}ptau
    ./submit_fma_ecloud_v1.0.sh --job_name=${job_name} --ptau_max=${ptau_max} \
                         --qx=${qx} --qy=${qy} --qprime=${qprime} --IMO=${IMO} --VRF=${VRF} \
                         --seyMB=${seyMB} --seyMQ=${seyMQ} ${MB} ${MQ} --intensity=${intensity}

    MB=     ; [[ "$MB" == "--MB" ]] && MB_name=_${seyMB}seyMB || MB_name=;
    MQ=--MQ ; [[ "$MQ" == "--MQ" ]] && MQ_name=_${seyMQ}seyMQ || MQ_name=;
    job_name=FMA_LHC_450GeV${MB_name}${MQ_name}_${intensity}e11ppb_${qx}qx_${qy}qy_${qprime}qprime_${IMO}IMO_${VRF}VRF_${ptau_max}ptau
    ./submit_fma_ecloud_v1.0.sh --job_name=${job_name} --ptau_max=${ptau_max} \
                         --qx=${qx} --qy=${qy} --qprime=${qprime} --IMO=${IMO} --VRF=${VRF} \
                         --seyMB=${seyMB} --seyMQ=${seyMQ} ${MB} ${MQ} --intensity=${intensity}
done

#for qx in 62.260 62.265 62.270 62.275 62.280 62.285 62.290 62.295 62.300 62.305 62.310 62.310 62.315 62.320 62.325 62.330 62.335
#do
#    for qy in 60.260 60.265 60.270 60.275 60.280 60.285 60.290 60.295 60.300 60.305 60.310 60.310 60.315 60.320 60.325 60.330 60.335
#    do
#        intensity=1.20
#    MB=--MB ; [[ "$MB" == "--MB" ]] && MB_name=_${seyMB}seyMB || MB_name=;
#    MQ=--MQ ; [[ "$MQ" == "--MQ" ]] && MQ_name=_${seyMQ}seyMQ || MQ_name=;
#    job_name=FMA_LHC_450GeV${MB_name}${MQ_name}_${intensity}e11ppb_${qx}qx_${qy}qy_${qprime}qprime_${IMO}IMO_${VRF}VRF_${ptau_max}ptau
#    ./submit_fma_ecloud_v1.0.sh --job_name=${job_name} --ptau_max=${ptau_max} \
#                         --qx=${qx} --qy=${qy} --qprime=${qprime} --IMO=${IMO} --VRF=${VRF} \
#                         --seyMB=${seyMB} --seyMQ=${seyMQ} ${MB} ${MQ} --intensity=${intensity}
#
#    done
#done
