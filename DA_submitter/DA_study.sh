#!/bin/bash

ptau_max=6.6e-4

qx=62.270
qy=60.295
qprime=20
IMO=40
VRF=8
intensity=1.20

seyMB=1.30
seyMQ=1.30
MB=--MB
MQ=--MQ
[[ "$MB" == "--MB" ]] && MB_name=_${seyMB}seyMB || MB_name=;
[[ "$MQ" == "--MQ" ]] && MQ_name=_${seyMQ}seyMQ || MQ_name=;

job_name=DA_LHC_450GeV${MB_name}${MQ_name}_${intensity}e11ppb_${qx}qx_${qy}qy_${qprime}qprime_${IMO}IMO_${VRF}VRF_${ptau_max}ptau

echo submit_da_ecloud_v1.0.sh --job_name=${job_name} --ptau_max=${ptau_max} \
                     --qx=${qx} --qy=${qy} --qprime=${qprime} --IMO=${IMO} --VRF=${VRF} \
                     --seyMB=${seyMB} --seyMQ=${seyMQ} ${MB} ${MQ} --intensity=${intensity}

./submit_da_ecloud_v1.0.sh --job_name=${job_name} --ptau_max=${ptau_max} \
                     --qx=${qx} --qy=${qy} --qprime=${qprime} --IMO=${IMO} --VRF=${VRF} \
                     --seyMB=${seyMB} --seyMQ=${seyMQ} ${MB} ${MQ} --intensity=${intensity}


