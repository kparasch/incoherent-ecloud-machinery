#!/bin/bash

qx=62.270
qy=60.295
qprime=15
IMO=40
VRF=6

#./do_all.sh ${qx} ${qy} ${qprime} ${IMO} ${VRF}    

for qprime in 0 3 5
do
    ./do_all.sh ${qx} ${qy} ${qprime} 0 ${VRF}    
    ./do_all.sh ${qx} ${qy} ${qprime} 10 ${VRF}    
    ./do_all.sh ${qx} ${qy} ${qprime} 20 ${VRF}    
    ./do_all.sh ${qx} ${qy} ${qprime} 40 ${VRF}    
done
