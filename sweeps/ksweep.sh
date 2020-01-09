#!/bin/bash
ZGN_N=2
ZGN_t1=1e4
ZGN_Npoints=10
ZGN_noisemax=5.0;
ZGN_C=0.95
filebase=data/knoise

mkdir -p $filebase
rm $filebase/*

for ZGN_temp in `seq 0 $ZGN_Npoints`; do
ZGN_delta=`bc -l <<< "$ZGN_noisemax/$ZGN_Npoints*$ZGN_temp"`
./kuramoto -C $ZGN_C -i $ZGN_delta -y 1 -c 0  -t $ZGN_t1 $filebase/uncorrelatedadd1 &
./kuramoto -C $ZGN_C -i $ZGN_delta -y 1 -c 1  -t $ZGN_t1 $filebase/correlatedadd1 &
done
wait

for ZGN_temp in `seq 0 $ZGN_Npoints`; do
ZGN_delta=`bc -l <<< "$ZGN_noisemax/$ZGN_Npoints*$ZGN_temp"`
./kuramoto -C $ZGN_C -i $ZGN_delta -y 2 -c 0  -t $ZGN_t1 $filebase/uncorrelatedadd2 &
./kuramoto -C $ZGN_C -i $ZGN_delta -y 2 -c 1  -t $ZGN_t1 $filebase/correlatedadd2 &
done
wait

for ZGN_temp in `seq 0 $ZGN_Npoints`; do
ZGN_delta=`bc -l <<< "$ZGN_noisemax/$ZGN_Npoints*$ZGN_temp"`
./kuramoto -C $ZGN_C -i $ZGN_delta -f 1,0 -y 0 -c 0  -t $ZGN_t1 $filebase/uncorrelatedmult &
./kuramoto -C $ZGN_C -i $ZGN_delta -f 1,0 -y 0 -c 1  -t $ZGN_t1 $filebase/correlatedmult &
done
wait
