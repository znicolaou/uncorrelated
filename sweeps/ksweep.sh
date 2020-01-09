#!/bin/bash
ZGN_N=2
ZGN_t1=1e4
ZGN_t2=1e3
ZGN_dt=1e-4
ZGN_dt2=1e-5
ZGN_Npoints=10
ZGN_noisemax=5.0;
ZGN_C=0.95
filebase=data/knoise

mkdir -p $filebase
rm $filebase/*

for ZGN_temp in `seq 0 $ZGN_Npoints`; do
ZGN_delta=`bc -l <<< "$ZGN_noisemax/$ZGN_Npoints*$ZGN_temp"`
./kuramoto -C $ZGN_C -i $ZGN_delta -y 1  -t $ZGN_t1 -o $ZGN_t2 -D $ZGN_dt -d $ZGN_dt2 $filebase/uncorrelatedadd1 &
./kuramoto -C $ZGN_C -i $ZGN_delta -y 1 -c -t $ZGN_t1  -o $ZGN_t2 -D $ZGN_dt -d $ZGN_dt2 $filebase/correlatedadd1 &
done
wait

for ZGN_temp in `seq 0 $ZGN_Npoints`; do
ZGN_delta=`bc -l <<< "$ZGN_noisemax/$ZGN_Npoints*$ZGN_temp"`
./kuramoto -C $ZGN_C -i $ZGN_delta -y 2 -t $ZGN_t1  -o $ZGN_t2 -D $ZGN_dt -d $ZGN_dt2 $filebase/uncorrelatedadd2 &
./kuramoto -C $ZGN_C -i $ZGN_delta -y 2 -c -t $ZGN_t1  -o $ZGN_t2 -D $ZGN_dt -d $ZGN_dt2 $filebase/correlatedadd2 &
done
wait

for ZGN_temp in `seq 0 $ZGN_Npoints`; do
ZGN_delta=`bc -l <<< "$ZGN_noisemax/$ZGN_Npoints*$ZGN_temp"`
./kuramoto -C $ZGN_C -i $ZGN_delta -f 1,0 -y 0 -t $ZGN_t1  -o $ZGN_t2 -D $ZGN_dt -d $ZGN_dt2 $filebase/uncorrelatedmult &
./kuramoto -C $ZGN_C -i $ZGN_delta -f 1,0 -y 0 -c -t $ZGN_t1  -o $ZGN_t2 -D $ZGN_dt -d $ZGN_dt2 $filebase/correlatedmult &
done
wait
