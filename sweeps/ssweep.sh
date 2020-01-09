#!/bin/bash
ZGN_N=100
ZGN_t1=1e4
ZGN_t2=0.5e3
ZGN_t3=1e3
ZGN_C=0.95
ZGN_delta1=0.0
ZGN_delta2=5.0
ZGN_dt=1e-2
ZGN_dt2=1e-4
ZGN_output=0

ZGN_Npoints=10

for ZGN_temp in `seq 0 $ZGN_Npoints`; do

jobs=`ps aux | grep stuartlandau | wc -l`
while [ $jobs -ge 2 ]; do
sleep 5
jobs=`ps aux | grep stuartlandau | wc -l`
done

ZGN_delta=`bc -l <<< "$ZGN_delta1+(${ZGN_delta2} - ${ZGN_delta1})/${ZGN_Npoints}*${ZGN_temp}"`

echo ${i} ${ZGN_temp} ${ZGN_seed}

./stuartlandau -i $sigma -t 1e3 -o 0.9e3 -s $seed -y 0 data/noise2/uncorrelatedmult &
./stuartlandau -i $sigma -t 1e3 -o 0.9e3 -s $seed -y 0 -c data/noise2/correlatedmult &

./stuartlandau 2 1.0,0.0 0.0,0.0 $ZGN_C $ZGN_delta 0 1 1 $ZGN_t1 $ZGN_t2 $ZGN_t3 $ZGN_dt $ZGN_dt2 $ZGN_seed rkf45 data/noise/correlatedmult1 $ZGN_output &
./stuartlandau 2 1.0,0.0 0.0,0.0 $ZGN_C $ZGN_delta 1 0 1 $ZGN_t1 $ZGN_t2 $ZGN_t3 $ZGN_dt $ZGN_dt2 $ZGN_seed rkf45 data/noise/uncorrelatedadd1 $ZGN_output &
./stuartlandau 2 1.0,0.0 0.0,0.0 $ZGN_C $ZGN_delta 1 1 1 $ZGN_t1 $ZGN_t2 $ZGN_t3 $ZGN_dt $ZGN_dt2 $ZGN_seed rkf45 data/noise/correlatedadd1 $ZGN_output &
./stuartlandau 2 1.0,0.0 0.0,0.0 $ZGN_C $ZGN_delta 2 0 1 $ZGN_t1 $ZGN_t2 $ZGN_t3 $ZGN_dt $ZGN_dt2 $ZGN_seed rkf45 data/noise/uncorrelatedadd2 $ZGN_output &
./stuartlandau 2 1.0,0.0 0.0,0.0 $ZGN_C $ZGN_delta 2 1 1 $ZGN_t1 $ZGN_t2 $ZGN_t3 $ZGN_dt $ZGN_dt2 $ZGN_seed rkf45 data/noise/correlatedadd2 $ZGN_output &
./stuartlandau 2 1.0,0.0 0.0,0.0 $ZGN_C $ZGN_delta 1 0 2 $ZGN_t1 $ZGN_t2 $ZGN_t3 $ZGN_dt $ZGN_dt2 $ZGN_seed rkf45 data/noise/uncorrelatedadd3 $ZGN_output &
./stuartlandau 2 1.0,0.0 0.0,0.0 $ZGN_C $ZGN_delta 1 1 2 $ZGN_t1 $ZGN_t2 $ZGN_t3 $ZGN_dt $ZGN_dt2 $ZGN_seed rkf45 data/noise/correlatedadd3 $ZGN_output &

done
