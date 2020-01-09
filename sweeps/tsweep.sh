#! /bin/bash
N=100
count=100
jobs=100
t1=6e2
t2=1e2

for seed in {1..10}; do
echo $seed
mkdir -p data/tnoise/$seed

for i in `seq 0 $count`; do
K=`bc -l <<< "5.0/${count}*${i}"`

js=`jobs | wc -l`
while [ $js -ge $jobs ]; do
sleep 1
js=`jobs | wc -l`
done

./twolayer -i 0 -t $t1 -o $t2 -a $t2 -s $seed -K $K data/tnoise/${seed}/noiseless$i  &
./twolayer -t $t1 -o $t2 -a $t2 -s $seed -K $K data/tnoise/${seed}/uncorrelated$i  &
./twolayer -c -t $t1 -o $t2 -a $t2 -s $seed -K $K data/tnoise/${seed}/correlated$i  &

done
wait
done
done
