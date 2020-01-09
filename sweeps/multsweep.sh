
for n in `seq 0 10`; do
sigma=`bc -l <<< "0.5*$n"`
for seed in `seq 1 seeds`; do
./stuartlandau_2 2 0,1 0,0 0.95 $sigma 0 0 1e3 0.9e3 1e3 1e-2 1e-4 $seed rkf45 data/noise2/uncorrelatedmult 0 &
./stuartlandau_2 2 0,1 0,0 0.95 $sigma 0 1 1e3 0.9e3 1e3 1e-2 1e-4 $seed rkf45 data/noise2/correlatedmult 0 &
done
done
