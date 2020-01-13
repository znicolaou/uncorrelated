# Files in the uncorrelated repository
This repository contains three C source files (kuramoto.c, stuartlandau.c, and twolayer.c), which contain code to simulate coupled oscillators, and three Mathematica notebooks (kuramoto.nb, stuartlandau.nb, and twolayer.nb), which plot results, and scripts for running batches.
# System requirements
Compiling requires the gnu scientific library: https://www.gnu.org/software/gsl/.
# Usage
Compile each .c file with gcc, as in
```
gcc -o kuramoto kuramoto.c -lgsl -lgslcblas -lm
```
Running with the `-h` flag gives usage. For example, `./kuramoto -h` outputs
```
usage: ./kuramoto [-n NUM] [-c K] [-f FREQUENCIES] [-C] [-t TIME] [-a ATIME] [-o OTIME] [-d DT] [-D DT2] [-i SIGMA] [-s SEED] [-y TYPE][-v] filebase 

optional arguments:
	 -a ATIME, --animate ATIME 		 Time to start outputting animation data. Default 1e1. 
	 -C, --common  		 Use common noise.
	 -c C, --coupling C 		 Coupling constant. Default 0.95.
	 -d DT, --dt DT 		 Integration timestep. Default 1e-3.
	 -D DT2, --noisestep DT2 		 Noise sampling rate. Default 1e-2.
	 -i SIGMA, --intensity SIGMA 		 Noise intensity. Default 1.
	 -n NUM, --num NUM 		 Number of grid points in each dimension. Default 1536.
	 -o OTIME, --order OTIME 		 Time to begin averaging order parameter. Default 5e2.
	 -s SEED, --seed SEED 		 Random seed. Default 1.
	 -t TIME, --time TIME 		 Total integration time. Default 1e1.
	 -v, --verbose 		 Verbose output.
	 -y, --type 		 Noise type. 0 for multiplicative gamma, 1 for additive constant phase sensitivity, 2 for additive trigonometric sensitivity. Default 1.
positional arguments:
	filebase 		 Base file name for output. filebaseout.dat contains time step data, outlast.dat contains the last state, outanimation.dat contains the states after ta.

```
___
# Output files
Each program produces a filebase.out file, containing information about the command that was run and the results, and a filebase.dat file, containing binary data of all the oscillators states starting at the time specified by the -a flag. 
