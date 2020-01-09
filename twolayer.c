//test
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <getopt.h>


struct parameters{int N; double *omega; double C; double K; double *noise; double noiseintensity; int noisetype; };

int func (double t, const double y[], double f[], void *params) {

  int N = ((struct parameters *)params)->N;
  int noisetype = ((struct parameters *)params)->noisetype;
  double* omega = ((struct parameters *)params)->omega;
  double C = ((struct parameters *)params)->C;
  double K = ((struct parameters *)params)->K;
  double *noise = ((struct parameters *)params)->noise;
  ((struct parameters *)params)->noiseintensity = 0;
  int j, i;

  double cossum=0, sinsum=0, R1=0, Theta1=0, R2=0, Theta2=0;
  for(j=0; j<2; j++) {
    cossum += 1.0/2*cos(y[j]);
    sinsum += 1.0/2*sin(y[j]);
  }
  R1 = sqrt(cossum*cossum+sinsum*sinsum);
  Theta1 = atan2(sinsum, cossum);

  cossum=0;
  sinsum=0;
  for(j=2; j<2+N; j++) {
    cossum += 1.0/N*cos(y[j]);
    sinsum += 1.0/N*sin(y[j]);
  }
  R2 = sqrt(cossum*cossum+sinsum*sinsum);
  Theta2 = atan2(sinsum, cossum);
  for(j = 0; j<2; j++) {
    f[j] = omega[j]+C*R1*sin(Theta1-y[j]);
    if(noisetype == 1){
      ((struct parameters *)params)->noiseintensity += (pow(noise[j],2.0) )/N;
      f[j] += noise[j];
    }
    else if(noisetype == 2) {
      ((struct parameters *)params)->noiseintensity += pow(cos(y[j])*noise[j],2.0)/N;
      f[j] += cos(y[j])*noise[j];
    }
    else{
      ((struct parameters *)params)->noiseintensity += (pow(f[j]*(noise[j]-1.0),2.0))/N;
      f[j] += f[j]*(noise[j]-1.0);
    }
  }

  for(j = 2; j<2+N; j++) {
    f[j] = omega[j]+K*R2*R1*R1*sin(Theta2-y[j]);
  }
  //filebases=(noiseless uncorrelated correlated); noises=(0.0 1.0 1.0); correlateds=(0 0 1); for j in {0..2}; do for i in {0..10}; do K1=`bc -l <<< "10.0/10*$i"`; ./noiseslave 100 0.95 $K1 0.0 ${noises[j]} ${correlateds[j]} 1e3 1e2 1e-2 1e-3 1 rkf45 data/noise/${filebases[j]}$i 1 &> $i.out & done; wait; done;


  return GSL_SUCCESS;
}

int main (int argc, char* argv[]) {
  struct timeval start,end;
  int i,j,k;

  int N=100;
  double C = 0.95;
  double K = 5.0;
  double sigma = 1.0;
  int correlated = 0;
  int noisetype = 1;
  double t = 0.;
  double tmax = 1e3;
  double ta = 5e2;
  double to = 5e2;
  double dt = 1e-2;
  double dt2 = 1e-4;
  int seed = 1;
  char* filebase = NULL;
  int verbose = 0;
  int opt;

  while (1)
  {
    static struct option long_options[] =
    {
      {"animate",  required_argument, 0, 'a'},
      {"coupling1",  required_argument, 0, 'C'},
      {"coupling2",  required_argument, 0, 'K'},
      {"common",  no_argument, 0, 'c'},
      {"dt",  required_argument, 0, 'd'},
      {"noisestep",  required_argument, 0, 'D'},
      {"help",  no_argument, 0, 'h'},
      {"intensity",  required_argument, 0, 'i'},
      {"num",  required_argument, 0, 'n'},
      {"order",  required_argument, 0, 'o'},
      {"seed",  required_argument, 0, 's'},
      {"time",  required_argument, 0, 't'},
      {"verbose",  no_argument, 0, 'v'},
      {"type",  required_argument, 0, 'y'},
      {0, 0, 0, 0}
    };
    int option_index = 0;

    opt = getopt_long (argc, argv, "n:a:cC:d:D:hi:K:o:s:t:vy:", long_options, &option_index);
    if (opt == -1)
    break;
    switch (opt)
    {
      case 'n':
      N=atoi(optarg);
      break;
      case 'a':
      ta=atof(optarg);
      break;
      case 'C':
      C=atof(optarg);
      break;
      case 'c':
      correlated=1;
      break;
      case 'd':
      dt2=atof(optarg);
      break;
      case 'D':
      dt=atof(optarg);
      break;

      case 'h':
      printf("usage: ./kuramoto [-n NUM] [-c K] [-f FREQUENCIES] [-C] [-t TIME] [-a ATIME] [-o OTIME] [-d DT] [-D DT2] [-i SIGMA] [-s SEED] [-y TYPE][-v] filebase \n\n");
      printf("optional arguments:\n");
      printf("\t -a ATIME, --animate ATIME \t\t Time to start outputting animation data. Default 1e1. \n");
      printf("\t -C, --common  \t\t Use common noise.\n");
      printf("\t -c C, --coupling C \t\t Coupling constant. Default 0.95.\n");
      printf("\t -d DT, --dt DT \t\t Integration timestep. Default 1e-3.\n");
      printf("\t -D DT2, --noisestep DT2 \t\t Noise sampling rate. Default 1e-2.\n");
      printf("\t -i SIGMA, --intensity SIGMA \t\t Noise intensity. Default 1.\n");
      printf("\t -n NUM, --num NUM \t\t Number of grid points in each dimension. Default 1536.\n");
      printf("\t -o OTIME, --order OTIME \t\t Time to begin averaging order parameter. Default 5e2.\n");
      printf("\t -s SEED, --seed SEED \t\t Random seed. Default 1.\n");
      printf("\t -t TIME, --time TIME \t\t Total integration time. Default 1e1.\n");
      printf("\t -v, --verbose \t\t Verbose output.\n");
      printf("\t -y, --type \t\t Noise type. 0 for multiplicative gamma, 1 for additive constant phase sensitivity, 2 for additive trigonometric sensitivity. Default 1.\n");

      printf("positional arguments:\n");

      printf("filebase \t\t Base file name for output. filebaseout.dat contains time step data, outlast.dat contains the last state, outanimation.dat contains the states after ta. \n");

      exit(0);
      case 'i':
      sigma=atof(optarg);
      break;
      case 'K':
      K=atof(optarg);
      break;
      case 'o':
      to=atoi(optarg);
      break;
      case 's':
      seed=atoi(optarg);
      break;
      case 't':
      tmax=atof(optarg);
      break;
      case 'v':
      verbose=1;
      break;
      case 'y':
      noisetype=atoi(optarg);
      break;
      case '?':
      abort ();
      default:
      abort ();
    }
  }

  if (optind == argc-1) {
    filebase=argv[optind];
  }
  else{
    printf("Specify the outpuf filebase, or use --help for usage. \n");
    return 0;
  }

  int Nt = (int)tmax/dt;
  int Nto = (int)to/dt;
  double order1 = 0, order2 = 0;
  double netnoiseintensity = 0;

  double *omega, *y, *yerr, *noise;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  FILE *out, *outsignal, *outfrequencies;
  char file[256];

  omega = calloc(2+N, sizeof(double));
  y=calloc(2+N,sizeof(double));
  yerr=calloc(2+N,sizeof(double));
  noise=calloc(2,sizeof(double));
  struct parameters params={N, omega, C, K, noise, 0, noisetype};

  //random ics
  gsl_rng_set(r,seed);
  for(j=0; j<2+N; j++) {
    y[j] = 2*3.14*gsl_rng_uniform(r);
    omega[j] = gsl_ran_cauchy(r,1);
  }
  omega[0]=-0.5;
  omega[1]=0.5;

  //open outputs
  strcpy(file,filebase);
  strcat(file, ".dat");
  outsignal=fopen(file,"w");
  strcpy(file,filebase);
  strcat(file, "frequencies.dat");
  outfrequencies=fopen(file,"w");
  if(verbose){
    strcpy(file,filebase);
    strcat(file, ".out");
    out=fopen(file,"w");

    fprintf(out, "%i %f %f %f %f %f\n", N, tmax-ta, dt, sigma, C, K);

    fprintf(out, "\n");
  }
  fflush(stdout);

  for(j = 0; j<2+N; j++) {
    fprintf(outfrequencies, "%f ", omega[j]);
  }
  fflush(outfrequencies);
  fclose(outfrequencies);
  
  //Set up integrator
  gsl_odeiv2_system sys = {func, NULL, 2+N, &params};
  gsl_odeiv2_step * step;
  step = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkf45, 2+N);

  gettimeofday(&start,NULL);

  //Do integration
  int count=0;
  double maxerr;
  while(count < Nt) {
    count++;
    maxerr=0;

    if (noisetype==0) {
      if(sigma == 0.0) {
        for(j=0; j<2; j++) {
          noise[j] = 1.0;
        }
      }
      else{
        if(correlated==0) {
          for(j=0; j<2; j++)
          noise[j] = gsl_ran_gamma(r,2*dt/(sigma*sigma), sigma*sigma/(2*dt));
        }
        else {
          double ran = gsl_ran_gamma(r,2*dt/(sigma*sigma), sigma*sigma/(2*dt));
          for(j=0; j<2; j++)
          noise[j] = ran;
        }
      }
    }
    else {
      if(correlated==0) {
        for(j=0; j<2; j++)
        noise[j] = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
      }
      else {
        double ran = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
        for(j=0; j<2; j++)
        noise[j] = ran;
      }
    }

    for(k=0; k<dt/dt2; k++) {
      t=count*dt+k*dt2;
      gsl_odeiv2_step_apply (step, count*dt+k*dt2, dt2, y, yerr, NULL, NULL, &sys);

      for(j=0; j<2+N; j++) {
        if(fabs(yerr[j]) > maxerr)
        maxerr = fabs(yerr[j]);
        if(isnan(y[j])) {
          printf("\nerror nan at %6f\n", t/tmax);
          return 1;
        }
      }
    }
    netnoiseintensity += params.noiseintensity;

    if (t>=to) {
      for (j=0; j<2; j++) {
        for (k=0; k<2; k++) {
          order1 += (cos(y[j])*cos(y[k])+sin(y[j])*sin(y[k]))/(2*2);
        }
      }
      for (j=2; j<N+2; j++) {
        for (k=2; k<N+2; k++) {
          order2 += (cos(y[j])*cos(y[k])+sin(y[j])*sin(y[k]))/(N*N);
        }
      }
    }
    if(t>ta){
      fwrite(y, sizeof(double), N+2, outsignal);
      fflush(outsignal);
    }
    if(verbose){
      gettimeofday(&end,NULL);
      printf("%6f\t%6f\t%6f\t%6e\t\r",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fprintf(out, "%f\t%f\t%f\t%6e\t\n",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fflush(stdout);
    }
  }

  gettimeofday(&end,NULL);
  printf("runtime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));

  //Output results
  fflush(outsignal);
  fclose(outsignal);

  if(verbose){
    fflush(out);
    fclose(out);

    printf("%f %f %f \n", order1/(Nt-Nto), order2/(Nt-Nto), sqrt(netnoiseintensity*2*dt/Nt));
  }

  strcpy(file,filebase);
  strcat(file, "meanphase.dat");
  out=fopen(file,"a+");
  fprintf(out, "%i %f %f %f %i %f %f %f\n", N, C, K, sigma, seed, order1/(Nt-Nto), order2/(Nt-Nto), sqrt(netnoiseintensity*2*dt/Nt));
  fflush(out);
  fclose(out);

  gsl_odeiv2_step_free(step);
  gsl_rng_free(r);
  free(omega);
  free(y);
  free(yerr);
  free(noise);
  return 0;
}
