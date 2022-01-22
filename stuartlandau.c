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

struct parameters{int N; double *alpha; double* gamma; double C; double *noise; double *noiseout; double noiseintensity; int noisetype; };

int func (double t, const double y[], double f[], void *params) {

  int N = ((struct parameters *)params)->N;
  int noisetype = ((struct parameters *)params)->noisetype;
  double* alpha = ((struct parameters *)params)->alpha;
  double* gamma = ((struct parameters *)params)->gamma;
  double C = ((struct parameters *)params)->C;
  double *noise = ((struct parameters *)params)->noise;
  double *noiseout = ((struct parameters *)params)->noiseout;
  ((struct parameters *)params)->noiseintensity = 0;
  int j, i;

  for(j = 0; j<N; j++) {
    //intrinsic dynamics
    f[2*j] = y[2*j]-alpha[j]*y[2*j+1]-(y[2*j]*y[2*j]+y[2*j+1]*y[2*j+1])*(y[2*j]-gamma[j]*y[2*j+1]);
    f[2*j+1] = y[2*j+1]+alpha[j]*y[2*j]-(y[2*j]*y[2*j]+y[2*j+1]*y[2*j+1])*(y[2*j+1]+gamma[j]*y[2*j]);

    //coupling
    for(i = 0; i<N; i++){
      f[2*j] -= C/N*(y[2*j]*y[2*i+1]-y[2*i]*y[2*j+1])*y[2*j+1];
      f[2*j+1] += C/N*(y[2*j]*y[2*i+1]-y[2*i]*y[2*j+1])*y[2*j];
    }

    //noise
    if(noisetype == 0){
      ((struct parameters *)params)->noiseintensity += pow((noise[j]-1.0),2.0)/N;
      noiseout[2*j] = f[2*j]*(noise[j]-1);
      noiseout[2*j+1] = f[2*j+1]*(noise[j]-1);
      // f[2*j] *= noise[j];
      // f[2*j+1] *= noise[j];
    }
    else if(noisetype == 1) {
      ((struct parameters *)params)->noiseintensity += pow(noise[j],2.0)/N;
      noiseout[2*j] = -y[2*j+1]*noise[j];
      noiseout[2*j+1] = y[2*j]*noise[j];
      // f[2*j] += -y[2*j+1]*noise[j];
      // f[2*j+1] += y[2*j]*noise[j];
    }
    else{
      ((struct parameters *)params)->noiseintensity += pow(noise[j],2.0)/N;
      noiseout[2*j] = 0;
      noiseout[2*j+1] = noise[j];
      // f[2*j+1] += noise[j];
    }
    f[2*j] += noiseout[2*j];
    f[2*j+1] += noiseout[2*j+1];

  }
  return GSL_SUCCESS;
}

int main (int argc, char* argv[]) {
  struct timeval start,end;
  int i,j,k;

  int N=2;
  double *alpha, *gamma;
  alpha = calloc(N, sizeof(double));
  gamma = calloc(N, sizeof(double));
  char string[1048576];
  strcpy(string, "1,0");
  char *split = strtok(string, ",");
  for(j=0; j<N; j++) {
    alpha[j] = atof(split);
    split=strtok(NULL,",");
  }
  strcpy(string, "0,0");
  split = strtok(string, ",");
  for(j=0; j<N; j++) {
    gamma[j] = atof(split);
    split=strtok(NULL,",");
  }

  double C = 0.95;
  double sigma = 1.0;
  int noisetype = 1;
  int correlated = 0;
  double t = 0.;
  double tmax = 1e3;
  double ta = 5e2;
  double to = 5e2;
  double dt = 1e-1;
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
        {"coupling",  required_argument, 0, 'C'},
        {"common",  no_argument, 0, 'c'},
        {"dt",  required_argument, 0, 'd'},
        {"noisestep",  required_argument, 0, 'D'},
        {"alpha",  required_argument, 0, 'f'},
        {"gamma",  required_argument, 0, 'g'},
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

    opt = getopt_long (argc, argv, "n:a:cC:d:D:f:g:hi:o:s:t:vy:", long_options, &option_index);
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
    case 'f':
      free(alpha);
      alpha = calloc(N, sizeof(double));
      strcpy(string, optarg);
      split = strtok(string, ",");
      for(j=0; j<N; j++) {
        alpha[j] = atof(split);
        split=strtok(NULL,",");
      }
      break;
    case 'g':
      free(gamma);
      gamma = calloc(N, sizeof(double));
      strcpy(string, optarg);
      split = strtok(string, ",");
      for(j=0; j<N; j++) {
        gamma[j] = atof(split);
        split=strtok(NULL,",");
      }
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
  double order = 0;
  double netnoiseintensity = 0;
  double *y, *yerr, *noise, *noiseout, *phase;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  FILE *out, *outsignal, *outnoise, *outnoise2;
  char file[1048576];

  y=calloc(2*N,sizeof(double));
  noiseout=calloc(2*N,sizeof(double));
  yerr=calloc(2*N,sizeof(double));
  noise=calloc(N,sizeof(double));
  phase=calloc(N,sizeof(double));
  struct parameters params={N, alpha, gamma, C, noise,noiseout, 0, noisetype};

  //create the random noise function.
  gsl_rng_set(r,seed);
  for(j=0; j<N; j++) {
    double theta = 4*asin(1.0)*gsl_rng_uniform(r);
    y[2*j] = cos(theta);
    y[2*j+1] = sin(theta);
  }

  //Open output files
  strcpy(file,filebase);
  strcat(file, ".dat");
  outsignal=fopen(file,"w");

  strcpy(file,filebase);
  strcat(file, "noise.dat");
  outnoise=fopen(file,"w");

  strcpy(file,filebase);
  strcat(file, "noise2.dat");
  outnoise2=fopen(file,"w");

  strcpy(file,filebase);
  strcat(file, ".out");
  out=fopen(file,"w");
  fprintf(out, "%i %f %f %f %f\n", N, tmax-ta, dt, sigma, C);
  for(int i=0; i<argc; i++){
    fprintf(out,"%s ",argv[i]);
  }
  fprintf(out,"\n");
  //Set up integrator
  gsl_odeiv2_system sys = {func, NULL, 2*N, &params};
  gsl_odeiv2_step * step;
  step = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkf45, 2*N);

  gettimeofday(&start,NULL);


  //Do integration
  for (j=0; j<N; j++)
  phase[j] = atan2(y[2*j+1], y[2*j]);

  int count=0;
  double maxerr, ran1, ran2;
  while(count < Nt) {
    count++;
    maxerr=0;

    if(noisetype == 0) {
      if(sigma == 0.0) {
        for (j=0; j<N; j++) {
          noise[j] = 1.0;
        }
      }
      else if (correlated == 0) {
        for (j=0; j<N; j++) {
          ran1 = gsl_ran_gamma(r,dt/(2*sigma*sigma), 2*sigma*sigma/dt);
          noise[j] = ran1;
        }
      }
      else if (correlated == 1){
        ran1 = gsl_ran_gamma(r,dt/(2*sigma*sigma), 2*sigma*sigma/dt);
        for (j=0; j<N; j++) {
          noise[j] = ran1;
        }
      }
    }
    else{
      if (correlated == 0) {
        for (j=0; j<N; j++) {
          ran1 = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
          noise[j] = ran1;
        }

      }
      else if (correlated == 1){
        ran1 = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
        for (j=0; j<N; j++) {
          noise[j] = ran1;
        }
      }
    }

    for(k=0; k<dt/dt2; k++) {
      t=count*dt+k*dt2;

      gsl_odeiv2_step_apply (step, count*dt+k*dt2, dt2, y, yerr, NULL, NULL, &sys);
      for(j=0; j<2*N; j++) {
        if(fabs(yerr[j]) > maxerr)
        maxerr = fabs(yerr[j]);
        if(isnan(y[j])) {
          printf("\nerror nan at %6f\n", t/tmax);
          return 1;
        }
      }
    }
    netnoiseintensity += params.noiseintensity;

    for (j=0; j<N; j++)
    phase[j] = atan2(y[2*j+1], y[2*j]);

    if (t>=to) {
      for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
          order += cos(phase[j]-phase[k])/(N*N);
        }
      }
    }

    if(t > ta) {
      fwrite(y, sizeof(double), 2*N, outsignal);
      fwrite(noiseout, sizeof(double), 2*N, outnoise);
      fwrite(noise, sizeof(double), N, outnoise2);
      fflush(outsignal);
      fflush(outnoise);
      fflush(outnoise2);
    }
    if(verbose){
      gettimeofday(&end,NULL);
      printf("%6f\t%6f\t%6f\t%6e\t\r",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fprintf(out, "%f\t%f\t%f\t%6e\t\n",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fflush(stdout);
    }

  }

  gettimeofday(&end,NULL);
  printf("\nruntime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));
  printf("%f %f \n", order/(Nt-Nto), sqrt(netnoiseintensity*2*dt/Nt));

  //Output results
  fflush(outsignal);
  fclose(outsignal);
  fflush(outnoise);
  fclose(outnoise);
  fflush(outnoise2);
  fclose(outnoise2);
  fprintf(out, "\nruntime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));
  fprintf(out, "%i %f %f %i %f %f\n", N, C, sigma, seed, order/(Nt-Nto), sqrt(netnoiseintensity*2*dt/Nt));
  fflush(out);
  fclose(out);


  free(alpha);
  free(gamma);
  free(y);
  free(yerr);
  free(noise);
  free(noiseout);
  free(phase);
  gsl_odeiv2_step_free(step);
  gsl_rng_free(r);

  return 0;
}
