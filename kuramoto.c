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


struct parameters{int N; double *omega; double C; double *noise; double noiseintensity; int noisetype; };

int func (double t, const double y[], double f[], void *params) {

  int N = ((struct parameters *)params)->N;
  int noisetype = ((struct parameters *)params)->noisetype;
  double* omega = ((struct parameters *)params)->omega;
  double C = ((struct parameters *)params)->C;
  double *noise = ((struct parameters *)params)->noise;
  ((struct parameters *)params)->noiseintensity = 0;
  int j, i;

  double cossum=0, sinsum=0, R=0, Theta=0;
  for(j=0; j<N; j++) {
    cossum += 1.0/N*cos(y[j]);
    sinsum += 1.0/N*sin(y[j]);
  }
  R = sqrt(cossum*cossum+sinsum*sinsum);
  Theta = atan2(sinsum, cossum);

  for(j = 0; j<N; j++) {
    //intrinsic dynamics and coupling
    f[j] = omega[j]+C*R*sin(Theta-y[j]);

    //noise
    if(noisetype == 0){
      ((struct parameters *)params)->noiseintensity += (pow(f[j]*(noise[j]-1.0),2.0))/N;
      f[j] += f[j]*(noise[j]-1.0);
    }
    else if(noisetype == 1) {
      ((struct parameters *)params)->noiseintensity += (pow(noise[j],2.0) )/N;
      f[j] += noise[j];
    }
    else{
      ((struct parameters *)params)->noiseintensity += pow(cos(y[j])*noise[j],2.0)/N;
      f[j] += cos(y[j])*noise[j];
    }
  }

  return GSL_SUCCESS;
}

int main (int argc, char* argv[]) {
  struct timeval start,end;
  int i,j,k;

  int N=2;
  double *omega;
  omega = calloc(N, sizeof(double));
  char string[1048576];
  strcpy(string, "-0.5,0.5");
  char *split = strtok(string, ",");
  for(j=0; j<N; j++) {
    omega[j] = atof(split);
    split=strtok(NULL,",");
  }

  double C = 0.95;
  double sigma = 1.0;
  int noisetype = 1;
  int correlated = 0;
  double t = 0.;
  double tmax = 1e3;
  double ta = 0;
  double to = 5e2;
  double dt = 1e-2;
  double dt2 = 1e-3;
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
        {"frequencies",  required_argument, 0, 'f'},
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

    opt = getopt_long (argc, argv, "n:a:cC:d:D:f:hi:o:s:t:vy:", long_options, &option_index);
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
      free(omega);
      omega = calloc(N, sizeof(double));
      strcpy(string, optarg);
      char *split = strtok(string, ",");
      for(j=0; j<N; j++) {
        omega[j] = atof(split);
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
  double ti = dt;
  double order = 0;
  double netnoiseintensity = 0;
  double *y, *yerr, *noise, *phase;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  FILE *out, *outsignal, *outphase, *outnoise;
  char file[256];

  y=calloc(N,sizeof(double));
  yerr=calloc(N,sizeof(double));
  noise=calloc(N,sizeof(double));
  phase=calloc(N,sizeof(double));
  struct parameters params={N, omega, C, noise, 0, noisetype};

  //create the random noise function.
  gsl_rng_set(r,seed);
  for(j=0; j<N; j++) {
    double theta=2*3.14*gsl_rng_uniform(r);
    y[2*j] = theta;
  }

  //Output the noise data
  if(verbose){
    strcpy(file,filebase);
    strcat(file, "noise.dat");
    outnoise=fopen(file,"w");

    strcpy(file,filebase);
    strcat(file, ".out");
    out=fopen(file,"w");

    strcpy(file,filebase);
    strcat(file, "phase.dat");
    outphase=fopen(file,"w");

    strcpy(file,filebase);
    strcat(file, ".dat");
    outsignal=fopen(file,"w");

    fprintf(out, "%i %f %f %f %f\n", N, tmax, dt, sigma, C);

    fprintf(out, "\n");
  }
  fflush(stdout);

  //Set up integrator
  gsl_odeiv2_system sys = {func, NULL, N, &params};
  gsl_odeiv2_step * step;
  step = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkf45, N);

  gettimeofday(&start,NULL);

  //Do integration
  int count=0;
  double maxerr;
  while(count < Nt) {
    count++;
    maxerr=0;

    if (noisetype==0) {
      if(sigma == 0.0) {
        for(j=0; j<N; j++) {
          noise[j] = 1.0;
        }
      }
      else{
        if(correlated==0) {
          for(j=0; j<N; j++)
          noise[j] = gsl_ran_gamma(r,2*dt/(sigma*sigma), sigma*sigma/(2*dt));
        }
        else {
          double ran = gsl_ran_gamma(r,2*dt/(sigma*sigma), sigma*sigma/(2*dt));
          for(j=0; j<N; j++)
          noise[j] = ran;
        }
      }
    }
    else {
      if(correlated==0) {
        for(j=0; j<N; j++)
        noise[j] = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
      }
      else {
        double ran = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
        for(j=0; j<N; j++)
        noise[j] = ran;
      }
    }

    for(k=0; k<dt/dt2; k++) {
      t=count*dt+k*dt2;

      gsl_odeiv2_step_apply (step, count*dt+k*dt2, dt2, y, yerr, NULL, NULL, &sys);
      for(j=0; j<N; j++) {
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
      for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
          order += cos(y[j]-y[k])/(N*N);
        }
      }
    }

    if(verbose){
      if(t>ta){
        fwrite(phase, sizeof(double), N, outphase);
        fflush(outphase);
        fwrite(noise,sizeof(double), N, outnoise);
        fflush(outnoise);
        fwrite(y, sizeof(double), N, outsignal);
        fflush(outsignal);
      }
      gettimeofday(&end,NULL);
      printf("%6f\t%6f\t%6f\t%6e\t\r",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fprintf(out, "%f\t%f\t%f\t%6e\t\n",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fflush(stdout);
    }

  }

  gettimeofday(&end,NULL);
  printf("runtime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));

  //Output results

  if(verbose){
    fflush(outphase);
    fflush(outsignal);
    fflush(outnoise);
    fflush(out);
    fclose(outphase);
    fclose(outsignal);
    fclose(outnoise);
    fclose(out);

    printf("%f %f \n", order/(Nt-Nto), sqrt(netnoiseintensity*2*dt/Nt));
  }

  strcpy(file,filebase);
  strcat(file, "meanphase.dat");
  out=fopen(file,"a+");
  fprintf(out, "%i %f %f %i %f %f\n", N, C, sigma, seed, order/(Nt-Nto), sqrt(netnoiseintensity*2*dt/Nt));
  fflush(out);
  fclose(out);

  gsl_odeiv2_step_free(step);
  gsl_rng_free(r);
  free(omega);
  free(y);
  free(yerr);
  free(noise);
  free(phase);
  return 0;
}
