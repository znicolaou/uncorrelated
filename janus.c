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


struct parameters{int N; double *omega; double C; double B; double *noise; int noisetype; };

int func (double t, const double y[], double f[], void *params) {

  int N = ((struct parameters *)params)->N;
  int noisetype = ((struct parameters *)params)->noisetype;
  double* omega = ((struct parameters *)params)->omega;
  double B = ((struct parameters *)params)->B;
  double C = ((struct parameters *)params)->C;
  double *noise = ((struct parameters *)params)->noise;
  int j, i;

  //note: MKL could significantly improve this
  for(j = 2; j<N; j+=2) {
    //intrinsic dynamics and coupling
    f[j] = omega[j]+B*sin(y[j+1]-y[j])+C*sin(y[j-1]-y[j]);
    f[j+1] = omega[j+1]+C*sin(y[j+2]-y[j+1])+B*sin(y[j]-y[j+1]);

    //noise
    else if(noisetype == 1) {
      f[j] += noise[j];
      f[j+1] += noise[j+1];
    }
    else{
      f[j] += cos(y[j])*noise[j];
      f[j+1] += cos(y[j+1])*noise[j+1];
    }
  }
  //intrinsic dynamics and coupling
  f[0] = omega[0]+B*sin(y[1]-y[0])+C*sin(y[N-1]-y[0]);
  f[1] = omega[1]+C*sin(y[2]-y[1])+B*sin(y[0]-y[1]);
  f[N-2] = omega[N-2]+B*sin(y[N-1]-y[N-2])+C*sin(y[N-3]-y[N-2]);
  f[N-1] = omega[N-1]+C*sin(y[0]-y[N-1])+B*sin(y[N-2]-y[N-1]);
  //noise
  else if(noisetype == 1) {
    f[0] += noise[1];
    f[1] += noise[0];
    f[N-2] += noise[N-2];
    f[N-1] += noise[N-1];
  }
  else{
    f[0] += cos(y[0])*noise[0];
    f[1] += cos(y[1])*noise[1];
    f[N-2] += cos(y[N-2])*noise[N-2];
    f[N-1] += cos(y[N-1])*noise[N-1];
  }

  return GSL_SUCCESS;
}

int main (int argc, char* argv[]) {
  struct timeval start,end;
  int i,j,k;

  int N=2;
  double C = 0.25;
  double B = 0.25;
  double sigma = 1.0;
  int noisetype = 1;
  int correlated = 0;
  double t = 0.;
  double tmax = 1e3;
  double ta = 5e2;
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
        {"icoupling",  required_argument, 0, 'B'},
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

    opt = getopt_long (argc, argv, "n:a:cC:B:d:D:hi:o:s:t:vy:", long_options, &option_index);
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
    case 'B':
      B=atof(optarg);
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
        printf("usage: ./janus [-n NUM] [-C K] [-B BETA] [-c] [-t TIME] [-a ATIME] [-o OTIME] [-d DT] [-D DT2] [-i SIGMA] [-s SEED] [-y TYPE][-v] filebase \n\n");
        printf("optional arguments:\n");
        printf("\t -a ATIME, --animate ATIME \t\t Time to start outputting animation data. Default 1e1. \n");
        printf("\t -c, --common  \t\t Use common noise.\n");
        printf("\t -C K, --coupling C \t\t Coupling constant. Default 0.25.\n");
        printf("\t -B BETA, --icoupling B \t\t Internal coupling constant. Default 0.25.\n");
        printf("\t -d DT, --dt DT \t\t Integration timestep. Default 1e-3.\n");
        printf("\t -D DT2, --noisestep DT2 \t\t Noise sampling rate. Default 1e-2.\n");
        printf("\t -i SIGMA, --intensity SIGMA \t\t Noise intensity. Default 1.\n");
        printf("\t -n NUM, --num NUM \t\t Number of grid points in each dimension. Default 1536.\n");
        printf("\t -o OTIME, --order OTIME \t\t Time to begin averaging order parameter. Default 5e2.\n");
        printf("\t -s SEED, --seed SEED \t\t Random seed. Default 1.\n");
        printf("\t -t TIME, --time TIME \t\t Total integration time. Default 1e1.\n");
        printf("\t -v, --verbose \t\t Verbose output.\n");
        printf("\t -y, --type \t\t Noise type. 1 for additive constant phase sensitivity, 2 for additive trigonometric sensitivity. Default 1.\n");

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
  double *y, *yerr, *noise;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  FILE *out, *outsignal;
  char file[256];

  double *omega;
  omega = calloc(N, sizeof(double));
  for(j=0; j<N; j++) {
    omega[j] = pow(-1,j)*0.5;
  }

  y=calloc(N,sizeof(double));
  yerr=calloc(N,sizeof(double));
  noise=calloc(N,sizeof(double));
  struct parameters params={N, omega, C, B, noise, noisetype};

  //create the random noise function.
  gsl_rng_set(r,seed);
  for(j=0; j<N; j++) {
    double theta=2*3.14*gsl_rng_uniform(r);
    y[2*j] = theta;
  }

  //Output the noise data
  strcpy(file,filebase);
  strcat(file, ".dat");
  outsignal=fopen(file,"w");
  strcpy(file,filebase);
  strcat(file, ".out");
  out=fopen(file,"w");
  for(int i=0; i<argc; i++){
    fprintf(out,"%s ",argv[i]);
  }
  fprintf(out,"\n");
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

    if(correlated==0) {
      for(j=0; j<N; j++)
      noise[j] = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
    }
    else {
      double ran = gsl_ran_gaussian(r,sigma/sqrt(2*dt));
      for(j=0; j<N; j++)
      noise[j] = ran;
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

    if (t>=to) {
      for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
          order += cos(y[j]-y[k])/(N*N);
        }
      }
    }

    if(t>ta){
      fwrite(y, sizeof(double), N, outsignal);
      fflush(outsignal);
    }
    if(verbose){
      gettimeofday(&end,NULL);
      printf("%6f\t%6f\t%6f\t%6e\t\r",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fflush(stdout);
    }

  }

  gettimeofday(&end,NULL);
  printf("\nruntime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));
  printf("%f \n", order/(Nt-Nto));

  //Output results
  fflush(outsignal);
  fclose(outsignal);
  fprintf(out, "runtime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));
  fprintf(out, "%i %f %f %f %f %i %f \n", N, tmax-ta, dt, C, sigma, seed, order/(Nt-Nto));
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
