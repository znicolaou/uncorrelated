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
#include <fftw3.h>



struct parameters{int N; double C; double* gamma; double *noise;  };

const double V = 10.5;
// const double V = 15;
const double R = 20.0;
const double Ch = 1600.0;
const double a = 0.3;
const double b = 0.00006;
const double c = 0.001;
double gam[2]={0.0101,0.01};

double getPhase(double theta, double* Sr, double* Si){
    double ret = theta;
    for(int i = 1; i < 50; i++) {
        ret += 2*Sr[i]*(cos(i*theta)-1)-2*Si[i]*sin(i*theta);
    }
    return ret;
}

int func (double t, const double y[], double f[], void *params) {
    int N = ((struct parameters *)params)->N;
    double C = ((struct parameters *)params)->C;
    double *gamma = ((struct parameters *)params)->gamma;
    double *noise = ((struct parameters *)params)->noise;
    int i,j;

    for(j = 0; j<N; j++) {
        f[2*j] = ((V-y[2*j])/R-(Ch*exp(0.5*y[2*j])/(1+Ch*exp(y[2*j]))+a*exp(y[2*j]))*(1-y[2*j+1]));
        f[2*j+1] = ((exp(0.5*y[2*j])*(1-y[2*j+1])/(1+Ch*exp(y[2*j]))-b*Ch*exp(2*y[2*j])*y[2*j+1]/(c*Ch+exp(y[2*j])))/gamma[j]);
        for(i = 0; i<N; i++){
            f[2*j] += C*(y[2*i]-y[2*j]);
        }

        f[2*j] += noise[j]/R;
    }

    return GSL_SUCCESS;
}

int main (int argc, char* argv[]) {
  struct timeval start,end;
  int i,j,k;

  double sigma = 1.0;
  int noisetype = 1;
  int correlated = 0;
  double t = 0.;
  double tmax = 1e3;
  double ta = 5e2;
  double to = 5e2;
  double dt = 1e-1;
  double dt2 = 1e-2;
  double C=1.0;
  int seed = 1;
  char* filebase = NULL;
  char* pfname = NULL;
  int verbose = 0;
  int opt;

  while (1)
  {
    static struct option long_options[] =
      {
        {"animate",  required_argument, 0, 'a'},
        {"params",  required_argument, 0, 'p'},
        {"coupling",  no_argument, 0, 'C'},
        {"common",  no_argument, 0, 'c'},
        {"dt",  required_argument, 0, 'd'},
        {"noisestep",  required_argument, 0, 'D'},
        {"help",  no_argument, 0, 'h'},
        {"intensity",  required_argument, 0, 'i'},
        {"order",  required_argument, 0, 'o'},
        {"seed",  required_argument, 0, 's'},
        {"time",  required_argument, 0, 't'},
        {"verbose",  no_argument, 0, 'v'},
        {"type",  required_argument, 0, 'y'},
        {0, 0, 0, 0}
      };
    int option_index = 0;

    opt = getopt_long (argc, argv, "a:cC:d:D:hi:o:p:s:t:vy:", long_options, &option_index);
    if (opt == -1)
            break;
    switch (opt)
    {
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
      dt=atof(optarg);
      break;
    case 'D':
      dt2=atof(optarg);
      break;

    case 'h':
        printf("usage: ./janus [-n NUM] [-C K] [-B BETA] [-c] [-t TIME] [-a ATIME] [-o OTIME] [-d DT] [-D DT2] [-i SIGMA] [-s SEED] [-y TYPE][-v] filebase \n\n");
        printf("optional arguments:\n");
        printf("\t -a ATIME, --animate ATIME \t\t Time to start outputting animation data. Default 1e1. \n");
        printf("\t -c, --common  \t\t Use common noise.\n");
        printf("\t -C K, --coupling C \t\t Coupling constant. Default 0.25.\n");
        printf("\t -B BETA, --icoupling B \t\t Internal coupling constant. Default 0.25.\n");
        printf("\t -d DT, --dt DT \t\t Output timesteps. Default 1e-1.\n");
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
      to=atof(optarg);
      break;
    case 'p':
      pfname=optarg;
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
  int N=2;
  int Nt = (int)tmax/dt;
  int Nto = (int)to/dt;
  int Nta = (int)ta/dt;
  double ti = dt;
  double order = 0;
  double *y, *yerr, *noise, *signal, *phase, *diff;//, *gamma;
  fftw_complex *fftw_y2, *fftw_F2;
  fftw_plan iplan2, fplan2;


  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  FILE *out, *outsignal, *outphase;
  char file[256];

  y=calloc(2*N,sizeof(double));
  yerr=calloc(2*N,sizeof(double));
  noise=calloc(N,sizeof(double));
  signal=calloc(2*N*Nt,sizeof(double));
  phase=calloc(N*(Nt-Nto),sizeof(double));
  diff=calloc((Nt-Nto),sizeof(double));

  fftw_y2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Nt-Nto));
  fftw_F2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Nt-Nto));

  fplan2 = fftw_plan_dft_1d(Nt-Nto, fftw_y2, fftw_F2, FFTW_FORWARD, FFTW_ESTIMATE);
  iplan2 = fftw_plan_dft_1d(Nt-Nto, fftw_F2, fftw_y2, FFTW_BACKWARD, FFTW_ESTIMATE);

  struct parameters params={N, C, gam, noise};

  //create the random noise function.
  gsl_rng_set(r,seed);
  for(j=0; j<N; j++) {
      y[2*j] = -2.0+gsl_ran_gaussian(r,0.1);
      y[2*j+1] = 0.5;
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
  gsl_odeiv2_system sys = {func, NULL, 2*N, &params};
  gsl_odeiv2_step * step;
  step = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkf45, 2*N);

  gettimeofday(&start,NULL);

  //Do integration
  int count=0;
  double maxerr;
  while(count < Nt) {
    maxerr=0;
    if(correlated==0) {
      for(j=0; j<N; j++)
      noise[j] = gsl_ran_gaussian(r,sigma/sqrt(2*dt2));
    }
    else {
      double ran = gsl_ran_gaussian(r,sigma/sqrt(2*dt2));
      for(j=0; j<N; j++)
      noise[j] = ran;
    }
    for(k=0; k<dt/dt2; k++) {
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
    for(j=0; j<N; j++){
        signal[2*j+2*N*count] = y[2*j];
        signal[2*j+1+2*N*count] = y[2*j+1];
    }
    count++;
    t=count*dt;



    if(t>ta){
      fwrite(y, sizeof(double), 2*N, outsignal);
      fflush(outsignal);
    }
    if(verbose){
      gettimeofday(&end,NULL);
      printf("%6f\t%6f\t%6f\t%6e\t\r",t/tmax,end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec), (end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec))/(t/tmax)*(1-t/tmax), maxerr);
      fflush(stdout);
    }
  }

  strcpy(file,filebase);
  strcat(file, "phase.dat");
  outphase=fopen(file,"w");

  for(j=0; j<N; j++){
      //Find phases and order parameter
      //I will use the current ifar as the signal
      for(i=0; i<(Nt-Nto); i++){
          fftw_y2[i][0] = (Ch*exp(0.5*signal[2*j+2*N*(i+Nto)])/(1+Ch*exp(signal[2*j+2*N*(i+Nto)]))+a*exp(signal[2*j+2*N*(i+Nto)]))*(1-signal[2*j+1+2*N*(i+Nto)]);
          fftw_y2[i][1] = 0;
      }

      fftw_execute(fplan2);

      for(i=0; i<=(Nt-Nto)/2; i++) {
        fftw_F2[i][0]/=(Nt-Nto);
        fftw_F2[i][1]/=(Nt-Nto);
      }
      fftw_F2[0][0] = 0;
      fftw_F2[0][1] = 0;
      for(i=1+(Nt-Nto)/2; i<(Nt-Nto); i++) {
          fftw_F2[i][0] = 0;
          fftw_F2[i][1] = 0;
      }
      fftw_execute(iplan2);
      for(i=0; i<(Nt-Nto); i++){
        phase[i] = atan2(fftw_y2[i][1], fftw_y2[i][0]);
      }
      double Sr[50];
      double Si[50];
      int n=0;
      for(n=1; n<50; n++){
        Si[n]=0;
        Sr[n]=0;
        for(i=0; i<(Nt-Nto); i++){
          Si[n]+=-cos(n*phase[i])/(n*(Nt-Nto));
          Sr[n]+=-sin(n*phase[i])/(n*(Nt-Nto));
        }
      }

      double accumulate=0;
      phase[0] = getPhase(phase[0],Sr,Si);
      for(i=1; i<(Nt-Nto); i++){
          phase[i] = accumulate+getPhase(phase[i],Sr,Si);
          if(fabs(phase[i-1]-phase[i]) > 2*asin(1.0)) {
              accumulate += 4*asin(1.0)*((phase[i-1]-phase[i])<0?-1:1);
              phase[i] += 4*asin(1.0)*((phase[i-1]-phase[i])<0?-1:1);
          }
          diff[i] += ((j%2)==0?1:-1)*phase[i];
      }
      fwrite(phase, sizeof(double), (Nt-Nto), outphase);
  }
  fclose(outphase);
  for(i=0; i<(Nt-Nto); i++) {
      order += pow(cos(diff[i]/2.0),2.0);
  }

  gettimeofday(&end,NULL);
  printf("\nruntime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));
  printf("%f \n", order/(Nt-Nto));

  //Output results
  fflush(outsignal);
  fclose(outsignal);
  fprintf(out, "runtime: %6f\n",end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec-start.tv_usec));
  fprintf(out, "%i %i %f %i %f %f \n", N, (Nt-Nta), tmax-ta, seed, sigma, order/(Nt-Nto));
  fflush(out);
  fclose(out);

  gsl_odeiv2_step_free(step);
  gsl_rng_free(r);
  fftw_destroy_plan(iplan2);
  fftw_destroy_plan(fplan2);
  free(y);
  free(yerr);
  free(noise);
  free(signal);
  free(phase);
  free(diff);
  return 0;
}
