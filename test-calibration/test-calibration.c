/*
Code by Will Percival & Dida Markovic

Uses stellar numbers from Marco Scodeggio

Not to be used except for internal Euclid work unless approved:

Send an email to Dida Markovic (dida.markovic@port.ac.uk)
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

const int EMAX = 4; // maximum number of dithers - test for this on file input

struct overlap_large {
  double area;             // area of overlap
  double ssq_calib;        // sigma^2 for calibrators in this exposure covering the overlap region
  int nexposure;           // number of exposures
  int iexposure[EMAX];     // integers of exposure IDs included in this overlap
  double flux_calib[EMAX]; // normalised total flux measured for calibrators in this overlap
                           // - each exposure has its own measurement
};
struct overlap_large *p_overlap;

// std error handler
void err_known(const char *error_text)
{
  fprintf(stderr,"Run-time error...\n%s\n",error_text);
  fprintf(stderr,"now exiting to system...\n");
  exit(1);
}

// numerical recipes routines used
void powell(double[],double**,int,double,int*,double*,double (*func)(double []));
float gasdev(long*);
float poidev(float,long*);
double **dmatrix(int,int,int,int);

// function to calculate chi^2
double calc_chisq(double*);

// some global variables that save effort
const double SIG_INIT = 0.04;                // 4% initial calibration
long nexposure, noverlap;

// calibration star densities and rms from Marco - simplified version of Dida's algorithm.
const int NBIN_STAR=7;
double rms_v[NBIN_STAR] = {0.00193, 0.00252, 0.00356, 0.00642, 0.01594, 0.04102, 0.11044};
double dens_v[NBIN_STAR] = {319.0, 602.0, 975.0, 1563.0, 2343.0, 3321.0, 4288.0};

int main() {

  const int bsz=800; char fname[bsz];

  // ********************************************************************************
  // read in survey mask element file
  long NX, NY, NDITH;

  FILE *fin_survey;
  sprintf(fname,"full-survey-overlaps.txt");
  if((fin_survey=fopen(fname,"r"))==NULL) err_known("cannot open input file");
  fscanf(fin_survey,"# %ld dithers\n",&NDITH);
  if(NDITH>EMAX) err_known("NDITH > EMAX");
  
  fscanf(fin_survey,"# NRA=%ld, NDEC=%ld\n",&NX,&NY);
  nexposure = 16.0*NDITH*NX*NY;
  fscanf(fin_survey,"# %ld overlap polygons\n",&noverlap);

   // memory for overlap polygons input
  if(!(p_overlap = (struct overlap_large*)malloc(noverlap*sizeof(struct overlap_large)))) 
    err_known("memory allocation problem for overlaps");
  
  for(int i=0;i<noverlap;i++) {
    fscanf(fin_survey,"%*d %lf %d :",&p_overlap[i].area,&p_overlap[i].nexposure);
    for(int ie=0;ie<p_overlap[i].nexposure;ie++) fscanf(fin_survey,"%d ",&p_overlap[i].iexposure[ie]);
    fscanf(fin_survey,"\n");
  }
  fclose(fin_survey);

  printf("Read in %ld overlaps in survey\n",noverlap);
  
  // ********************************************************************************
  // set up initial calibration - all calibrators have these fluctuations in addition
  // their intrinsic measurement noise
  long seed=-1;
  double old_calib[nexposure+1], new_calib[nexposure+1];
  for(int i=1;i<=nexposure;i++) old_calib[i] = SIG_INIT*gasdev(&seed);

  // slightly faster if you start with new_calib as the perfect answer (=-old_calib[i])
  // but this doesn't seem fair!
  for(int i=1;i<=nexposure;i++) new_calib[i]=0.0;
 
  // ********************************************************************************
  // set up mock flux measurements for each exposure of calibrators in overlaps
  for(int i=0;i<noverlap;i++) {

    // each star gives rise to a measurement with variance (rms_v[istar])^2
    // loop over magnitude bins, and Poisson sampling each (to contain nstar stars),
    // gives a variance of the final variance of the calibration measurement
    // of 1/(sum_bins nstar/var)
    double sigma = 0.0, nstar_tot=0.0;
    for(int istar=0;istar<NBIN_STAR;istar++) {
      double nstar = poidev(p_overlap[i].area*dens_v[istar],&seed);
      sigma += nstar/(rms_v[istar]*rms_v[istar]);
      nstar_tot += nstar;
    }
    if(nstar_tot>0) sigma = 1./sqrt(sigma);
    else            p_overlap[i].nexposure=0;

    // actual flux measurement for each exposure in overlap, drawn from the same
    // distribution for the intrinsic noise, plus the initial calibration of that exposure
    for(int ie=0;ie<p_overlap[i].nexposure;ie++) {
      if(nstar_tot>0) p_overlap[i].flux_calib[ie] = sigma*gasdev(&seed) + old_calib[p_overlap[i].iexposure[ie]+1];
      else            p_overlap[i].flux_calib[ie] = 0.0;
    }
    
    // variance of distribution of measured fluxes around mean
    // the (N-1)/N factor allows for decreased scatter around measured mean
    p_overlap[i].ssq_calib = (p_overlap[i].nexposure-1.)/p_overlap[i].nexposure * (sigma*sigma+SIG_INIT*SIG_INIT);
    
  }

  // exit(0);
  
  // ********************************************************************************
  // perform minimisation

  // internal parameters required by powell
  double **xi = dmatrix(1,nexposure,1,nexposure);
  for(int i=1;i<=nexposure;i++)
    for(int j=1;j<=nexposure;j++) 
      xi[i][j]=(i==j?0.1*SIG_INIT:0.0);
  double endval=0.0;
  int itmp;
  powell(new_calib,xi,nexposure,1.0e-3,&itmp,&endval,calc_chisq);

  // ********************************************************************************
  // now test post-ubercal calibrations - these are the initial calibrations
  // for each exposure as stored in old_calib + the best-fit new calibrations in new_calib
  double sigmasq_init=0.0, sigmasq_final=0.0;
  double mean_init=0.0, mean_final=0.0;
  for(int i=1;i<=nexposure;i++) {
    printf("Exposure %d, calib old %g, new %g\n",i,old_calib[i],old_calib[i]+new_calib[i]);
    mean_init  += old_calib[i];
    mean_final += old_calib[i]+new_calib[i]; 
    sigmasq_init  += old_calib[i]*old_calib[i];
    sigmasq_final += (old_calib[i]+new_calib[i])*(old_calib[i]+new_calib[i]);
  }
  mean_init  /= (double)nexposure;
  mean_final /= (double)nexposure;

  printf("Initial calibration %g, final calibration %g\n",
	 sqrt(fabs(sigmasq_init/(double)nexposure-mean_init*mean_init)),
	 sqrt(fabs(sigmasq_final/(double)nexposure-mean_final*mean_final)));

  exit(0);
}

// this is the main function calculating chi^2
double calc_chisq(double *new_calib) {

  // this is the prior - we only expect calibrations of order the initial ones
  double chisq=0.0;
  for(int i=1;i<=nexposure;i++) chisq += new_calib[i]*new_calib[i];
  chisq /= (SIG_INIT*SIG_INIT);

  // this holds the new flux measurements for each calibrator in each overlap
  double *new_flux = malloc(EMAX*sizeof(double));
  
  // this is the big loop over overlaps - matching fluxes in each is ubercal!
  for(int i=0;i<noverlap;i++) if(p_overlap[i].nexposure>1) {

      // new flux calibrations in overlaps
      for(int ie=0;ie<p_overlap[i].nexposure;ie++)
	new_flux[ie] = p_overlap[i].flux_calib[ie]+new_calib[p_overlap[i].iexposure[ie]+1];
      
      // mean flux in each overlap
      double mean=0.0;
      for(int ie=0;ie<p_overlap[i].nexposure;ie++) mean += new_flux[ie];
      mean /= (double)p_overlap[i].nexposure;
      
      // add to chi^2 for this overlap region from differences with mean
      // relies on ssq.calib being variance of exposure calibration measurements
      // around mean
      for(int ie=0;ie<p_overlap[i].nexposure;ie++) {
	double diff = new_flux[ie]-mean;
	chisq += diff*diff / p_overlap[i].ssq_calib;
      }
      
    }
  free(new_flux);
  
  return chisq;
}

