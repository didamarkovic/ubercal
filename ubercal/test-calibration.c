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

const int VERB = 0;
FILE *FLIKE;

struct overlap_large {
  double area;             // area of overlap
  double ssq_calib;        // sigma^2 for calibrators in this exposure covering the overlap region
  int nexposure;           // number of exposures
  int iexposure[4];     // integers of exposure IDs included in this overlap
  double flux_calib[4]; // normalised total flux measured for calibrators in this overlap
                           // - each exposure has its own measurement
}; // TODO: now it is ndetector & idetector, not exposure!
const int EMAX = 4;  // maximum number of dithers initiate in the overlap_large struct - test for this on file input

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
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

// function to calculate chi^2
double calc_chisq(double*);

// some global variables that save effort
const double SIG_INIT = 0.16;                // 4% initial exposure calibration (16% detector)
const double SIG_FINAL = 0.0;                // 0% expected final calibration (?)
long nexposure, noverlap;

// buffer size
const int bsz=800; 

// calibration star densities and rms from Marco - simplified version of Dida's algorithm.
//const int NBIN_STAR_DEF=7;
//double rms_v_DEF[7] = {0.00193, 0.00252, 0.00356, 0.00642, 0.01594, 0.04102, 0.11044};
//double dens_v_DEF[7] = {31.90, 60.20, 97.50, 156.30, 234.30, 332.10, 428.80};
double *rms_v, *dens_v;

int main(int argc, char *argv[]) {

  char fname[bsz], fpath[bsz], fstar[bsz]; long seed = -1;
  int NDETX;

  // ********************************************************************************
  // read input parameter for random seed
  if (argc > 1){
    sscanf(argv[1], "%s", fpath);
    sscanf(argv[2], "%ld", &seed);
  } else {
    sscanf("./", "%s", fpath);
  }
  if (argc > 3) {
    sscanf(argv[3], "%s", fstar); 
  } else {
    sscanf("./stars.dat", "%s", fstar);
  }
  if(VERB>1) printf("Results in %s.\n",fpath);
  if(VERB>0) printf("random seed = %ld\nStars in %s.\n",seed,fstar);
  fflush(stdout);  

  // ********************************************************************************
  // read in stellar densities and RMS if file input
  long NBIN_STAR;

  FILE *fin_stars;
  strcpy(fname, fstar);
  if((fin_stars=fopen(fname,"r"))==NULL) {
    strcat(fname, " <- test-calibration cannot open star file");
    err_known(fname);
  }
  fscanf(fin_stars,"# %ld star types\n",&NBIN_STAR);
  if(VERB>1) printf("Read in %ld stellar types.\n",NBIN_STAR);

  // memory for stars input
  if(!(rms_v = (double*)malloc(NBIN_STAR*sizeof(double))))
    err_known("memory allocation problem for RMS");
  if(!(dens_v = (double*)malloc(NBIN_STAR*sizeof(double))))
    err_known("memory allocation problem for density");
  
  // Read stars from file
  for(int i=0;i<NBIN_STAR;i++) {
    if (fscanf(fin_stars,"%*d %lf %lf \n", &dens_v[i], &rms_v[i]) != 2){
      strcat(fname, " <- test-calibration cannot read star file");
      err_known(fname);
    }
    if(VERB>1) printf("RMS=%f, dens=%f\n", rms_v[i], dens_v[i]);
  }
  fclose(fin_stars);
  fname[0] = '\0';

  // ********************************************************************************
  // read in survey mask element file
  long NX, NY, NDITH;

  FILE *fin_survey;
  strcpy(fname, fpath);
  strcat(fname, "full-survey-overlaps.txt");
  if((fin_survey=fopen(fname,"r"))==NULL) {
    strcat(fname, " <- test-calibration cannot open input file");
    err_known(fname);
    }
  fscanf(fin_survey,"# %ld dithers, %d detectors squared\n",&NDITH,&NDETX);
  if(NDITH>EMAX) err_known("NDITH > EMAX");
  
  fscanf(fin_survey,"# NRA=%ld, NDEC=%ld\n",&NX,&NY);
  nexposure = NDETX*NDETX*NDITH*NX*NY;
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
  fname[0] = '\0';

  if(VERB>0) printf("Read in %ld overlaps in survey in %ld parent polygons.\n",noverlap,nexposure);
  
  // ********************************************************************************
  // set up initial calibration - all calibrators have these fluctuations in addition
  // their intrinsic measurement noise
  double old_calib[nexposure+1], new_calib[nexposure+1];
  for(int i=1;i<=nexposure;i++) old_calib[i] = SIG_INIT*gasdev(&seed);
  
  // ********************************************************************************
  // set up mock flux measurements for each exposure of calibrators in overlaps
  for(int i=0;i<noverlap;i++) {

    // each star gives rise to a measurement with variance (rms_v[istar])^2
    // loop over magnitude bins, and Poisson sampling each (to contain nstar stars),
    // gives a variance of the final variance of the calibration measurement
    // of 1/(sum_bins nstar/var)
    double nstar_tot=0.0, sigma = 0.0;
    for(int istar=0;istar<NBIN_STAR;istar++) {
      double nstar = poidev(p_overlap[i].area*dens_v[istar], &seed);
      sigma += nstar/(rms_v[istar]*rms_v[istar]);
      nstar_tot += nstar;
    }
    if(nstar_tot>0 && sigma<1.0e100) sigma = 1./sqrt(sigma);
    else            p_overlap[i].nexposure=0;

    if (VERB>1) printf("In overlap %d, the number of stars = %f and their RMS = %f.\n",i, nstar_tot, sigma);

    // actual flux measurement for each exposure in overlap, drawn from the same
    // distribution for the intrinsic noise, plus the initial calibration of that exposure
    for(int ie=0;ie<p_overlap[i].nexposure;ie++) {
      if(nstar_tot>0) {
        p_overlap[i].flux_calib[ie] = sigma*gasdev(&seed) + old_calib[p_overlap[i].iexposure[ie]+1];
      }
      else            p_overlap[i].flux_calib[ie] = 0.0;
      if(VERB>2) printf("The %g stars contribute to the f_meas = %g.\n", nstar_tot, p_overlap[i].flux_calib[ie]);
    }
    // variance of distribution of measured fluxes around mean
    // the (N-1)/N factor allows for decreased scatter around measured mean
    p_overlap[i].ssq_calib = (sigma*sigma+SIG_FINAL*SIG_FINAL) * (p_overlap[i].nexposure-1)/p_overlap[i].nexposure;
  }

  // exit(0);
  
  // ********************************************************************************
  // perform minimisation

  // slightly faster if you start with new_calib as the perfect answer (=-old_calib[i])
  // but this doesn't seem fair!
  for(int i=1;i<=nexposure;i++) new_calib[i]=0.0;

  // open file to write the minimisation steps
  if(VERB>1){
    strcpy(fname, fpath);
    strcat(fname,"likelihoods.txt");
    if((FLIKE=fopen(fname,"w"))==NULL) {
      strcat(fname, " <- c-e-p cannot open output file");
      err_known(fname);
    }
    fprintf(FLIKE,"# last line would be perfect, last two columns are the total chisq and the 'prior'\n");
  }

  // internal parameters required by powell
  double **xi = dmatrix(1,nexposure,1,nexposure);
  for(int i=1;i<=nexposure;i++)
    for(int j=1;j<=nexposure;j++) 
      xi[i][j]=(i==j?0.1*SIG_INIT:0.0);
  double endval=0.0;
  int itmp;

  // minimisation
  double ftol = 1e-7;
  powell(new_calib,xi,nexposure,ftol,&itmp,&endval,calc_chisq);
  free_dmatrix(xi,1,nexposure,1,nexposure);
  if (VERB>0) printf("After %i iterations, the final chi^2 = %g.\n", itmp, endval);

  // print truth and close file
  if(VERB>1){
    for(int i=1;i<=nexposure;i++) fprintf(FLIKE," %+1.3f", -old_calib[i]);
    fprintf(FLIKE," 0.0 0.0\n");
    fprintf(FLIKE,"\n");
    fclose(FLIKE);
  }

  // ********************************************************************************
  // write the results of the minimisation to a file
  FILE *fout;
  strcpy(fname, fpath);
  strcat(fname,"calibrations.txt");
  if((fout=fopen(fname,"w"))==NULL) {
    strcat(fname, " <- c-e-p cannot open output file");
    err_known(fname);
    }
  fprintf(fout,"# %f initial scatter\n", SIG_INIT);
  fprintf(fout,"# exposure-id initial-zero-point calibration-correction calibrated-zero-point\n");
  for(int i=1;i<=nexposure;i++){
      fprintf(fout,"%i %g %g %g\n", i, old_calib[i], new_calib[i], old_calib[i]+new_calib[i]);
  }
  fprintf(fout,"\n");
  fclose(fout);

  // ********************************************************************************
  // now test post-ubercal calibrations - these are the initial calibrations
  // for each exposure as stored in old_calib + the best-fit new calibrations in new_calib
  double mean_init=0.0, sigmasq_init=0.0;
  double mean_final=0.0, sigmasq_final=0.0;
  double mean_cal=0.0, sigmasq_cal=0.0;
  for(int i=1;i<=nexposure;i++) {
    if(VERB>1) printf("Exposure %d, old flux %g, calibration %g, overlap mean %g\n",i,old_calib[i],new_calib[i], old_calib[i]+new_calib[i]);
    mean_init  += old_calib[i]; 
    mean_final += old_calib[i]+new_calib[i]; 
    mean_cal += new_calib[i];
    sigmasq_init  += old_calib[i]*old_calib[i];
    sigmasq_final += (old_calib[i]+new_calib[i])*(old_calib[i]+new_calib[i]);
    sigmasq_cal += new_calib[i]*new_calib[i];
  }
  mean_init  /= (double)nexposure;
  mean_final /= (double)nexposure;
  mean_cal /= (double)nexposure;
  if(VERB>0) printf("Improvement due to Ubercal, q = %g and\n", sqrt(sigmasq_init/sigmasq_final)); 
  sigmasq_init  = sigmasq_init/(double)nexposure;
  sigmasq_final = sigmasq_final/(double)nexposure;
  sigmasq_cal = sigmasq_cal/(double)nexposure;

  // Now calculate the true covariance between the true initial zero-points and the calibrations
  if(VERB>0){
    double cov_init_cal_diag=0.0;
    for(int i=1;i<=nexposure;i++) {
      cov_init_cal_diag += (new_calib[i]-mean_cal)*(old_calib[i]-mean_init);
    }
    cov_init_cal_diag/=(double)nexposure;
    printf("True correlation between the true zero-point and the final calibration correction = %g\n",
     cov_init_cal_diag/sqrt((sigmasq_cal-mean_cal*mean_cal)*(sigmasq_init-mean_init*mean_init)));
  }

  if(VERB>0) printf("Initial 0-point %g, final 0-point %g, calibrated by %g\n", mean_init, mean_final, mean_cal);

  printf("Initial scatter %g, final scatter %g, calibration scatter %g\n",
   sqrt(sigmasq_init), sqrt(sigmasq_final), sqrt(sigmasq_cal));

  exit(0);
}

// this is the main function calculating chi^2
double calc_chisq(double *new_calib) {

  if(VERB>2) printf("\n--- NEW ITERATION ---\n"); // diagnostics

  // this is the prior - we only expect calibrations of order the initial ones
  double chisq0=0.0, chisq=0.0;
  for(int i=1;i<=nexposure;i++){
    chisq += new_calib[i]*new_calib[i];
    if(VERB>1) fprintf(FLIKE," %+1.3f", new_calib[i]); // diagnostics
  }
  chisq /= (SIG_INIT*SIG_INIT);
  if(VERB>1) chisq0 = chisq; // diagnostics

  if(VERB>2) printf("\tchisq_0 = %g\n", chisq); // diagnostics

  // this holds the new flux measurements for each calibrator in each overlap
  double *new_flux = malloc(EMAX*sizeof(double));

  // this is the big loop over overlaps - matching fluxes in each is ubercal!
  for(int i=0;i<noverlap;i++) if(p_overlap[i].nexposure>1) {

      // new flux calibrations in overlaps
      for(int ie=0;ie<p_overlap[i].nexposure;ie++)
        new_flux[ie] = p_overlap[i].flux_calib[ie]+new_calib[p_overlap[i].iexposure[ie]+1];
      
      // mean flux in each overlap//
      double mean=0.0;
      for(int ie=0;ie<p_overlap[i].nexposure;ie++) {
        mean += new_flux[ie];
      }
      mean /= (double)p_overlap[i].nexposure;
      if(VERB>2) printf("\tmean in overlap %d = %g\n", i, mean); // diagnostics
      
      // add to chi^2 for this overlap region from differences with mean
      // relies on ssq.calib being covariance of exposure calibration measurements
      // around their sample mean
      for(int je=0;je<p_overlap[i].nexposure;je++) { 
        double diff_j = new_flux[je]-mean;

        for(int ie=0;ie<p_overlap[i].nexposure;ie++) {

          // covariance, which becomes variance*Bessel factor for ie=je
          double cov_ieje = -p_overlap[i].ssq_calib/p_overlap[i].nexposure;
          if(ie==je) cov_ieje += p_overlap[i].ssq_calib;

          double diff_i = new_flux[ie]-mean;
          chisq += diff_i*diff_j / cov_ieje;
        }

      }
      
    }
  free(new_flux);

  // diagnostics
  if(VERB>1) {
    fprintf(FLIKE," %1.2e", chisq);
    fprintf(FLIKE," %1.2e", chisq0);
    fprintf(FLIKE,"\n");
  }
  if(VERB>2) printf("\tchisq_tot = %g at sig_i = %g\n", chisq, SIG_INIT);

  return chisq;
}

