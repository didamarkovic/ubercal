#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

const int NX = 20;   // number of exposures in x - direction
const int NY = 20;   // number of exposures in y - direction
const int NDITH = 4; // number of dithers

struct overlap {
  double area;       // area of overlap
  int nexposure;
  int iexposure[NDITH]; // integers of exposures included in this overlap - max of one from each dither
};

// std error handler
void err_known(const char *error_text)
{
  fprintf(stderr,"Run-time error...\n%s\n",error_text);
  fprintf(stderr,"now exiting to system...\n");
  exit(1);
}

int main() {

  const int bsz=800; char fname[bsz];

  // ********************************************************************************
  // memory for overlap polygons as input
  const int NSMALL_MAX = 30;
  struct overlap *p_small;
  if(!(p_small = (struct overlap*)malloc(NSMALL_MAX*sizeof(struct overlap)))) 
    err_known("memory allocation problem for polygons");

  // ********************************************************************************
  // read in survey mask element file
  
  FILE *fin_overlaps;
  sprintf(fname,"full-pointing-overlaps.txt");
  if((fin_overlaps=fopen(fname,"r"))==NULL)
    err_known("cannot open input file");
  int Ndither_tmp, nsmall;
  fscanf(fin_overlaps,"# %d dithers\n",&Ndither_tmp);
  if(Ndither_tmp!=NDITH) err_known("dither number mismatch");

  fscanf(fin_overlaps,"# %d overlap polygons\n",&nsmall);
  if(nsmall>NSMALL_MAX) err_known("input data size");
  
  for(int ismall=0;ismall<nsmall;ismall++) {

    fscanf(fin_overlaps,"%*d %lf %d :",&p_small[ismall].area,&p_small[ismall].nexposure);
    for(int ie=0;ie<p_small[ismall].nexposure;ie++) fscanf(fin_overlaps,"%d ",&p_small[ismall].iexposure[ie]);
    fscanf(fin_overlaps,"\n");
  }
  fclose(fin_overlaps);
  printf("Read in %d overlaps in patch\n",nsmall);
  
  // ********************************************************************************
  // build up large matrix - this is N_overlaps x N_exposures big
  // However, it's sparse, so store number and indices of exposures,
  // rather than filling the full matrix

  // ********************************************************************************
  // memory for overlap polygons as input
  struct overlap *p_large;
  if(!(p_large = (struct overlap*)malloc(NX*NY*nsmall*sizeof(struct overlap)))) 
    err_known("memory allocation problem for polygons");

  // ********************************************************************************
  // make list of all overlaps
  int ioverlap_large=0;
  for(int iy=0;iy<NY;iy++) for(int ix=0;ix<NX;ix++) 
      for(int ismall=0;ismall<nsmall;ismall++) {

	// remove overlaps that will be double counted
	// these extend to patches that will be counted again within the pattern
	int good_overlap=1, nexposure=0;
	p_large[ioverlap_large].area = p_small[ismall].area;
		
	for(int idither=0;idither<NDITH;idither++)
	  for(int jpoint=-1;jpoint<=1;jpoint++)
            for(int ipoint=-1;ipoint<=1;ipoint++) {
	      
	      int ipoly_small = idither*9 + (jpoint+1)*3 + (ipoint+1);
	      int jy = iy + jpoint;
	      int jx = ix + ipoint;
	      int ipoly_large = idither*NX*NY + jy*NX + jx;

	      for(int ie=0;ie<p_small[ismall].nexposure;ie++)
		if(p_small[ismall].iexposure[ie] == ipoly_small) {

		  // remove overlaps looking forward
		  //these will be covered from another large loop
		  if( (jx<NX && jx>ix) || (jy<NY && jy>iy) ) good_overlap=0;
		  
		  // only count large polygons within survey
		  if(jx>=0 && jx<NX && jy>=0 && jy<NY) {
		    p_large[ioverlap_large].iexposure[nexposure]=ipoly_large;
		  nexposure++;
		  }
		}
	    }
	p_large[ioverlap_large].nexposure = nexposure;
	
	if(nexposure>1 && good_overlap) ioverlap_large++;
      }
		  
  // ********************************************************************************
  // print out final list of overlaps for "full survey"

  FILE *fout_survey;
  sprintf(fname,"full-survey-overlaps.txt");
  if((fout_survey=fopen(fname,"w"))==NULL)
    err_known("cannot open output file");
  fprintf(fout_survey,"# %d dithers\n",NDITH);
  fprintf(fout_survey,"# NX=%d, NY=%d\n",NX,NY);
  fprintf(fout_survey,"# %d overlap polygons\n",ioverlap_large);
  for(int i=0;i<ioverlap_large;i++) {
    fprintf(fout_survey,"%d %g %d :",i,p_large[i].area,p_large[i].nexposure);
    for(int ie=0;ie<p_large[i].nexposure;ie++) fprintf(fout_survey,"%d ",p_large[i].iexposure[ie]);
    fprintf(fout_survey,"\n");
  }
  fclose(fout_survey);
  
  exit(0);
}


