#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

const double npix=2000.0;        // number of pixels along side of chip
const double pix_size=1.8e-5;    // pixel size in m
const double xgap_size=3.0e-3;   // x-direction chip gap size in m
const double ygap_size=6.0e-3;   // y-direction chip gap size in m
const double pix_scale=0.3;      // pixel scale in arcsec
const double xstart=180.*3600.;  // initial point for x in arcsec

const double XBIT=0.1; // offet allowed when polygon matching in arcsec - corresponds to 1/3 of a pixel;

const double mpi = 3.1415926536;

struct poly {
  double area;
  double weight;
  double xlow, ylow, xhigh, yhigh;
  int iexposure;
};

// std error handler
void err_known(const char *error_text)
{
  fprintf(stderr,"Run-time error...\n%s\n",error_text);
  fprintf(stderr,"now exiting to system...\n");
  exit(1);
}

// the six numbers that give the dither strategy
const double dither[6]={50.0,100.0,0.0,100.0,0.0,100.0};
// const double dither[6]={275.0,0.0,0.0,350.0,-275.0,0.0};
// const double dither[6]={1599.0,0.0,0.0,1774.0,-1599.0,0.0};

int main() {

  const int bsz=800; char buf[bsz], fname[bsz];

  // ********************************************************************************
  // memory for polygons - this is overkill for final ones, as we don't store vertices
  long NPOLY_MAX = 10000;
  struct poly *p_init;
  if(!(p_init = (struct poly*)malloc(NPOLY_MAX*sizeof(struct poly)))) 
    err_known("memory allocation problem for polygons");

  struct poly *p_final;
  if(!(p_final = (struct poly*)malloc(NPOLY_MAX*sizeof(struct poly)))) 
    err_known("memory allocation problem for polygons");

  // ********************************************************************************
  // need to work out relative positions of all CCDs in relevent exposures
  // this is where we put in the dither information - number, offsets
  const int Ndither = 4;

  // set up dither offsets
  double xdith[Ndither], ydith[Ndither];
  for(int idither=0;idither<Ndither;idither++) {
    if(idither==0) { xdith[idither]=0.0; ydith[idither]=0.0; }
    if(idither==1) { xdith[idither]=xdith[idither-1]+dither[0]; ydith[idither]=ydith[idither-1]+dither[1]; }
    if(idither==2) { xdith[idither]=xdith[idither-1]+dither[2]; ydith[idither]=ydith[idither-1]+dither[3]; }
    if(idither==3) { xdith[idither]=xdith[idither-1]+dither[4]; ydith[idither]=ydith[idither-1]+dither[5]; }
  }

  // shifts in x and y marking CCD edges
  double xp_shift = 4.0*(npix + xgap_size/pix_size)*pix_scale;
  double yp_shift = 4.0*(npix + ygap_size/pix_size)*pix_scale;

  int ipoly=0;    
  for(int idither=0;idither<Ndither;idither++) {
       
    // loop through neigbouring exposures
    // each pointing has 9 of these, 3 in ipoint and 3 in jpoint
    for(int jpoint=-1;jpoint<=1;jpoint++) {
      
      double ypoint = (double)jpoint * yp_shift;
      
      for(int ipoint=-1;ipoint<=1;ipoint++) {
	
	double xpoint = (double)ipoint * xp_shift;

	// loop over CCDs in each exposure - polynomial for each
	for(int i=0;i<4;i++)	    
	  for(int j=0;j<4;j++) {
	    p_init[ipoly].xlow      = xstart + xpoint + xdith[idither] + (double)i*(npix + xgap_size/pix_size)*pix_scale;
	    p_init[ipoly].xhigh     = p_init[ipoly].xlow + npix*pix_scale;	    
	    p_init[ipoly].ylow      = ypoint + ydith[idither] + (double)j*(npix + ygap_size/pix_size)*pix_scale;
	    p_init[ipoly].yhigh     = p_init[ipoly].ylow + npix*pix_scale;	      
	    p_init[ipoly].area      = 0.0;
	    p_init[ipoly].iexposure = idither*9 + (jpoint+1)*3 + (ipoint+1);
	    p_init[ipoly].weight    = (float)p_init[ipoly].iexposure;
	    if(++ipoly>NPOLY_MAX) err_known("too many initial polynomials for memory");
	  }
	
      }
    }
  }
  long npoly_init = ipoly;

  // ********************************************************************************
  // output file of polygons
  
  FILE *fout_mask;
  sprintf(fname,"full-pointing.vrt");
  if((fout_mask=fopen(fname,"w"))==NULL)
    err_known("cannot open output file");
  
  fprintf(fout_mask,"%ld polygons corresponding to the focal plane of NISP\n",npoly_init);
  fprintf(fout_mask,"unit s\n"); // arcsec

  for(int ipoly=0;ipoly<npoly_init;ipoly++) {
    fprintf(fout_mask,"vertices 0 ( 4 vertices, %ld weight, %lf %lf mid):\n",
	    (long)p_init[ipoly].weight,
	    0.5*(p_init[ipoly].xlow+p_init[ipoly].xhigh),
	    0.5*(p_init[ipoly].ylow+p_init[ipoly].yhigh)); // arcsec
    fprintf(fout_mask,"  %lf %lf\t%lf %lf\t%lf %lf\t%lf %lf\n",
	    p_init[ipoly].xlow,p_init[ipoly].yhigh,
	    p_init[ipoly].xlow,p_init[ipoly].ylow,
	    p_init[ipoly].xhigh,p_init[ipoly].ylow,
	    p_init[ipoly].xhigh,p_init[ipoly].yhigh);
  }
  fclose(fout_mask);
  
  // ********************************************************************************

  // now run mangle to Balkanise mask
  system("./mangle2.2/bin/poly2poly full-pointing.vrt full-pointing.pol");
  system("./mangle2.2/bin/pixelize full-pointing.pol full-pointing.pol");
  system("./mangle2.2/bin/snap full-pointing.pol full-pointing.pol");
  system("./mangle2.2/bin/balkanize -Ba full-pointing.pol full-pointing.pol");
  
  // Convert output back to vertex file
  system("./mangle2.2/bin/poly2poly -ovs full-pointing.pol full-pointing.vrt");

  // exit(0);
  
  // ********************************************************************************
  // now read in, and sum over multiple regions
  FILE *fp;
  if((fp=fopen(fname,"r"))==NULL)
    err_known("cannot open input file");

  // read in total number of polygons
  long npoly_final=0;
  fscanf(fp,"%ld polygons\n",&npoly_final);
  printf("%ld polygons in balkanised file\n",npoly_final);
  if(npoly_final>NPOLY_MAX) err_known("too many polygons in file");
  fgets(buf,bsz,fp);

  // this is the vestor we use to match final "small" polygons to the exposures, of which we consider Ndither*9 
  int poly_match[npoly_final][Ndither*9];
  for(int i=0;i<npoly_final;i++) for(int j=0;j<Ndither*9;j++) poly_match[i][j]=0;

  for(int ipoly=0;ipoly<npoly_final;ipoly++) {

    // read in head for each polygon
    fgets(buf,bsz,fp);

    int nvert;
    if(sscanf(buf,"vertices %*d ( %d vertices",&nvert)!=1) {
      printf("%s",buf);
      err_known("polygon input problem");
    }
    // if(nvert!=4) err_known("polygon size problem");

    // read in vertices for each polygon
    double x[nvert],y[nvert];
    for(int iv=0;iv<nvert;iv++) fscanf(fp,"%lf %lf",&x[iv],&y[iv]);
    fgets(buf,bsz,fp);

    // test if all final polygon vertices are within initial parent polygon
    // and flag polygons in Balkanised mask that are fully inside one of the initial polygons
    for(int j=0;j<npoly_init;j++) {
      int parent=1;
      for(int iv=0;iv<nvert;iv++) 
	if(p_init[j].xhigh < (x[iv]-XBIT) ||
	   p_init[j].xlow  > (x[iv]+XBIT) ||
	   p_init[j].yhigh < (y[iv]-XBIT) || 
	   p_init[j].ylow  > (y[iv]+XBIT) ) parent=0;
      
      if(parent) poly_match[ipoly][p_init[j].iexposure]=1;
    }
    
    // calculate area in deg^2
    p_final[ipoly].area = 0.0;
    int j=nvert-1;
    for(int iv=0;iv<nvert;iv++) {
      p_final[ipoly].area += (x[iv]+x[j]) * (y[iv]-y[j]);
      j=iv;
    }
    p_final[ipoly].area *= 0.5/(3600.*3600.);

  }
  fclose(fp);

  // exit(0);
  
  // ********************************************************************************
  // sort final polygons to reduce the number to unique sets of overlaps
  // after this the positions of the vertices become unimportant
  for(int i=npoly_final;i>=0;i--)
    for(int j=0;j<i;j++) {
      int match=1;
      for(int k=0;k<Ndither*9;k++) if(poly_match[i][k]!=poly_match[j][k]) match=0;
      if(match) {
	p_final[j].area+=p_final[i].area;
	p_final[i].area=0.0;
	break;
      }
    }

  // ********************************************************************************
  // final output of reduced set of overlap areas
  long noverlap=0, noverlap_central=0;
  for(int i=0;i<npoly_final;i++) if(p_final[i].area>0.0) {
      noverlap++;
      if(poly_match[i][4]==1 || poly_match[i][13]==1 || poly_match[i][22]==1 || poly_match[i][31]==1)
	noverlap_central++;
    }
  printf("%ld overlap regions, %ld using central exposures\n",
	 noverlap,noverlap_central);

  // check total areas agree with expectation
  for(int i=0;i<Ndither;i++) {
    int cen_exposure = 4 + 9*i;
    double tot_area=0.0;
    for(int i=0;i<npoly_final;i++) if(poly_match[i][cen_exposure]==1) tot_area += p_final[i].area;
    printf("exposure %d, total area %g deg^2\n",cen_exposure,tot_area);
  }


  FILE *fout_overlaps;
  sprintf(fname,"full-pointing-overlaps.txt");
  if((fout_overlaps=fopen(fname,"w"))==NULL)
    err_known("cannot open output file");
  fprintf(fout_overlaps,"# %d dithers\n",Ndither);
  fprintf(fout_overlaps,"# %ld overlap polygons\n",noverlap_central);

  //find total_area of central region
  double tot_area = xp_shift * yp_shift / (3600.*3600.);
  double stack_area[5];
  for(int i=1;i<5;i++) stack_area[i]=0.0;

  for(int l=0,i=0;i<npoly_final;i++)
    if(p_final[i].area>0.0)
      if(poly_match[i][4]==1 || poly_match[i][13]==1 || poly_match[i][22]==1 || poly_match[i][31]==1)
	{
	  int count=0;
	  for(int j=0;j<Ndither*9;j++) if(poly_match[i][j]) count++;
	  fprintf(fout_overlaps,"%d %g %d :",l++,p_final[i].area,count);
	  for(int j=0;j<Ndither*9;j++) if(poly_match[i][j]) fprintf(fout_overlaps,"%d ",j);
	  fprintf(fout_overlaps,"\n");

	  for(int j=0;j<Ndither*9;j++) printf("%d ",poly_match[i][j]);
	  printf("\n");
	
	  double exterior_poly=1.0;
	  for(int j=0;j<Ndither*9;j++) if(poly_match[i][j] && (j!=4 && j!=13 && j!=22 && j!=31)) exterior_poly++;
	  stack_area[count]+=(1./exterior_poly)*p_final[i].area;
	}
  for(int i=1;i<=4;i++) printf("%d exposures over area %g / deg^2, fraction %g\n",
			      i,stack_area[i],stack_area[i]/tot_area);
  double tot_stack=0.0;
  for(int i=1;i<=4;i++) tot_stack += stack_area[i];
  printf("fraction of area covered = %g\n",tot_stack/tot_area);

  fclose(fout_overlaps);
  exit(0);
}


