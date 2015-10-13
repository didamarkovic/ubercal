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

const int NRA = 2;   // number of exposures in x - direction
const int NDEC = 2;   // number of exposures in y - direction
const int NDITH = 4;

struct poly {
  double area;
  double corners[4][2];
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

int find_rot(double[2][3],double[2][3],double[3][3]);

// the six numbers that give the dither strategy
// const double dither[6]={0.0,0.0,0.0,0.0,0.0,0.0};
const double dither[6]={50.0,100.0,0.0,100.0,0.0,100.0}; // J-pattern
// const double dither[6]={50.0,100.0,50.0,0.0,50.0,0.0}; // step-pattern
// const double dither[6]={50.0,100.0,0.0,100.0,50.0,100.0}; // S-pattern
// const double dither[6]={50.0, 100.0, 50.0, 0.0, 50.0, 100.0}; // N-pattern
// const double dither[6]={0.0, 0.0, 318.0, 331.0, 0.0, 0.0}; // X-pattern
// const double dither[6]={50.0,0.0,0.0,100.0,-50.0,0.0}; // box-pattern
// const double dither[6]={275.0,0.0,0.0,350.0,-275.0,0.0};
// const double dither[6]={1599.0,0.0,0.0,1774.0,-1599.0,0.0};

int main() {

  const int bsz=800; char buf[bsz], fname[bsz];

  // ********************************************************************************
  // memory for polygons - this is overkill for final ones, as we don't store vertices
  long NPOLY_MAX = 10000000;
  struct poly *p_init;
  if(!(p_init = (struct poly*)malloc(NPOLY_MAX*sizeof(struct poly)))) 
    err_known("memory allocation problem for polygons");

  struct poly *p_final;
  if(!(p_final = (struct poly*)malloc(NPOLY_MAX*sizeof(struct poly)))) 
    err_known("memory allocation problem for polygons");
  for(int i=0;i<NPOLY_MAX;i++) {
    p_final[i].nexposure=0;
    for(int j=0;j<NDITH;j++) p_final[i].iexposure[j]=0;
  }
  
  // ********************************************************************************
  // need to work out relative positions of all CCDs in relevent exposures
  // this is where we put in the dither information - number, offsets

  // set up dither offsets
  double xdith[NDITH], ydith[NDITH];
  for(int idither=0;idither<NDITH;idither++) {
    if(idither==0) { xdith[idither]=0.0; ydith[idither]=0.0; }
    if(idither==1) { xdith[idither]=xdith[idither-1]+dither[0]; ydith[idither]=ydith[idither-1]+dither[1]; }
    if(idither==2) { xdith[idither]=xdith[idither-1]+dither[2]; ydith[idither]=ydith[idither-1]+dither[3]; }
    if(idither==3) { xdith[idither]=xdith[idither-1]+dither[4]; ydith[idither]=ydith[idither-1]+dither[5]; }
  }

  // shifts in x and y marking CCD edges
  double xp_shift = 4.0*(npix + xgap_size/pix_size)*pix_scale;
  double yp_shift = 4.0*(npix + ygap_size/pix_size)*pix_scale;

  int ipoly=0, iexpos=0, idetect=0;    
  for(int idither=0;idither<NDITH;idither++) {
       
    for(int jpoint=0;jpoint<NDEC;jpoint++) {
      
      double ypoint = (double)jpoint * yp_shift;
      
      for(int ipoint=0;ipoint<NRA;ipoint++) {
  
  double xpoint = (double)ipoint * xp_shift;
  
  // loop over CCDs in each exposure - polynomial for each
  for(int i=0;i<4;i++)      
    for(int j=0;j<4;j++) {

      // coords of polygon in ra, (dec+pi/2) / arcsec
      double poly_ralow = xstart + xpoint + xdith[idither] + (double)i*(npix + xgap_size/pix_size)*pix_scale;
      double poly_rahig = poly_ralow + npix*pix_scale;
      double poly_delow = ypoint + ydith[idither] + (double)j*(npix + ygap_size/pix_size)*pix_scale;
      double poly_dehig = poly_delow + npix*pix_scale;

      p_init[ipoly].corners[0][1] = poly_dehig;
      p_init[ipoly].corners[0][0] = poly_ralow;
      p_init[ipoly].corners[1][1] = poly_delow;
      p_init[ipoly].corners[1][0] = poly_ralow;
      p_init[ipoly].corners[2][1] = poly_delow;
      p_init[ipoly].corners[2][0] = poly_rahig;
      p_init[ipoly].corners[3][1] = poly_dehig;
      p_init[ipoly].corners[3][0] = poly_rahig;
      p_init[ipoly].area          = 0.0;
      p_init[ipoly].nexposure     = 1;
      p_init[ipoly].iexposure[0]  = idetect;      
      // p_init[ipoly].iexposure[0]  = iexpos;      
      if(++ipoly>NPOLY_MAX) err_known("too many initial polynomials for memory");

      idetect++;
    }

  iexpos++;
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
    fprintf(fout_mask,"vertices 0 ( 4 vertices, 1 weight, %lf %lf mid):\n",
      0.5*(p_init[ipoly].corners[1][0]+p_init[ipoly].corners[3][0]),
      0.5*(p_init[ipoly].corners[1][1]+p_init[ipoly].corners[3][1])); // arcsec
    fprintf(fout_mask,"  %lf %lf\t%lf %lf\t%lf %lf\t%lf %lf\n",
      p_init[ipoly].corners[0][0],p_init[ipoly].corners[0][1],
      p_init[ipoly].corners[1][0],p_init[ipoly].corners[1][1],
      p_init[ipoly].corners[2][0],p_init[ipoly].corners[2][1],
      p_init[ipoly].corners[3][0],p_init[ipoly].corners[3][1]);
  }
  fclose(fout_mask);

  // exit(0);
  
  // ********************************************************************************

  // now run mangle to Balkanise mask
  // doesn't work without the pixelization step - ???
  system("./mangle2.2/bin/pixelize full-pointing.vrt full-pointing.pol");
  system("./mangle2.2/bin/snap full-pointing.pol full-pointing.pol");
  system("./mangle2.2/bin/balkanize -Ba -ovs full-pointing.pol full-pointing.vrt");
  
  // ********************************************************************************
  // now read in, and sum over multiple regions
  FILE *fp;
  
  if((fp=fopen(fname,"r"))==NULL)
    err_known("cannot open input file");

  // read in total number of polygons
  long npoly_final=0;
  double tot_area=0.0;
  fscanf(fp,"%ld polygons\n",&npoly_final);
  printf("%ld polygons in balkanised file\n",npoly_final);
  if(npoly_final>NPOLY_MAX) err_known("too many polygons in file");
  fgets(buf,bsz,fp);

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
  if((p_init[j].corners[3][0] < (x[iv]-XBIT)) ||
     (p_init[j].corners[1][0] > (x[iv]+XBIT)) ||
     (p_init[j].corners[3][1] < (y[iv]-XBIT)) || 
     (p_init[j].corners[1][1] > (y[iv]+XBIT)) ) parent=0;
      if(parent) {
  p_final[ipoly].iexposure[p_final[ipoly].nexposure]=p_init[j].iexposure[0];  
  if( (p_final[ipoly].nexposure++) > NDITH) err_known("poly match");
      }
    }
    
    // calculate area in deg^2
    p_final[ipoly].area = 0.0;
    int j=nvert-1;
    for(int iv=0;iv<nvert;iv++) {
      p_final[ipoly].area += (x[iv]+x[j]) * (y[iv]-y[j]);
      j=iv;
    }
    p_final[ipoly].area *= 0.5/(3600.*3600.);

    tot_area += p_final[ipoly].area;
    
    // output polygons in mask
    //fprintf(fout_mask,"vertices 0 ( %d vertices, %d weight, 1 1 mid):\n",nvert,p_final[ipoly].nexposure);
    //for(int iv=0;iv<nvert;iv++) fprintf(fout_mask,"%g %g ",x[iv],y[iv]);
    //fprintf(fout_mask,"\n");
    
  }
  fclose(fp);
  //fclose(fout_mask);


  printf("tot area = %g\n",tot_area);
  // exit(0);

  /*
  for(int i=0;i<npoly_final;i++) {
    printf("%d %d :",i,p_final[i].nexposure);
    for(int ie=0;ie<p_final[i].nexposure;ie++) printf("%d ",p_final[i].iexposure[ie]);
    printf("\n");
  }
  */
  
  // ********************************************************************************
  // sort final polygons to reduce the number to unique sets of overlaps
  // after this the positions of the vertices become unimportant
  for(int i=npoly_final;i>=0;i--)
    for(int j=0;j<i;j++) {
      int match=1;
      for(int k=0;k<NDITH;k++) if(p_final[i].iexposure[k]!=p_final[j].iexposure[k]) { match=0; break; }
      if(match) {
  p_final[j].area+=p_final[i].area;
  p_final[i].area=0.0;
  break;
      }
    }

  // find number of overlaps
  int noverlap=0;
  for(int i=0;i<npoly_final;i++) if(p_final[i].area>0.0) noverlap++;
  printf("%d distinct overlap regions\n",noverlap);

  // ********************************************************************************
  // print out final list of overlaps for "full survey"
  FILE *fout_survey;
  sprintf(fname,"full-survey-overlaps.txt");
  if((fout_survey=fopen(fname,"w"))==NULL)
    err_known("cannot open output file");

  fprintf(fout_survey,"# %d dithers\n",NDITH);
  fprintf(fout_survey,"# NRA=%d, NDEC=%d\n",NRA,NDEC);
  fprintf(fout_survey,"# %d overlap polygons\n",noverlap);
  for(int i=0;i<npoly_final;i++) if(p_final[i].area>0.0) {
    fprintf(fout_survey,"%d %g %d :",i,p_final[i].area,p_final[i].nexposure);
    printf("%d %g %d :",i,p_final[i].area,p_final[i].nexposure);
    for(int ie=0;ie<p_final[i].nexposure;ie++) {
      fprintf(fout_survey,"%d ",p_final[i].iexposure[ie]);
      printf("%d ",p_final[i].iexposure[ie]);
    }
    printf("\n");
    fprintf(fout_survey,"\n");
  }
  fclose(fout_survey);

  exit(0);
}

int find_rot(double a[2][3], double b[2][3], double rot[3][3]) {

  int ncross_product(double[3],double[3],double[3]);

  double aa[3][3], bb[3][3];

  // find 3 orthogonal vectors of original pointing
  for(int i=0;i<3;i++) aa[0][i] = a[0][i];
  ncross_product(aa[0],a[1],aa[1]);
  ncross_product(aa[0],aa[1],aa[2]);

  // for(int i=0;i<3;i++) for(int j=0;j<3;j++) printf("aa[%d][%d]=%g\n",i,j,aa[i][j]);
  
  // find 3 orthogonal vectors of the rotated pointing
  for(int i=0;i<3;i++) bb[0][i] = b[0][i];
  ncross_product(bb[0],b[1],bb[1]);
  ncross_product(bb[0],bb[1],bb[2]);

  // for(int i=0;i<3;i++) for(int j=0;j<3;j++) printf("bb[%d][%d]=%g\n",i,j,bb[i][j]);
  
  // From Triad method, http://en.wikipedia.org/wiki/Triad_method, 
  // matrix is (aa[0]:aa[1];aa[2])(bb[0]:bb[1]:bb[2])^T
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++) {
      rot[i][j] = 0.0;
      for(int k=0;k<3;k++) rot[i][j] += aa[k][i]*bb[k][j];
    }
  
  return 1;
}

int ncross_product(double a[3], double b[3], double c[3]) {

  // normalised cross product of 2 vectors
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]-b[0];

  double c_amp = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
  for(int i=0;i<3;i++) c[i] /= c_amp;
  
  return 1;
}
