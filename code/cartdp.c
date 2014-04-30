/* Routines to transform a given set of points to a Gastner-Newman cartogram
 *
 * Written by Mark Newman
 *
 * See http://www.umich.edu/~mejn/ for further details.
 */

/* Inclusions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

#include "cart.h"

/* Constants */

#define INITH 0.001          // Initial size of a time-step
#define TARGETERROR 0.01     // Desired accuracy per step in pixels
#define MAXERROR 0.5         // Maximum permissible error
#define MAXRATIO 4.0         // Max ratio to increase step size by
#define EXPECTEDTIME 1.0e8   // Guess as to the time it will take, used to
                             // estimate completion

#define PI 3.1415926535897932384626

/* Globals */

double *rhot;          // Pop density at time t
double *fftrho;        // FT of initial density
double *fftexpt;       // FT of density at time t

double **vxt[6];       // x-velocity at time t
double **vyt[6];       // y-velocity at time t

double *expky;         // Array needed for the Gaussian convolution

fftw_plan rhotplan; // Plan for rho(t) back-transform at time t


/* Function to make space for the density array.  This is done in such a
 * way as to allow rho to be accessed either as a single block (for FFTW)
 * or as a normal double-indexed array (for regular use) */

double** cart_dmalloc(int xsize, int ysize)
{
  int ix;
  double **userrho;

  userrho = malloc(xsize*sizeof(double*));
  *userrho = fftw_malloc(xsize*ysize*sizeof(double));
  for (ix=1; ix<xsize; ix++) userrho[ix] = *userrho + ix*ysize;

  return userrho;
}


/* Function to free space for the density array */

void cart_dfree(double **userrho)
{
  fftw_free(*userrho);
  free(userrho);
}


/* Function to allocate space for the global arrays */

void cart_makews(int xsize, int ysize)
{
  int s,i;

  /* Space for the FFT arrays is allocated single blocks, rather than using
   * a true two-dimensional array, because libfftw demands that it be so */

  rhot = fftw_malloc(xsize*ysize*sizeof(double));
  fftrho = fftw_malloc(xsize*ysize*sizeof(double));
  fftexpt = fftw_malloc(xsize*ysize*sizeof(double));

  for (s=0; s<6; s++) {
    vxt[s] = malloc((xsize+1)*sizeof(double*));
    for (i=0; i<=xsize; i++) vxt[s][i] = malloc((ysize+1)*sizeof(double));
  }
  for (s=0; s<6; s++) {
    vyt[s] = malloc((xsize+1)*sizeof(double*));
    for (i=0; i<=xsize; i++) vyt[s][i] = malloc((ysize+1)*sizeof(double));
  }

  expky = malloc(ysize*sizeof(double));

  /* Make a plan for the back transform */

  rhotplan = fftw_plan_r2r_2d(xsize,ysize,fftexpt,rhot,
				   FFTW_REDFT01,FFTW_REDFT01,FFTW_MEASURE);
}


/* Function to free up space for the global arrays and destroy the FFT
 * plans */

void cart_freews(int xsize, int ysize)
{
  int s,i;

  fftw_free(rhot);
  fftw_free(fftrho);
  fftw_free(fftexpt);

  for (s=0; s<6; s++) {
    for (i=0; i<=xsize; i++) free(vxt[s][i]);
    free(vxt[s]);
  }
  for (s=0; s<6; s++) {
    for (i=0; i<=xsize; i++) free(vyt[s][i]);
    free(vyt[s]);
  }

  free(expky);

  fftw_destroy_plan(rhotplan);
}


/* Function to calculate the discrete cosine transform of the input data.
 * assumes its input is an fftw_malloced array in column-major form with
 * size xsize*ysize */

void cart_forward(double *rho, int xsize, int ysize)
{
  fftw_plan plan;

  plan = fftw_plan_r2r_2d(xsize,ysize,rho,fftrho,
                          FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}


/* Function to calculate the discrete cosine transform of the input data.
 * This function is just a wrapper for forward(), so the user doesn't
 * need to see the fftw-format density array */

void cart_transform(double **userrho, int xsize, int ysize)
{
  cart_forward(*userrho, xsize, ysize);
}


/* Function to calculate the population density at arbitrary time by back-
 * transforming and put the result in the rhot snapshot array.
 * Calculates unnormalized densities, since FFTW gives unnormalized back-
 * transforms, but this doesn't matter because the cartogram method is
 * insensitive to variation in the density by a multiplicative constant */

void cart_density(double t, int xsize, int ysize)
{
  int ix,iy;
  double kx,ky;
  double expkx;

  /* Calculate the expky array, to save time in the next part */

  for (iy=0; iy<ysize; iy++) {
    ky = PI*iy/ysize;
    expky[iy] = exp(-ky*ky*t);
  }

  /* Multiply the FT of the density by the appropriate factors */

  for (ix=0; ix<xsize; ix++) {
    kx = PI*ix/xsize;
    expkx = exp(-kx*kx*t);
    for (iy=0; iy<ysize; iy++) {
      fftexpt[ix*ysize+iy] = expkx*expky[iy]*fftrho[ix*ysize+iy];
    }
  }

  /* Perform the back-transform */

  fftw_execute(rhotplan);
}


/* Function to calculate the velocity at all integer grid points for a
 * specified snapshot */

void cart_vgrid(int s, int xsize, int ysize)
{
  int ix,iy;
  double r00,r10;
  double r01,r11;
  double mid;

  /* Do the corners */

  vxt[s][0][0] = vyt[s][0][0] = 0.0;
  vxt[s][xsize][0] = vyt[s][xsize][0] = 0.0;
  vxt[s][0][ysize] = vyt[s][0][ysize] = 0.0;
  vxt[s][xsize][ysize] = vyt[s][xsize][ysize] = 0.0;

  /* Do the top border */

  r11 = rhot[0];
  for (ix=1; ix<xsize; ix++) {
    r01 = r11;
    r11 = rhot[ix*ysize];
    vxt[s][ix][0] = -2*(r11-r01)/(r11+r01);
    vyt[s][ix][0] = 0.0;
  }

  /* Do the bottom border */

  r10 = rhot[ysize-1];
  for (ix=1; ix<xsize; ix++) {
    r00 = r10;
    r10 = rhot[ix*ysize+ysize-1];
    vxt[s][ix][ysize] = -2*(r10-r00)/(r10+r00);
    vyt[s][ix][ysize] = 0.0;
  }

  /* Left edge */

  r11 = rhot[0];
  for (iy=1; iy<ysize; iy++) {
    r10 = r11;
    r11 = rhot[iy];
    vxt[s][0][iy] = 0.0;
    vyt[s][0][iy] = -2*(r11-r10)/(r11+r10);
  }

  /* Right edge */

  r01 = rhot[(xsize-1)*ysize];
  for (iy=1; iy<ysize; iy++) {
    r00 = r01;
    r01 = rhot[(xsize-1)*ysize+iy];
    vxt[s][xsize][iy] = 0.0;
    vyt[s][xsize][iy] = -2*(r01-r00)/(r01+r00);
  }

  /* Now do all the points in the middle */

  for (ix=1; ix<xsize; ix++) {
    r01 = rhot[(ix-1)*ysize];
    r11 = rhot[ix*ysize];
    for (iy=1; iy<ysize; iy++) {
      r00 = r01;
      r10 = r11;
      r01 = rhot[(ix-1)*ysize+iy];
      r11 = rhot[ix*ysize+iy];
      mid = r10 + r00 + r11 + r01;
      vxt[s][ix][iy] = -2*(r10-r00+r11-r01)/mid;
      vyt[s][ix][iy] = -2*(r01-r00+r11-r10)/mid;
    }
  }
}



/* Function to calculate the velocity at an arbitrary point from the grid
 * velocities for a specified snapshot by interpolating between grid
 * points.  If the requested point is outside the boundaries, we
 * extrapolate (ensures smooth flow back in if we get outside by mistake,
 * although we should never actually do this because function cart_twosteps()
 * contains code to prevent it) */

void cart_velocity(double rx, double ry, int s, int xsize, int ysize,
		   double *vxp, double *vyp)
{
  int ix,iy;
  double dx,dy;
  double dx1m,dy1m;
  double w11,w21,w12,w22;

  /* Deal with the boundary conditions */

  ix = rx;
  if (ix<0) ix = 0;
  else if (ix>=xsize) ix = xsize - 1;

  iy = ry;
  if (iy<0) iy = 0;
  else if (iy>=ysize) iy = ysize - 1;

  /* Calculate the weights for the bilinear interpolation */

  dx = rx - ix;
  dy = ry - iy;

  dx1m = 1.0 - dx;
  dy1m = 1.0 - dy;

  w11 = dx1m*dy1m;
  w21 = dx*dy1m;
  w12 = dx1m*dy;
  w22 = dx*dy;

  /* Perform the interpolation for x and y components of velocity */

  *vxp = w11*vxt[s][ix][iy] + w21*vxt[s][ix+1][iy] +
         w12*vxt[s][ix][iy+1] + w22*vxt[s][ix+1][iy+1];
  *vyp = w11*vyt[s][ix][iy] + w21*vyt[s][ix+1][iy] +
         w12*vyt[s][ix][iy+1] + w22*vyt[s][ix+1][iy+1];
}


/* Function to integrate 2h time into the future two different ways using
 * four-order Runge-Kutta and compare the differences for the purposes of
 * the adaptive step size.  Parameters are:
 *   *pointx = array of x-coords of points
 *   *pointy = array of y-coords of points
 *   npoints = number of points
 *   t = current time, i.e., start time of these two steps
 *   h = delta t
 *   s = snapshot index of the initial time
 *   xsize, ysize = size of grid
 *   *errorp = the maximum integration error found for any polygon vertex for
 *             the complete two-step process
 *   *drp = maximum distance moved by any point
 *   *spp = the snapshot index for the final function evaluation
 */

void cart_twosteps(double *pointx, double *pointy, int npoints,
		   double t, double h, int s, int xsize, int ysize,
		   double *errorp, double *drp, int *spp)
{
  int s0,s1,s2,s3,s4,s5;
  int p;
  double x,y;
  double v1x,v1y;
  double v2x,v2y;
  double v3x,v3y;
  double v4x,v4y;
  double v5x,v5y;
  double v6x,v6y;
  double v7x,v7y;
  double k1x,k1y;
  double k2x,k2y;
  double k3x,k3y;
  double k4x,k4y;
  double k5x,k5y;
  double k6x,k6y;
  double k7x,k7y;
  double dx4th, dy4th;
  double dx5th, dy5th;
  double ex,ey;
  double esq,esqmax;
  double drsq,drsqmax;

  s0 = s;
  s1 = (s+1)%6;
  s2 = (s+2)%6;
  s3 = (s+3)%6;
  s4 = (s+4)%6;
  s5 = (s+5)%6;

  /* Calculate the density field for the new time slices,
     and the resulting velocity grids */

  cart_density(t+0.4*h,xsize,ysize);
  cart_vgrid(s1,xsize,ysize);

  cart_density(t+0.6*h,xsize,ysize);
  cart_vgrid(s2,xsize,ysize);

  cart_density(t+1.6*h,xsize,ysize);
  cart_vgrid(s3,xsize,ysize);

  cart_density(t+(16.0/9)*h,xsize,ysize);
  cart_vgrid(s4,xsize,ysize);

  cart_density(t+2.0*h,xsize,ysize);
  cart_vgrid(s5,xsize,ysize);

  /* Do the DOPRI steps for each point in turn */

  esqmax = drsqmax = 0.0;

  for (p=0; p<npoints; p++) {

    x = pointx[p];
    y = pointy[p];

    cart_velocity(x, y, s0, xsize, ysize, &v1x, &v1y);
    k1x = 2*h*v1x;
    k1y = 2*h*v1y;

    cart_velocity(x + 0.2*k1x, y + 0.2*k1y, s1, xsize, ysize, &v2x, &v2y);
    k2x = 2*h*v2x;
    k2y = 2*h*v2y;

    cart_velocity(x + (3.0/40)*k1x + (9.0/40)*k2x,
                  y + (3.0/40)*k1y + (9.0/40)*k2y,
                  s2, xsize, ysize, &v3x, &v3y);
    k3x = 2*h*v3x;
    k3y = 2*h*v3y;

    cart_velocity(x + (44.0/45)*k1x + (-56.0/15)*k2x + (32.0/9)*k3x,
                  y + (44.0/45)*k1y + (-56.0/15)*k2y + (32.0/9)*k3y,
                  s3, xsize, ysize, &v4x, &v4y);
    k4x = 2*h*v4x;
    k4y = 2*h*v4y;

    cart_velocity(x + (19372.0/6561)*k1x + (-25360.0/2187)*k2x + (64448.0/6561)*k3x + (-212.0/729)*k4x,
                  y + (19372.0/6561)*k1y + (-25360.0/2187)*k2y + (64448.0/6561)*k3y + (-212.0/729)*k4y,
                  s4, xsize, ysize, &v5x, &v5y);
    k5x = 2*h*v5x;
    k5y = 2*h*v5y;

    cart_velocity(x + (9017.0/3168)*k1x + (-355.0/33)*k2x + (46732.0/5247)*k3x + (49.0/176)*k4x + (-5103.0/18656)*k5x,
                  y + (9017.0/3168)*k1y + (-355.0/33)*k2y + (46732.0/5247)*k3y + (49.0/176)*k4y + (-5103.0/18656)*k5y,
                  s5, xsize, ysize, &v6x, &v6y);
    k6x = 2*h*v6x;
    k6y = 2*h*v6y;

    cart_velocity(x + (35.0/384)*k1x + (500.0/1113)*k3x + (125.0/192)*k4x + (-2187.0/6784)*k5x + (11.0/84)*k6x,
                  y + (35.0/384)*k1y + (500.0/1113)*k3y + (125.0/192)*k4y + (-2187.0/6784)*k5y + (11.0/84)*k6y,
                  s5, xsize, ysize, &v7x, &v7y);
    k7x = 2*h*v7x;
    k7y = 2*h*v7y;


    dx4th = (5179.0/57600) * k1x + (7571.0/16695) * k3x + (393.0/640) * k4x
          + (-92097.0/339200) * k5x + (187.0/2100) * k6x + (1.0/40) * k7x;
    dy4th = (5179.0/57600) * k1y + (7571.0/16695) * k3y + (393.0/640) * k4y
          + (-92097.0/339200) * k5y + (187.0/2100) * k6y + (1.0/40) * k7y;
    
    dx5th = (35.0/384) * k1x + (500.0/1113) * k3x + (125.0/192) * k4x
          + (-2187.0/6784) * k5x + (11.0/84) * k6x;
    dy5th = (35.0/384) * k1y + (500.0/1113) * k3y + (125.0/192) * k4y
          + (-2187.0/6784) * k5y + (11.0/84) * k6y;


    /* Calculate the (squared) error */

    ex = dx5th - dx4th;
    ey = dy5th - dy4th;
    esq = ex*ex + ey*ey;
    if (esq>esqmax) esqmax = esq;

    drsq = dx5th*dx5th + dy5th*dy5th;
    if (drsq>drsqmax) drsqmax = drsq;

    x += dx5th;
    y += dy5th;

    if (x<0) x = 0;
    else if (x>xsize) x = xsize;
    if (y<0) y = 0;
    else if (y>ysize) y = ysize;

    pointx[p] = x;
    pointy[p] = y;

  }

  *errorp = sqrt(esqmax);
  *drp =  sqrt(drsqmax);
  *spp = s5;
}


/* Function to estimate the percentage completion */

int cart_complete(double t)
{
  int res;

  res = 100*log(t/INITH)/log(EXPECTEDTIME/INITH);
  if (res>100) res = 100;

  return res;
}


/* Function to do the transformation of the given set of points
 * to the cartogram */

void cart_makecart(double *pointx, double *pointy, int npoints,
		   int xsize, int ysize, options_t *options)
{
  int i;
  int s,sp;
  int step;
  int done;
  double t,h,prev_h;
  double error,dr;
  double desiredratio, chosenratio;
  double *pointx_copy, *pointy_copy;

  pointx_copy = malloc(npoints * sizeof(double));
  pointy_copy = malloc(npoints * sizeof(double));

  /* Calculate the initial density and velocity for snapshot zero */

  cart_density(0.0,xsize,ysize);
  cart_vgrid(0,xsize,ysize);
  s = 0;

  /* Now integrate the points in the polygons */

  step = 0;
  t = 0.5*options->blur*options->blur;
  h = INITH;

  do {

    /* Do a combined (triple) integration step */

    memcpy(pointx_copy, pointx, npoints * sizeof(double));
    memcpy(pointy_copy, pointy, npoints * sizeof(double));
    cart_twosteps(pointx, pointy, npoints, t, h, s, xsize, ysize, &error, &dr, &sp);

    while(error > MAXERROR) {
      h /= 2;
      if (options->progress_mode == DETAILED)
        fprintf(stderr, "dr = %g and error = %g, so retrying with h = %g\n", dr, error, h);
    
      memcpy(pointx, pointx_copy, npoints * sizeof(double));
      memcpy(pointy, pointy_copy, npoints * sizeof(double));
      cart_twosteps(pointx,pointy,npoints,t,h,s,xsize,ysize,&error,&dr,&sp);
    }

    /* Increase the time by 2h and rotate snapshots */

    t += 2.0*h;
    step += 2;
    s = sp;

    /* Adjust the time-step.  Factor of 2 arises because the target for
     * the two-step process is twice the target for an individual step */

    desiredratio = pow(2 * TARGETERROR / error, 0.2);
    if (desiredratio > MAXRATIO) chosenratio = MAXRATIO;
    else chosenratio = desiredratio;
    
    prev_h = h;
    if (h * chosenratio <= options->max_h)
      h *= chosenratio;
    else
      fprintf(stderr, "h * chosenratio = %g, which is > max_h = %g\n", h * chosenratio, options->max_h);

    done = cart_complete(t);
    switch (options->progress_mode) {
      case NONE:
        break;
      case NORMAL:
        fprintf(stderr,"  %3i%%  |",done);
        for (i=0; i < done/2; i++) fprintf(stderr,"=");
        for (i=done/2; i < 50; i++) fprintf(stderr," ");
        fprintf(stderr,"|\r");
        break;
      case PERCENT:
        fprintf(stderr, "%i\n",done);
        break;
      case DETAILED:
        fprintf(stderr, "step=%d, t=%g, h=%g, dr=%g, error=%g, done=%d%%\n", step, t, h, dr, error, done);
        break;
    }
    
    /* If no point moved then we are finished */

  } while (dr > 0.0);

  switch (options->progress_mode) {
    case PERCENT:
      fprintf(stderr, "\n");
      break;
    case NORMAL:
      fprintf(stderr,"  100%%  |==================================================|\n");
      break;
    default:
      break;
  }
  
  free(pointx_copy);
  free(pointy_copy);
}
