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

double *rhot[5];       // Pop density at time t (five snaps needed)
double *temp;

double **vxt[5];       // x-velocity at time t
double **vyt[5];       // y-velocity at time t

double *expky;         // Array needed for the Gaussian convolution


/* Function to make space for the density array.  This is done in such a
 * way as to allow rho to be accessed either as a single block
 * or as a normal double-indexed array */

double** cart_dmalloc(int xsize, int ysize)
{
  int ix;
  double **userrho;

  userrho = malloc(xsize*sizeof(double*));
  *userrho = malloc(xsize*ysize*sizeof(double));
  for (ix=1; ix<xsize; ix++) userrho[ix] = *userrho + ix*ysize;

  return userrho;
}


/* Function to free space for the density array */

void cart_dfree(double **userrho)
{
  free(*userrho);
  free(userrho);
}


/* Function to allocate space for the global arrays */

void cart_makews(int xsize, int ysize)
{
  int s,i;

  for (s=0; s<5; s++) rhot[s] = malloc(xsize*ysize*sizeof(double));

  for (s=0; s<5; s++) {
    vxt[s] = malloc((xsize+1)*sizeof(double*));
    for (i=0; i<=xsize; i++) vxt[s][i] = malloc((ysize+1)*sizeof(double));
  }
  for (s=0; s<5; s++) {
    vyt[s] = malloc((xsize+1)*sizeof(double*));
    for (i=0; i<=xsize; i++) vyt[s][i] = malloc((ysize+1)*sizeof(double));
  }

  expky = malloc(ysize*sizeof(double));
  temp = malloc(xsize * ysize * sizeof(double));
}


/* Function to free up space for the global arrays */

void cart_freews(int xsize, int ysize)
{
  int s,i;

  for (s=0; s<5; s++) free(rhot[s]);

  for (s=0; s<5; s++) {
    for (i=0; i<=xsize; i++) free(vxt[s][i]);
    free(vxt[s]);
  }
  for (s=0; s<5; s++) {
    for (i=0; i<=xsize; i++) free(vyt[s][i]);
    free(vyt[s]);
  }

  free(expky);
  free(temp);
}


/* Copy the densities to snapshot zero. */

void cart_transform(double **userrho, int xsize, int ysize)
{
  int i;
  
  for (i=0; i < xsize * ysize; i++) {
    rhot[0][i] = (*userrho)[i];
  }
}


void cart_density_stripe(double *from, double w, double *to, int start, int step, int n)
{
  int k, i, ix;
  double a, b, c, d, e;
  double eps, sum;
  double left, right;
  
  // w = k + eps, where k is a non-negative integer and eps < 1
  k = (int) w;
  eps = w - k;
  
  sum = 0.0;
  if (k > 0) {
    sum = k - 0.5;
    for (int i=0; i<k && i<n; i++) {
      sum += from[start + i*step];
    }
    if (k > n) sum += k - n;
    sum += 0.5 * (k < n ? from[start + k*step] : 1.0);
  }
  
  for (i=0; i < n; i++) {
    ix = start + i * step;
    a = (i < k+1 ? 1.0 : from[ix - (k+1)*step]);
    b = (i < k ? 1.0 : from[ix - k*step]);
    c = (i < k-1 || i - (k-1) >= n ? 1.0 : from[ix - (k-1)*step]);
    d = (i + k < n ? from[ix + k*step] : 1.0);
    e = (i + k+1 < n ? from[ix + (k+1)*step] : 1.0);
    
    left = eps * ( (1 - eps/2) * b + (eps/2) * a );
    right = eps * ( (1 - eps/2) * d + (eps/2) * e );
    to[ix] = (left + sum + right) / (2 * w);
    
    sum += ((d + e) - (b + c)) / 2;
  }
}

/** Calculate the population density at arbitrary time, and put the result
 in a particular rhot[] snapshot array. Use a box blur with the same variance
 as the Gaussian. */
void cart_density(int from_s, double delta_t, int to_s, int xsize, int ysize)
{
  double w;
  int x, y;
  
  /* delta_t is σ^2 / 2, and the variance of a box filter with half-width w is w^2 / 3 */
  w = sqrt(6.0 * delta_t);
  
  for (x=0; x < xsize; x++) {
    cart_density_stripe(rhot[from_s], w, temp, x * ysize, 1, ysize);
  }
  for (y=0; y < ysize; y++) {
    cart_density_stripe(temp, w, rhot[to_s], y, ysize, xsize);
  }
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

  r11 = rhot[s][0];
  for (ix=1; ix<xsize; ix++) {
    r01 = r11;
    r11 = rhot[s][ix*ysize];
    vxt[s][ix][0] = -2*(r11-r01)/(r11+r01);
    vyt[s][ix][0] = 0.0;
  }

  /* Do the bottom border */

  r10 = rhot[s][ysize-1];
  for (ix=1; ix<xsize; ix++) {
    r00 = r10;
    r10 = rhot[s][ix*ysize+ysize-1];
    vxt[s][ix][ysize] = -2*(r10-r00)/(r10+r00);
    vyt[s][ix][ysize] = 0.0;
  }

  /* Left edge */

  r11 = rhot[s][0];
  for (iy=1; iy<ysize; iy++) {
    r10 = r11;
    r11 = rhot[s][iy];
    vxt[s][0][iy] = 0.0;
    vyt[s][0][iy] = -2*(r11-r10)/(r11+r10);
  }

  /* Right edge */

  r01 = rhot[s][(xsize-1)*ysize];
  for (iy=1; iy<ysize; iy++) {
    r00 = r01;
    r01 = rhot[s][(xsize-1)*ysize+iy];
    vxt[s][xsize][iy] = 0.0;
    vyt[s][xsize][iy] = -2*(r01-r00)/(r01+r00);
  }

  /* Now do all the points in the middle */

  for (ix=1; ix<xsize; ix++) {
    r01 = rhot[s][(ix-1)*ysize];
    r11 = rhot[s][ix*ysize];
    for (iy=1; iy<ysize; iy++) {
      r00 = r01;
      r10 = r11;
      r01 = rhot[s][(ix-1)*ysize+iy];
      r11 = rhot[s][ix*ysize+iy];
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
  int s0,s1,s2,s3,s4;
  int p;
  double rx1,ry1;
  double rx2,ry2;
  double rx3,ry3;
  double v1x,v1y;
  double v2x,v2y;
  double v3x,v3y;
  double v4x,v4y;
  double k1x,k1y;
  double k2x,k2y;
  double k3x,k3y;
  double k4x,k4y;
  double dx1,dy1;
  double dx2,dy2;
  double dx12,dy12;
  double dxtotal,dytotal;
  double ex,ey;
  double esq,esqmax;
  double drsq,drsqmax;

  s0 = s;
  s1 = (s+1)%5;
  s2 = (s+2)%5;
  s3 = (s+3)%5;
  s4 = (s+4)%5;

  /* Calculate the density field for the four new time slices */

  cart_density(s0, 0.5*h, s1, xsize, ysize);
  cart_density(s1, 0.5*h, s2, xsize, ysize);
  cart_density(s2, 0.5*h, s3, xsize, ysize);
  cart_density(s3, 0.5*h, s4, xsize, ysize);

  /* Calculate the resulting velocity grids */

  cart_vgrid(s1,xsize,ysize);
  cart_vgrid(s2,xsize,ysize);
  cart_vgrid(s3,xsize,ysize);
  cart_vgrid(s4,xsize,ysize);

  /* Do all three RK steps for each point in turn */

  esqmax = drsqmax = 0.0;

  for (p=0; p<npoints; p++) {

    rx1 = pointx[p];
    ry1 = pointy[p];

    /* Do the big combined (2h) RK step */

    cart_velocity(rx1,ry1,s0,xsize,ysize,&v1x,&v1y);
    k1x = 2*h*v1x;
    k1y = 2*h*v1y;
    cart_velocity(rx1+0.5*k1x,ry1+0.5*k1y,s2,xsize,ysize,&v2x,&v2y);
    k2x = 2*h*v2x;
    k2y = 2*h*v2y;
    cart_velocity(rx1+0.5*k2x,ry1+0.5*k2y,s2,xsize,ysize,&v3x,&v3y);
    k3x = 2*h*v3x;
    k3y = 2*h*v3y;
    cart_velocity(rx1+k3x,ry1+k3y,s4,xsize,ysize,&v4x,&v4y);
    k4x = 2*h*v4x;
    k4y = 2*h*v4y;

    dx12 = (k1x+k4x+2.0*(k2x+k3x))/6.0;
    dy12 = (k1y+k4y+2.0*(k2y+k3y))/6.0;

    /* Do the first small RK step.  No initial call to cart_velocity() is done
     * because it would be the same as the one above, so there's no need
     * to do it again */

    k1x = h*v1x;
    k1y = h*v1y;
    cart_velocity(rx1+0.5*k1x,ry1+0.5*k1y,s1,xsize,ysize,&v2x,&v2y);
    k2x = h*v2x;
    k2y = h*v2y;
    cart_velocity(rx1+0.5*k2x,ry1+0.5*k2y,s1,xsize,ysize,&v3x,&v3y);
    k3x = h*v3x;
    k3y = h*v3y;
    cart_velocity(rx1+k3x,ry1+k3y,s2,xsize,ysize,&v4x,&v4y);
    k4x = h*v4x;
    k4y = h*v4y;

    dx1 = (k1x+k4x+2.0*(k2x+k3x))/6.0;
    dy1 = (k1y+k4y+2.0*(k2y+k3y))/6.0;

    /* Do the second small RK step */

    rx2 = rx1 + dx1;
    ry2 = ry1 + dy1;

    cart_velocity(rx2,ry2,s2,xsize,ysize,&v1x,&v1y);
    k1x = h*v1x;
    k1y = h*v1y;
    cart_velocity(rx2+0.5*k1x,ry2+0.5*k1y,s3,xsize,ysize,&v2x,&v2y);
    k2x = h*v2x;
    k2y = h*v2y;
    cart_velocity(rx2+0.5*k2x,ry2+0.5*k2y,s3,xsize,ysize,&v3x,&v3y);
    k3x = h*v3x;
    k3y = h*v3y;
    cart_velocity(rx2+k3x,ry2+k3y,s4,xsize,ysize,&v4x,&v4y);
    k4x = h*v4x;
    k4y = h*v4y;

    dx2 = (k1x+k4x+2.0*(k2x+k3x))/6.0;
    dy2 = (k1y+k4y+2.0*(k2y+k3y))/6.0;

    /* Calculate the (squared) error */

    ex = (dx1+dx2-dx12)/15;
    ey = (dy1+dy2-dy12)/15;
    esq = ex*ex + ey*ey;
    if (esq>esqmax) esqmax = esq;

    /* Update the position of the vertex using the more accurate (two small
     * steps) result, and deal with the boundary conditions.  This code
     * does 5th-order "local extrapolation" (which just means taking
     * the estimate of the 5th-order term above and adding it to our
     * 4th-order result get a result accurate to the next highest order) */

    dxtotal = dx1 + dx2 + ex;   // Last term is local extrapolation
    dytotal = dy1 + dy2 + ey;   // Last term is local extrapolation
    drsq = dxtotal*dxtotal + dytotal*dytotal;
    if (drsq>drsqmax) drsqmax = drsq;

    rx3 = rx1 + dxtotal;
    ry3 = ry1 + dytotal;

    if (rx3<0) rx3 = 0;
    else if (rx3>xsize) rx3 = xsize;
    if (ry3<0) ry3 = 0;
    else if (ry3>ysize) ry3 = ysize;

    pointx[p] = rx3;
    pointy[p] = ry3;

  }

  *errorp = sqrt(esqmax);
  *drp =  sqrt(drsqmax);
  *spp = s4;
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

  /* Calculate the initial velocity for snapshot zero */

  cart_vgrid(0, xsize, ysize);
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
