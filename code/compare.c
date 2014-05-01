#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <fftw3.h>

#include "bessel.h"

#define PI 3.1415926535897932384626

#ifndef OFFSET
#  define OFFSET 0.005
#endif



double h_values[] = {
    0.001, 0.004, 0.016, 0.00848687, 0.00707059, 0.00996513, 0.00679333, 0.00474243,
    0.00705422, 0.0142413, 0.0143158, 0.0285653, 0.0981817, 0.248204, 0.553491,
    1.09554, 1.70166, 3.08705, 5.17658, 10.171, 22.5224, 46.0513, 88.0177, 177.542,
    357.122, 677.612, 1277.94, 2515.53, 5142.73, 10414.7, 22117.5, 50782.3, 120012,
    363376, 759655, 1.89866e+06, 7.59463e+06, 3.03785e+07, 1.21514e+08, 0
};

double *fftrho, *fftexpt, *expky, *fft_rhot, *clt_rhot, *clt_temp;
double **fftvxt, **fftvyt, **cltvxt, **cltvyt;
fftw_plan rhotplan;

double** cart_dmalloc(int xsize, int ysize)
{
    int ix;
    double **userrho;

    userrho = malloc(xsize*sizeof(double*));
    *userrho = fftw_malloc(xsize*ysize*sizeof(double));
    for (ix=1; ix<xsize; ix++) userrho[ix] = *userrho + ix*ysize;

    return userrho;
}

void cart_dfree(double **userrho)
{
    fftw_free(*userrho);
    free(userrho);
}

int readdensity(FILE *stream, double **rho, int xsize, int ysize)
{
    int ix,iy;
    int n;

    for (iy=0; iy < ysize * 3; iy++) {
        for (ix=0; ix < xsize * 3; ix++) {
            rho[ix][iy] = 1.0;
        }
    }

    for (iy=ysize; iy < ysize * 2; iy++) {
        for (ix=xsize; ix < xsize * 2; ix++) {
            n = fscanf(stream, "%lf", &rho[ix][iy]);
            if (n != 1) return 1;
            rho[ix][iy] += OFFSET;
        }
    }

    return 0;
}

void set_clt_rhot(double *rho, int xsize, int ysize) {
    int ix, iy;

    for (ix=0; ix<xsize; ix++) {
        for (iy=0; iy<ysize; iy++) {
            clt_rhot[ix*ysize + iy] = rho[ix*ysize + iy];
        }
    }
}

void cart_forward(double **rho, int xsize, int ysize)
{
    fftw_plan plan;

    set_clt_rhot(*rho, xsize, ysize);

    plan = fftw_plan_r2r_2d(xsize, ysize, *rho, fftrho,
                          FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

void cart_density_fft(double t, int xsize, int ysize)
{
    int ix, iy;
    double kx, ky;
    double expkx;

    /* Calculate the expky array, to save time in the next part */

    for (iy=0; iy<ysize; iy++) {
        ky = PI*iy/ysize;
        expky[iy] = exp(-ky*ky*t);
    }

    /* Multiply the FT of the density by the appropriate factors */

    for (ix=0; ix < xsize; ix++) {
        kx = PI*ix / xsize;
        expkx = exp(-kx*kx*t);
        for (iy=0; iy < ysize; iy++) {
            fftexpt[ix*ysize + iy] = expkx * expky[iy] * fftrho[ix*ysize + iy];
        }
    }

    /* Perform the back-transform */

    fftw_execute(rhotplan);
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

void cart_density_stripe2(double *from, double a, double b, double *to, int start, int step, int n)
{
  //double b = (1.0 - a) / 2;
  int i;
  
  for (i=1; i < n-1; i++) {
    int ix = start + i * step;
    int prev = ix - step;
    int next = ix + step;

    to[ix] = a * from[ix] + b * (from[prev] + from[next]);
  }
  to[start] = to[start + (n-1)*step] = 1.0;
}

void cart_density_clt(double delta_t, int xsize, int ysize)
{
  double w;
  int x, y;
  
  /* delta_t is σ^2 / 2, and the variance of a box filter with half-width w is w^2 / 3 */
  w = sqrt(6.0 * delta_t);
  if (w < 1) {
    double a = exp(-delta_t) * bessi(0, delta_t);
    double b = exp(-delta_t) * bessi(1, delta_t);
    printf("w = %g, a=%g / %g\n", w, a, a + 2*b);
  
    for (x=0; x < xsize; x++) {
      cart_density_stripe2(clt_rhot, a, b, clt_temp, x * ysize, 1, ysize);
    }
    for (y=0; y < ysize; y++) {
      cart_density_stripe2(clt_temp, a, b, clt_rhot, y, ysize, xsize);
    }
    return;
  }
  
  for (x=0; x < xsize; x++) {
    cart_density_stripe(clt_rhot, w, clt_temp, x * ysize, 1, ysize);
  }
  for (y=0; y < ysize; y++) {
    cart_density_stripe(clt_temp, w, clt_rhot, y, ysize, xsize);
  }
}

void cart_vgrid(int xsize, int ysize, double *rhot, double **vxt, double **vyt)
{
  int ix,iy;
  double r00,r10;
  double r01,r11;
  double mid;

  /* Do the corners */

  vxt[0][0] = vyt[0][0] = 0.0;
  vxt[xsize][0] = vyt[xsize][0] = 0.0;
  vxt[0][ysize] = vyt[0][ysize] = 0.0;
  vxt[xsize][ysize] = vyt[xsize][ysize] = 0.0;

  /* Do the top border */

  r11 = rhot[0];
  for (ix=1; ix<xsize; ix++) {
    r01 = r11;
    r11 = rhot[ix*ysize];
    vxt[ix][0] = -2*(r11-r01)/(r11+r01);
    vyt[ix][0] = 0.0;
  }

  /* Do the bottom border */

  r10 = rhot[ysize-1];
  for (ix=1; ix<xsize; ix++) {
    r00 = r10;
    r10 = rhot[ix*ysize+ysize-1];
    vxt[ix][ysize] = -2*(r10-r00)/(r10+r00);
    vyt[ix][ysize] = 0.0;
  }

  /* Left edge */

  r11 = rhot[0];
  for (iy=1; iy<ysize; iy++) {
    r10 = r11;
    r11 = rhot[iy];
    vxt[0][iy] = 0.0;
    vyt[0][iy] = -2*(r11-r10)/(r11+r10);
  }

  /* Right edge */

  r01 = rhot[(xsize-1)*ysize];
  for (iy=1; iy<ysize; iy++) {
    r00 = r01;
    r01 = rhot[(xsize-1)*ysize+iy];
    vxt[xsize][iy] = 0.0;
    vyt[xsize][iy] = -2*(r01-r00)/(r01+r00);
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
      vxt[ix][iy] = -2*(r10-r00+r11-r01)/mid;
      vyt[ix][iy] = -2*(r01-r00+r11-r10)/mid;
    }
  }
}

double compare(int xsize, int ysize) {
    int ix, iy;
    double max_sq_err;

    cart_vgrid(xsize, ysize, fft_rhot, fftvxt, fftvyt);
    cart_vgrid(xsize, ysize, clt_rhot, cltvxt, cltvyt);

    max_sq_err = 0;
    for (ix=0; ix<=xsize; ix++) {
        for (iy=0; iy<=ysize; iy++) {
            double dx = fftvxt[ix][iy] - cltvxt[ix][iy];
            double dy = fftvyt[ix][iy] - cltvyt[ix][iy];

            double sq_err = dx*dx + dy*dy;
            if (sq_err > max_sq_err) max_sq_err = sq_err;
        }
    }

    return sqrt(max_sq_err);
}


int main(int argc, char **argv) {
    int xsize, ysize;
    double **rho_0;
    FILE *infp;
    int i;
    double t;
    double cumulative_error;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s xsize ysize density-grid\n", argv[0]);
        exit(1);
    }

    xsize = atoi(argv[1]);
    ysize = atoi(argv[2]);

    infp = fopen(argv[3], "r");
    if (infp == NULL) {
        fprintf(stderr, "%s: unable to open file '%s': ", argv[0], argv[3]);
        perror(NULL);
        exit(4);
    }

    /* Read in the density data and transform it */
    rho_0 = cart_dmalloc(xsize * 3, ysize * 3);
    if (readdensity(infp, rho_0, xsize, ysize)) {
        fprintf(stderr, "%s: density file contains too few or incorrect data\n", argv[0]);
        exit(6);
    }
    fclose(infp);

    fftrho = fftw_malloc(xsize*3 * ysize*3 * sizeof(double));
    clt_rhot = malloc(xsize*3 * ysize*3 * sizeof(double));
    cart_forward(rho_0, xsize*3, ysize*3);

    fftexpt = fftw_malloc(xsize*3 * ysize*3 * sizeof(double));
    fft_rhot = fftw_malloc(xsize*3 * ysize*3 * sizeof(double));
    rhotplan = fftw_plan_r2r_2d(xsize*3, ysize*3, fftexpt, fft_rhot,
                                FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);
    expky = malloc(ysize*3 * sizeof(double));
    clt_temp = malloc(xsize*3 * ysize*3 * sizeof(double));

    fftvxt = malloc((xsize*3 + 1) * sizeof(double*));
    cltvxt = malloc((xsize*3 + 1) * sizeof(double*));
    fftvyt = malloc((xsize*3 + 1) * sizeof(double*));
    cltvyt = malloc((xsize*3 + 1) * sizeof(double*));
    for (i=0; i<=xsize*3; i++) {
        fftvxt[i] = malloc((ysize*3 + 1) * sizeof(double));
        cltvxt[i] = malloc((ysize*3 + 1) * sizeof(double));
        fftvyt[i] = malloc((ysize*3 + 1) * sizeof(double));
        cltvyt[i] = malloc((ysize*3 + 1) * sizeof(double));
    }

    t = 0;
    cumulative_error = 0;

    for (i=0; h_values[i] != 0; i++) {
        double h, err;
        
        h = h_values[i];
        t += 2*h;

        cart_density_fft(t, xsize*3, ysize*3);

        if (t < 1) {
            set_clt_rhot(*rho_0, xsize*3, ysize*3);
            cart_density_clt(0.5*t, xsize*3, ysize*3);
            cart_density_clt(0.5*t, xsize*3, ysize*3);
            cart_density_clt(0.5*t, xsize*3, ysize*3);
            cart_density_clt(0.5*t, xsize*3, ysize*3);
            // cart_density_clt(2.0*t, xsize*3, ysize*3);
        } else {
            cart_density_clt(0.5*h, xsize*3, ysize*3);
            cart_density_clt(0.5*h, xsize*3, ysize*3);
            cart_density_clt(0.5*h, xsize*3, ysize*3);
            cart_density_clt(0.5*h, xsize*3, ysize*3);
        }

        err = compare(xsize*3, ysize*3);
        cumulative_error += h*err;
        printf("[%d] t=%g, h=%g, error=%g\n", i, t, h, h*err);
    }
    printf("\nCumulative error = %g\n", cumulative_error);

    cart_dfree(rho_0);
    free(clt_rhot);
    free(clt_temp);

    for (i=0; i<=xsize; i++) {
        free(fftvxt[i]);
        free(cltvxt[i]);
    }
    free(fftvxt);
    free(cltvxt);

    for (i=0; i<=xsize; i++) {
        free(fftvyt[i]);
        free(cltvyt[i]);
    }
    free(fftvyt);
    free(cltvyt);

    free(expky);
    fftw_free(fftrho);
    fftw_free(fftexpt);
    fftw_free(fft_rhot);
    fftw_destroy_plan(rhotplan);
    return 0;
}
