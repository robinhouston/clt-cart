/* From http://www.astro.rug.nl/~gipsy/sub/bessel.c */

#include <math.h>

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10


static double bessi0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}




static double bessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}




/*
#>            bessi.dc2

Function:     bessi

Purpose:      Evaluate Modified Bessel function of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          #include "bessel.h"
              double   result; 
              result = bessi( int n,
                              double x )


              bessi    Return the Modified  Bessel function Iv(x) of 
                       integer order for input value x.
              n        Integer order of Bessel function.
              x        Double at which the function is evaluated.

                      
Description:  bessy evaluates at x the Modified Bessel function of 
              integer order n.
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/



double bessi( int n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) for n >= 0*/
/*------------------------------------------------------------*/
{
   int j;
   double bi,bim,bip,tox,ans;


   if (n < 0)
       n = -n;
   if (n == 0)
      return( bessi0(x) );
   if (n == 1)
      return( bessi1(x) );


   if (x == 0.0)
      return 0.0;
   else {
      tox=2.0/fabs(x);
      bip=ans=0.0;
      bi=1.0;
      for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
         bim=bip+j*tox*bi;
         bip=bi;
         bi=bim;
         if (fabs(bi) > BIGNO) {
            ans *= BIGNI;
            bi *= BIGNI;
            bip *= BIGNI;
         }
         if (j == n) ans=bip;
      }
      ans *= bessi0(x)/bi;
      return  x < 0.0 && n%2 == 1 ? -ans : ans;
   }
}

