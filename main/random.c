#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/************* random number gengerator (period >2e18)  ***********/

double  generate_random(long *idum)

{

const long   IM1=2147483563;
const long   IM2=2147483399;
const double AM=(1.0/IM1);
const long   IMM1=(IM1-1);
const long   IA1=40014;
const long   IA2=40692;
const long   IQ1=53668;
const long   IQ2=52774;
const long   IR1=12211;
const long   IR2=3791; 
const long   NTAB=32;
const double NDIV=(1+IMM1/NTAB);
const double EPS=1.2e-7;
const double RNMX=(1.0-EPS);

int j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[32];     /* NTAB=32 */
double temp;

 if ((*idum)<= 0) 
    {
     if(-(*idum)<1) *idum=1;
     else *idum=-(*idum);
   
     idum2 =(*idum);

     for (j=NTAB+7;j>0;j--) 
         {
          k=(*idum)/IQ1;
          *idum=IA1*(*idum-k*IQ1)-k*IR1;
          if (*idum<0) *idum +=IM1;
          if (j<NTAB) iv[j]=*idum;
          }
     iy=iv[0];
    }  

k=(*idum)/IQ1;
*idum=IA1*(*idum-k*IQ1)-k*IR1;

if (*idum<0) *idum += IM1;
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;

if (idum2<0) idum2 += IM2;

j=iy/NDIV;
iy=iv[j]-idum2;
iv[j]=*idum;

if (iy<1) iy += IMM1;

if((temp =AM*iy)> RNMX) return RNMX;

else 

return temp;

}
