/*
   calculates the table with the frequency-probability distribution for 
   the reemited photon packets

   The line index is also the temperature (assumes integer temperature)

   compile: cc -g -o pemit pemit.c opacity.c  ../main/inout.c ../main/bb.c  -lm

   input : [opacity.dat]
           mopacity.dat 

   output: ptable.dat (table with the reemittion frequencies)
           ftable.dat (table that has the frequencies of the reemitted photons)


   NOTE: MAY NEED TO BE CHANGED FOR DIFFERENT OPACITY 

   EFFICIENCY NOTE: The probability is too small sometimes but it is  being
                    recorded... then spend too much time trying to find
                    the correct freq to reemit photons

*/ 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../opacity/params.opa.h"
#include "../params.h"
#include "../main/structs.h"
#include "../main/functions.h"
#include "../main/opacity.h"
#include "../main/inout.h" 

double simp_kdb(double, double, double, double);

struct opacity_type opacity;
struct reemit_type reemit;

int main()

{
int i,k;
double temp,aux;
FILE *fp;

/* read the opacity table (opacity should be read in cgs units !!!) */
/* opacity table should be lambda, abs, scat, increasing lambda     */

read_opacity(&opacity); 

/* make the table of the frequencies you want to reemit the photons */

make_freq_table(&reemit);


fp = fopen("ftable.dat", "w");

for (k=0; k<NFREQ ; k++)  fprintf(fp,"%.4e\n", reemit.freq[k]);  

fclose(fp);

/* calculate the reemition probability */

/* NOTE: ASSUMES THAT THE PHOTONS WILL BE REEMITTED IN THE CHOSEN F INTERVAL ! */

aux=1.*(MAX_TEMP-1)/(TEMPDIM-1);

for(k=1;k<TEMPDIM;k++)

  {
   temp=aux*k;
   reemit.prob[k][0]=0.;
  
   for (i=1;i<NFREQ;i++)
     
       reemit.prob[k][i]=reemit.prob[k][i-1]+
                     (simp_kdb(reemit.freq[i-1],reemit.freq[i],temp,4));

   for  (i=1;i<NFREQ;i++)

      reemit.prob[k][i]=reemit.prob[k][i]/reemit.prob[k][NFREQ-1]; 
      
/* doing some diagnositcs */

   printf("\nTEMP=%4.2f\n",temp);

   for (i=0;i<NFREQ;i=i+100) printf("%.4e\n",reemit.prob[k][i]);
   
   }

/* write ptable.dat (table of reemission probabilities) */ 

fp = fopen("ptable.dat", "w");

for (i=1;i<TEMPDIM; i++)
       {
        for (k=0; k<NFREQ ; k++) 

            fprintf(fp,"%.2e ", reemit.prob[i][k]); 
     
        fprintf(fp,"\n");
       }

fclose(fp);

 
return (0);

}

/******* Function calculate integral(k dblack) (using Simpson's rule) ********/
 
double simp_kdb(a,b,temp,n) 
 
     double a; /* initial frequency */
     double b; /* final frequency   */
     double temp; 
     double n;
{ 
  double sum=0.;
  double space; 
  int counter;  
 
  
space=(b-a)/pow(2.,n);
 
sum=dblack(a,temp)*k_abs(a,opacity);  

for (counter=1 ; counter<=pow(2.,n)-2 ; counter=counter+2.) 
        
sum=sum+4.*(dblack(a+space*counter,temp)*k_abs(a+space*counter,opacity))
    +2.*(dblack(a+space*(counter+1),temp)*k_abs(a+space*(counter+1),opacity));

sum=sum
+4.*(dblack(a+space*(pow(2.,n)-1),temp)*k_abs(a+space*(pow(2.,n)-1),opacity))+
+dblack(b,temp)*k_abs(b,opacity);
 
sum=sum*(space/3.); 

return(sum);
 
} 

/*********** function make frequency table  ************************/
/* freq in cgs */

make_freq_table(reemit_ptr)

struct reemit_type *reemit_ptr;

{
int k;
double dlambda;
double l1=L1;
double l2=L2;

dlambda=log10(l2/l1)/NFREQ;

    (*reemit_ptr).freq[NFREQ-1]=l1;

    for (k=2;k<=NFREQ+1;k++)

      (*reemit_ptr).freq[NFREQ-k]=
      pow(10,(log10((*reemit_ptr).freq[NFREQ-k+1])+dlambda));
 
 
 printf("lambda2=%.4e\n",(*reemit_ptr).freq[0]);
 printf("lambda1=%.4e\n",(*reemit_ptr).freq[NFREQ-1]);

    for (k=0;k<NFREQ;k++)

     (*reemit_ptr).freq[k]=c_CONST/((*reemit_ptr).freq[k]*1.e-4);

 printf("f1=%.4e\n",(*reemit_ptr).freq[0]);
 printf("f2=%.4e\n",(*reemit_ptr).freq[NFREQ-1]);


return;

}


