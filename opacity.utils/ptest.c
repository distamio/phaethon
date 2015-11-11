/* 
  Basically reads in ptable.dat and outputs and sm friendly ptable_test.dat 
  (columns of reemission probabilities for selected temperatures 
   100K,200K,..and frequencies).

  use this to plot cumulative probability

  compile: cc -g -o ptest ptest.c opacity.c ../main/inout.c ../main/bb.c -lm

  input: ptable.dat [in ../ directory !]
         ftable.dat
     
  output:ptable_test.dat
         ftalbe_test.dat
*/

#define MAX_LINE     400

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "structs.h"
#include "params.opa.h"
#include "params.h"
#include "phaethon.params.h"
#include "constants.h"
#include "functions.h"
#include "inout.h"
#include "opacity.h"


void write_ptable2(struct reemit_type);

struct reemit_type reemit;

int main()

{
int i,k;
FILE *fp;

/* read the table that has the reemission frequencies */

read_ftable(&reemit);

/* read the table that has the reemission probabilities */

read_ptable2(&reemit);

 
   if (( fp = fopen("ptable_test.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN ptable_test.dat FILE. \n");
         exit(0);
         }

   /*  Note: the first row in ptable.dat is for T=1 */

   for (k=0; k<NFREQ ; k++)
       {
        for  (i=1;i<200; i=i+10)

            fprintf(fp,"%.4e ", reemit.prob[i][k]); 
     
        fprintf(fp,"\n");
       }

   fclose(fp);


  
   if (( fp = fopen("ftable_test.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN ftable_test.dat FILE. \n");
         exit(0);
         }

   for (k=0; k<NFREQ ; k=k++)
       {
            fprintf(fp,"%.4e\n ",reemit.freq[k]); 
       }

   fclose(fp);


return 0;

}


