/** Compute PE demand for each time step based on given 
    monthly PE and PEADJ (may be adjusted by pec and peac) **/

/* 
   INPUT: y -- year
	  m -- month
	  d -- date
	  npix -- number of pixels of whole basin
	  dtm -- time step in minutes
	  pe -- monthly PE (array of npix)
	  pe_adj -- monthly PE Adjustment factor (array of npix)
  OUTPUT: ped -- PE demand for each time step (array of npix)
*/

#include <stdlib.h>
#include "com_header.h"
#include "models.h"

extern void get_daily_pe(int, int, int, int, float *, float *, float *);

void get_ped(int y, int m, int d, int npix, int dtm, 
	     float *pe, float *pe_adj, float *ped)
{
   int i;
   float *pe_day;

   pe_day = malloc(sizeof(float)*npix);

 /*  printf("from get_ped(): y=%d, m=%d, d=%d\n",y,m,d); */

   /* calculate daily values based on monthly info by taking 16th day of
      each month and then interpolating between linearily */

   /*** DEBUG
   for (i=0; i<npix*12; i++) {
       printf("from get_ped(): pe[%d]=%7.3f\n",i,pe[i]);
   }
    exit(1);
   ****/

   /** commented on 3/18/2002  -- take a lot of time on get_ped
   for (i=0; i<npix; i++) {
       get_daily_pe(y,m,d,npix,pe,pe_adj,pe_day);
   }
   **/
       get_daily_pe(y,m,d,npix,pe,pe_adj,pe_day);

	/*** DEBUG
     	printf("ok in get_ped() after get_daily_pe: d, pe_day[1] %d %7.3f\n",d,pe_day[1]);
	******/

   /* calculate hourly values by averaging daily values */
   for (i=0; i<npix; i++) {
       ped[i] = pe_day[i] / (24.*60./dtm);  /* distribute over a time step */
   }

  free(pe_day);


}  /* End of program */
