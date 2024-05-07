/* Created by Ziya Zhang on Nov. 2000 */
/* Compute daily PE based on monthly info */
/* Take the value on 16th day (adjusted by pe_adj) 
   of each month and interpolate between */ 
/* Time (year, month, and date) needs to be given */

/* INPUT: yr -- year
	  mon -- month
	  day -- day
	  npix -- number of hrap pixels for selected basins
	  pe -- grided PE values (1-D array) (12 months)
	  pe_adj -- grided PE adjustment values (1-D array) (12 months)
  OUTPUT: pe_day -- daily PE in 1_D array
*/

#include "models.h"

extern int days_of_month(int, int);

void get_daily_pe(int yr, int mon, int day, int npix,
		  float *pe, float *pe_adj, float *pe_day)
{

 int i, j;
 int d1, d2, d3;
 int days_tot;
 float pe_tmp, pe_adj_tmp, pe_day_tmp;
 float pe_prev, pe_this, pe_next, pe_diff;
 float pe_adj_prev, pe_adj_this, pe_adj_next;


  /*  printf("from get_daily_pe(): yr=%d, mon=%d\n",yr,mon);  */
   if (mon == 1) d1 = days_of_month(12, yr);
   else d1 = days_of_month(mon-1, yr);

   d2 = days_of_month(mon, yr);

   /** commented on 3/18/2002 
   if (mon == 12) d3 = days_of_month(1, yr);
   else d3 = days_of_month(mon+1, yr);
   **/


 /* determine value for days in the first half month */
 /* take this month's value as the value of 16th day of thiis month
    and take the 16th day of the previous month and interpolate
    between to get daily values */

 /*** DEBUG
   for (i=0; i<npix; i++) {
       printf("pe[%d]=%f7.3\n",i,pe[i]);
   }
 ****/

   if (day < 16) {

      for (i=0; i<npix; i++) {
	  if (mon == 1)  {
	      pe_prev = pe[i+11*npix];
              pe_adj_prev = pe_adj[i+11*npix];
	  }
	  else { 
	      pe_prev = pe[i+(mon-2)*npix];
              pe_adj_prev = pe_adj[i+(mon-2)*npix];
	  }
          pe_adj_this = pe_adj[i+(mon-1)*npix];
	  pe_prev = pe_prev * pe_adj_prev;
	  pe_this = pe[i+(mon-1)*npix] * pe_adj_this;
	  pe_diff = pe_this - pe_prev;
	  days_tot = d1;  /* same as days_tot = 16 + (d1 - 15) -1; */
          pe_tmp = pe_this - (16-day)*pe_diff / days_tot;
	  pe_day[i] = pe_tmp;
      }

   }  /* end of if (day < 16) */

 /* determine value for days in the second half month */
   else {          /* if (day >= 16) {  */

      for (i=0; i<npix; i++) {
          pe_adj_this = pe_adj[i+(mon-1)*npix];
	  pe_this = pe[i+(mon-1)*npix];
	  if(mon == 12) {
	      pe_next = pe[i];
              pe_adj_next = pe_adj[i];
	  }
	  else { 
	      pe_next = pe[i+mon*npix];
              pe_adj_next = pe_adj[i+mon*npix];
	  }
	  pe_this = pe_this * pe_adj_this;
	  pe_next = pe_next * pe_adj_next;
	  pe_diff = pe_this - pe_next;
	  days_tot = d2;  /* same as days_tot = d2-15+16-1; */
          pe_tmp = pe_this - (day-16)*pe_diff / days_tot;
	  pe_day[i] = pe_tmp;
      }

   }   /* end of else ---if (day >= 16) */


}  /* End of program */
