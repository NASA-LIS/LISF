/* Created by Ziya Zhang on Nov. 2000 */
/* Determine how many days that each month has for a given year */
/* Return number of days in integer for a given month and year */

/*
  INPUT: m -- month (integer)
	 y -- year (integer)
  OUTPUT: return number of days as integer
*/

#include <stdio.h>
#include <stdlib.h>
#include "com_header.h"
#include "models.h"

extern int is_leap_year(int);

int days_of_month(int m, int y)
{
  int num_days[12]={31,28,31,30,31,30,31,31,30,31,30,31};

     if( m <= 0 || m > 12) {
       printf("Month [1-12] and year [yyyy] need to be given correctly.\n");
       printf("Month %d, Year %d\n",m,y);
       exit(1);
     }
     if( is_leap_year(y) ) num_days[1] = 29;
     return num_days[m-1];

    /****** old version--works
	switch(m) {
	  case 1: return 31;
	  case 2:
	     if( is_leap_year(y) ) return 29;
	     else
	     return 28;
	  case 3: return 31;
	  case 4: return 30;
	  case 5: return 31;
	  case 6: return 30;
	  case 7: return 31;
	  case 8: return 31;
	  case 9: return 30;
	  case 10: return 31;
	  case 11: return 30;
	  case 12: return 31;
	  default: {
	     printf("Month [1-12] and year [yyyy] need to be given correctly.\n");
	  }
	} 
	end of old version ********/
}

