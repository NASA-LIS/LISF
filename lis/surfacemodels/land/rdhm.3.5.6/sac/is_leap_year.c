/* Created by Ziya Zhang on Nov. 2000 */
/* Check if a given year is a leap year or not */
/* Returns true (1) if it is, return false (0) otherwise */

/*
  INPUT:  yr -- year (integer)
  OUTPUT: return: 1 (true), 0 (false)
*/

#include "com_header.h"
#include "models.h"

int is_leap_year(int yr)
{
	if ( ((yr % 4 == 0) && (yr % 100 != 0)) || (yr % 400 == 0) )
	   return 1;  /* is a leap year */
	else
	   return 0;  /* is NOT a leap year */
}
