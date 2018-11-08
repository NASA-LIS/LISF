c Modified Ziya Zhang's C code for Greenness on May 2009
c Compute daily GRND based on monthly info 
c/* Take the value on 16th day of each month and interpolate between */ 
c/* Time (year, month, and date) needs to be given */
C
c INPUT: yr -- year
c	  mon -- month
c	  day -- day
c	  npix -- number of hrap pixels for selected basins
c	  grn -- grided grn values (1-D array) (12 months)
c  OUTPUT: daily grnd 
C
c#include "models.h"

      subroutine get_daily_grn(yr,mon,day,npix,grn,grnday)      


      integer yr,mon,day
      integer d1, d2, d3
      integer days_tot
      integer ndays(12)/31,28,31,30,31,30,31,31,30,31,30,31/
      real grn(12),grnday(1)
      real grn_tmp, grn_prev, grn_this, grn_next, grn_diff

      if (mon .eq. 1) then
       d1 = ndays(mon)
       if(mod(yr,4) .eq. 0 .and. mon .eq. 2) d1=d1+1
      else 
       d1 = ndays(mon-1)
       if(mod(yr,4) .eq. 0 .and. mon-1 .eq. 2) d1=d1+1
      endif
      d2 = ndays(mon)
      if(mod(yr,4) .eq. 0 .and. mon .eq. 2) d2=d2+1

c /* determine value for days in the first half month */
c /* take this month's value as the value of 16th day of this month
c    and take the 16th day of the previous month and interpolate
c    between to get daily values */

      if (day .lt. 16) then
       do i=1,npix 
	if (mon .eq. 1) then 
	 grn_prev = grn(i+11*npix)
        else 
	 grn_prev = grn(i+(mon-2)*npix)
	endif
        grn_this = grn(i+(mon-1)*npix)
	grn_diff = grn_this - grn_prev
	days_tot = d1  
        grn_tmp = grn_this - (16-day)*grn_diff / days_tot
	grnday(i) = grn_tmp;
       enddo
c/* end of if (day < 16) */

c /* determine value for days in the second half month */
      else 
c            /* if (day >= 16) {  */
       do i=1,npix
	grn_this = grn(i+(mon-1)*npix)
	if(mon .eq. 12) then
	 grn_next = grn(i)
	else
	 grn_next = grn(i+mon*npix)
	endif 
	grn_diff = grn_this - grn_next
	days_tot = d2;  
        grn_tmp = grn_this - (day-16)*grn_diff / days_tot
	grnday(i) = grn_tmp
       enddo
      endif
c   /* end of else ---if (day >= 16) */

      return
      end
      

