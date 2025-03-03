!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: AGRMET_sfcalc_barnes
! \label{AGRMET_sfcalc_barnes}
!
! !REVISION HISTORY: 
!     28 feb 99  initial version..............capt andrus, mr moore/dnxm
!     22 jul 99  ported to IBM SP-2. added intent attributes to
!                arguments.  made all grid specific variables 
!                dynamically allocatable...................mr gayno/dnxm
!
! !INTERFACE: 
subroutine AGRMET_sfcalc_barnes( obscnt, obs, ri, rj, fguess, radius, cparam, &
     land, minwnd, isize, imax, jmax )
  
  use LIS_logMod, only : LIS_logunit
  implicit none
! !ARGUMENTS: 
  integer,       intent(in)       :: isize
  integer,       intent(in)       :: imax
  integer,       intent(in)       :: jmax
  integer,       intent(in)       :: obscnt
  real,          intent(in)       :: obs     ( isize )
  real,          intent(in)       :: ri      ( isize )
  real,          intent(in)       :: rj      ( isize )
  real,          intent(inout)    :: fguess  ( imax, jmax )
  integer,       intent(in)       :: radius  ( imax, jmax )
  character*2,   intent(in)       :: cparam
  real,          intent(in)       :: land    ( imax, jmax )  
  real,          intent(in)       :: minwnd

!
! !DESCRIPTION:   
!   barnes optimal interpolation technique subroutine
!
!   adjust the first guess with observations using the barnes
!   optimal analysis technique. \newline
!   \textbf{Method} \newline
!     - set limits for tossing bad observations (there are different
!       limits depending upon the type of data.) \newline
!     - allocate grid specific arrays. \newline
!     - initialize weight and adjustment arrays to zero. \newline
!     - loop through the observations.              \newline
!       - for good observations (not previously flagged as bad in
!         routine getsfc). \newline
!         - bilinearly interpolate the first guess value to
!           the observation location and determine the difference
!           between the observation and the interpolated first guess
!           value (the residual). \newline
!         - for residuals within the acceptable limit, set the scan
!           radius around the current observation. \newline
!           - loop through all grid points within the scan radius
!             around the observation. \newline
!             - for land points, calculate the observation weight.
!               calculate the observation's adjustment of the 
!               first guess and add it the total adjustment from
!               all observations. \newline
!     - loop through all grid points \newline
!       - if observations have adjusted the first guess, calculate
!         the final blended analysis. \newline
!       - range check the final analysis. \newline
!     - deallocate grid specific arrays. \newline
!
! The arguments and variables are:
!  \begin{description}
!  \item[obscnt]
!   Number of observations that have passed all
!   quality control checks
!  \item[obs]
!   array of observations
!  \item[ri]
!    Array of observation locations on the AGRMET grid
!   (i dimension)
!  \item[rj]
!    Array of observation locations on the AGRMET grid
!   (j dimension)
!  \item[fguess]
!   on input  - the first guess fields
!   on output - the final barnes analysis  
!  \item[radius]
!   radius of influence at each agrmet grid point
!  \item[cparam]
!    parameter to be processed (2 character code)
!    ws - wind speed
!    tp - temperature
!    rh - relative humidity
!  \item[land]
!    land/sea mask of agrmet grid
!  \item[minwnd]
!    minimum allowable windspeed on the agrmet grid   
!  \item[isize]
!   Max number of observations allowed for 
!  \item[imax]
!    east/west dimension of agrmet grid
!  \item[jmax]
!    north/south dimension of agrmet grid
!  \item[a,b,c,d]
!   variables that hold the values of the first guess
!   at the 4 agrmet grid points which surround the 
!   observation
!  \item[del]
!   sum of the product of the observation weights
!   and the residuals at an agrmet grid point
!  \item[diff]
!   array of residuals, or the difference between the
!   observation and the first guess at the
!   observation point
!  \item[diflim]
!   maximum allowable residual that is included
!   in the barnes analysis
!  \item[e,f]
!   intermediate values in the bilinear interpolation
!   of the first guess to the observation point
!  \item[maxi]
!    upper bound (i coordinate) of box in
!    which barnes analysis is performed for one
!    observation
!  \item[maxj]
!    upper bound (j coordinate) of box in
!    which barnes analysis is performed for one
!    observation
!  \item[mini]
!    lower bound (i coordinate) of box in
!    which barnes analysis is performed for one
!    observation
!  \item[minj]
!    lower bound (j coordinate) of box in
!    which barnes analysis is performed for one
!    observation
!  \item[ptrad]
!    distance between observations and agrmet 
!    grid point in grid lengths
!  \item[sr]
!    radius of influence at an agrmet grid point
!  \item[sumsqr]
!    function which calculates the sum of the square
!    of the distance between an observation and the
!    agrmet grid point in grid lengths
!  \item[wt]
!    array of observation weights
!  \item[wgt]
!    sum of the observation weights at an agrmet grid point
!  \item[x]
!    the fractional portion of variable ri
!  \item[y]
!    the fractional portion of variable rj
!   \end{description}
!
!EOP     
  integer                         :: i
  integer                         :: ii
  integer                         :: ip1
  integer                         :: is
  integer                         :: j
  integer                         :: jp1
  integer                         :: js
  integer                         :: k 
  integer                         :: maxi
  integer                         :: maxj
  integer                         :: mini
  integer                         :: minj
  real                            :: a
  real                            :: b
  real                            :: c
  real                            :: d
  real,          allocatable      :: del     ( : , : )
  real                            :: diff    ( isize )
  real                            :: diflim
  real                            :: e
  real                            :: f
  real                            :: ptrad
  real                            :: sr
  real                            :: sumsqr
  real                            :: wt      ( 0:80 )
  real,          allocatable      :: wgt     ( : , : )
  real,          allocatable      :: wgtmax  ( : , : )
  real                            :: x
  real                            :: y
  
  sumsqr(x,i,y,j) = (x-i) * (x-i) + (y-j) * (y-j)
  
  data wt &
       /1.000,0.990,0.980,0.970,0.960,0.950,0.930,0.910,0.890,0.870 &
       ,0.850,0.830,0.810,0.790,0.760,0.730,0.700,0.670,0.640,0.610 &
       ,0.580,0.550,0.520,0.490,0.460,0.430,0.400,0.380,0.360,0.340 &
       ,0.320,0.300,0.285,0.270,0.260,0.250,0.240,0.230,0.220,0.211 &
       ,0.202,0.194,0.186,0.179,0.172,0.165,0.159,0.153,0.147,0.141 &
       ,0.135,0.130,0.125,0.120,0.115,0.111,0.106,0.102,0.098,0.094 &
       ,0.091,0.087,0.084,0.081,0.077,0.074,0.071,0.069,0.066,0.063 &
       ,0.061,0.058,0.056,0.054,0.050,0.045,0.040,0.035,0.030,0.025 &
       ,0.020/  
  
!     ------------------------------------------------------------------
!     executable code starts here...set difference limit and gross error
!     limits based on which parameter is being processed. these limits
!     are selected such that every observation should be used.
!     ------------------------------------------------------------------

  if( cparam .eq. 'ws' )then
     diflim = 75.0
  elseif( cparam .eq. 'tp' )then
     diflim = 25.0 
  elseif( cparam .eq. 'rh' )then
     diflim = 0.90
  endif

!     ------------------------------------------------------------------
!     allocate local grid specific arrays.
!     ------------------------------------------------------------------
 
  allocate ( del    (imax,jmax) )
  allocate ( wgt    (imax,jmax) )
  allocate ( wgtmax (imax,jmax) )

!     ------------------------------------------------------------------
!     initialize calculation/accumulation arrays to zero.
!     ------------------------------------------------------------------

  del    = 0.0
  wgt    = 0.0
  wgtmax = 0.0

!     ------------------------------------------------------------------
!     for each observed value determine its difference from an
!     interpolated first guess value (residual).  observations that
!     failed the quality control tests in getsfc were assigned a minus 1
!     and will be ignored.
!     ------------------------------------------------------------------

  do k = 1, obscnt
     
     if (obs(k) .ge. 0.0) then
        
        i   = max(int(ri(k)),1)
        j   = max(int(rj(k)),1)
        ip1 = min((i + 1),imax)
        jp1 = min((j + 1),jmax)

!     ------------------------------------------------------------------
!         obtain an bi-linearly interpolated 1st guess value
!         and calculate the difference between observation and
!         first guess (the residual).
!     ------------------------------------------------------------------

        x = rj(k) - float(j)
        y = ri(k) - float(i)
        if((i.ge.1.and.i.le.imax).and.&
             (j.ge.1.and.j.le.jmax).and.&
             (ip1.ge.1.and.ip1.le.imax).and.&
             (jp1.ge.1.and.jp1.le.jmax)) then 
        
           a = fguess(i,j)
           b = fguess(i,jp1)
           c = fguess(ip1,j)
           d = fguess(ip1,jp1)

!     ------------------------------------------------------------------
!         for non land points, routine sfcval sets all first guess
!         values to minus 1.  if all four surrounding points are land,
!         use a bi-linear interpolation.  along coasts, the 
!         interpolation is simplified.  if there are 2 land points, 
!         perform a linear interpolation.  if there is only one land
!         point, don't interpolate, simply use first guess value 
!         at that point for calculating the residual.  If no land
!         points surround the observation, set residual to 9999 so
!         it is ignored in the analysis below.
!     ------------------------------------------------------------------

           if ((a .ge. 0.0) .and. (b .ge. 0.0) .and. &
                (c .ge. 0.0) .and. (d .ge. 0.0)) then
              e = (b - a) * x + a
              f = (d - c) * x + c
              diff(k) = obs(k) - ((f - e) * y + e)
           elseif ((a .ge. 0.0) .and. (b .ge. 0.0)) then
              diff(k) = obs(k) - ((b - a) * x + a)
           elseif ((c .ge. 0.0) .and. (d .ge. 0.0)) then
              diff(k) = obs(k) - ((d - c) * x + c)
           elseif (a .ge. 0.0) then
              diff(k) = obs(k) - a
           elseif (b .ge. 0.0) then
              diff(k) = obs(k) - b
           elseif (c .ge. 0.0) then
              diff(k) = obs(k) - c
           elseif (d .ge. 0.0) then
              diff(k) = obs(k) - d
           else
              diff(k) = 9999.0
           endif
        
!     ------------------------------------------------------------------
!         if obs-fguess difference is in reasonable range, process.
!     ------------------------------------------------------------------
        
           if( abs(diff(k)) .le. diflim )then

!     ------------------------------------------------------------------
!           calc the i- and j-coords of the "box" of points to consider
!           for the analysis.
!     ------------------------------------------------------------------
           
              sr   = float( radius( nint( ri(k) ), nint( rj(k) ) ) )
              mini = int(max(ri(k)-sr,1.0))
              minj = int(max(rj(k)-sr,1.0))
              maxi = int(min(ri(k)+sr,float(imax)))
              maxj = int(min(rj(k)+sr,float(jmax)))
           
!     ------------------------------------------------------------------
!           loop thru all land points within the "box."
!     ------------------------------------------------------------------
              
              do js = minj, maxj
                 do is = mini, maxi
                    
                    if( land(is,js) .gt. 0 )then

!     ------------------------------------------------------------------
!                 calculate the weights of the observation at each
!                 grid point based on distance.  sum the weights of
!                 all observations at a grid point.  add the 
!                 contribution of this observation to the total 
!                 adjustment and store in variable del.
!     ------------------------------------------------------------------

                       ptrad = sqrt( sumsqr(ri(k),is,rj(k),js) )
                    
                       if( ptrad.le.sr )then
                          if (sr .eq. 0) then
                             write (LIS_logunit,*) ri(k), rj(k), cparam, obs(k)
                             ii            = 0
                          else
                             ii            = int(80.0 * (ptrad / sr))
                             wgt(is,js)    = wgt(is,js) + wt(ii)
                             wgtmax(is,js) = max(wgtmax(is,js),wt(ii))
                             del(is,js)    = del(is,js) + diff(k) * wt(ii)
                          endif
                       endif
                    
                    endif
                    
                 enddo
              enddo
           
           endif
        endif
     endif
  enddo

!     ------------------------------------------------------------------
!     loop thru all grid points with positive wgtmax values.
!     ------------------------------------------------------------------

  do j = 1, jmax
     do i = 1, imax

        if( wgtmax(i,j).gt.0.0 )then

!     ------------------------------------------------------------------
!           calculate the final blended analysis.
!     ------------------------------------------------------------------
       
           fguess(i,j) = fguess(i,j) + (del(i,j) / wgt(i,j)) * &
                wgtmax(i,j)

!     ------------------------------------------------------------------
!           limit the new value to pre-determined max and min values,
!           depending on the variable.
!     ------------------------------------------------------------------

           if( cparam.eq.'ws' )then
              fguess(i,j) = min( max(fguess(i,j),minwnd), 75.0 )
           elseif( cparam.eq.'tp' )then
              fguess(i,j) = min( max(fguess(i,j),200.0), 350.0 )
           elseif( cparam.eq.'rh' )then
              fguess(i,j) = min( max(fguess(i,j),0.0), 1.0 )
            end if  
            
         endif
         
      enddo
   enddo

!     ------------------------------------------------------------------
!     deallocate local grid specific arrays.
!     ------------------------------------------------------------------
 
   deallocate ( del    )
   deallocate ( wgt    )
   deallocate ( wgtmax )
   
   return
 end subroutine AGRMET_sfcalc_barnes
 
