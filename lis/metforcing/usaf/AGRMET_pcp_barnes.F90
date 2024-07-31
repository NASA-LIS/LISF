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
! !ROUTINE: AGRMET_pcp_barnes
! \label{AGRMET_pcp_barnes}
!
! !REVISION HISTORY:
!
!
!     15 jun 96  initial version..................mr moore sysm(agromet)
!     10 apr 97  brought up to software standards.  combined barn1, 
!                barn2, and blend into one routine.  made use of source 
!                type by using different spread weights and radius based
!                upon source type.......................ssgt miller/sysm
!     02 may 97  changed order of array indices to prevent 'thrashing'  
!                in memory. updated prolog and brought up to standards. 
!                ..............................capt andrus/sysm(agromet)
!      7 oct 99  ported to ibm sp-2, updated prolog, incorporated 
!                FORTRAN 90 features.................capt hidalgo/agrmet
!     31 aug 01  corrected error which assigned wrong weights and
!                radii to the wrong precipitation source.  geoprecip
!                and climo were reversed as were pres/past wx and
!                ssm/i.  added code to use addrad(7) and (8) which were
!                not being used.  made corresponding fixes to 
!                control.spread file.......................mr gayno/dnxm
!     3 nov 05 Sujay Kumar, Initial Code
!     30 AUG 2010 Modified to include observations that, due to land mask
!                 resolution limits, were not included as they fell over 
!                 grid cells flagged as "water" in coarser resolutions
!                 ......................................Michael Shaw/WXE
!     13 May 2013 Modified to not spread obs if background is CMORPH
!                 ..............................Ryan Ruhge/16WS/WXE/SEMS
!
! !INTERFACE:    
subroutine AGRMET_pcp_barnes( n, mrgp, radius, addrad, maxrad, minrad, &
     src, srcwts, varrad, tmp, wgt, cmorphpixel, imax, jmax)
! !USES: 
  use LIS_coreMod, only : LIS_rc
  use LIS_LMLCMod, only : LIS_LMLC
  use LIS_logMod,  only : LIS_logunit

  implicit none

! !ARGUMENTS: 
  integer,    intent(in)        :: n 
  integer,    intent(in)        :: addrad(8) 
  integer,    intent(in)        :: maxrad(6)
  integer,    intent(in)        :: minrad(6)  
  real,       intent(in)        :: radius(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer,    intent(in)        :: src(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer,    intent(in)        :: varrad(4) 
  real,      intent(inout)      :: mrgp(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,      intent(in)         :: srcwts(8) 
  real,      intent(out)        :: tmp(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,      intent(inout)      :: wgt(LIS_rc%lnc(n),LIS_rc%lnr(n))
  logical,   intent(in)         :: cmorphpixel(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer,   intent(in)         :: imax
  integer,   intent(in)         :: jmax
!
! !DESCRIPTION:
!
!     to perform a barnes analysis on the precipitation data.
!
!     \textbf{Method} \newline
!
!     1. loop thru pts in the hemisphere \newline
!       a. if pt's parsed real precip amt is valid ... \newline
!         1) set a local scan area \newline
!         2) calc a local mean parsed real precip amt \newline
!         3) eval amt diff btwn local mean and current pt
!            adjust initial spread radius accordingl
!       b. set final spreading radius for each precip source \newline
!     2. initialize variables to zero \newline
!     3. loop thru pts in the hemisphere \newline
!       a. set the spreading radius in grid dist's. \newline
!       b. determine the surrounding pts to be affected by
!          the amt at the pt (the spread pt). \newline
!       c. for each valid surrounding pt... \newline
!         1) calc the pt's dist from the spread pt. \newline
!         2) if the pt is within the spread radius, calc
!            and accum the spread wgt and a weighted 3-hrly
!            merged precip amts. \newline
!     4. again, loop thru the pts in the hemisphere... \newline
!        if a precip amt was spread to this pt, calc a new
!        3-hrly merged precip value for the pt using the
!        accumulated spread wgt and the accum weighted
!        3-hrly merged precip amt. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[a]         dummy variable  
!   \item[addrad]    added radius value for barnes analysis  
!   \item[b]         dummy variable  
!   \item[cmorphpixel] indicates if pixel is fed with cmorph data
!   \item[cnt]       cnumber of points to which values were spread 
!   \item[idis]      integer distance from pt to surrounding pts   
!   \item[iend]      end point for ii loop   
!   \item[ii]        loop index for local spread area
!   \item[imax]      number of gridpoints in east-west direction
!   \item[irad]      radius of the observation points
!   \item[is]        loop index for supergrid
!   \item[istrt]     start point for ii loop 
!   \item[jend]      end point for jj loop   
!   \item[jj]        loop index for local spread area
!   \item[jmax]      number of gridpoints in north-south direction
!   \item[js]        loop index for supergrid
!   \item[jstrt]     start point for jj loop 
!   \item[land]      supergrid point processing switches 
!   \item[maxi]      upper i coordinate value
!   \item[maxj]      upper j coordinate value
!   \item[maxrad]    maximum radius value
!   \item[meanv]     mean value of other (local) precip amts
!   \item[mini]      lower i coordinate value
!   \item[minj]      lower j coordinate value
!   \item[minrad]    minimum radius value
!   \item[mrgp]      merged precip amount
!   \item[numval]    counter 
!   \item[rad]       distance from (i,j) to (is,js)  
!   \item[radius]    supergrid integer radius for obs spreading  
!   \item[ri]        real version of i   
!   \item[rj]        real version of j   
!   \item[rrad]      supergrid real scan radius for obs spreading
!   \item[rsqd]      distance squared from (i,j) to (is,js)  
!   \item[sigmav]    std deviation of surrounding obs values 
!   \item[src]       source amount gridpoints
!   \item[srcdis]    source distance of the points  
!   \item[srcwt]     weighted source 
!   \item[srcwts]    value weight placed on each source  
!   \item[sumsqr]    sum of squares statement function   
!   \item[sumval]    accumulation variable   
!   \item[tmp]       temporary array to hold weighted values 
!   \item[varrad]    radius change amount due to the real precip variance
!   \item[vdiff]     abs value of diff btwn srel and meanv   
!   \item[wgt]       weighted array point values 
!   \item[wt]        1.0 to .018 weight applied to i/j   
!   \item[x]         dummy variable  
!   \item[y]         dummy variable  
!   \item[ztest]     variable that holds proper spread weight
!  \end{description}
!EOP
  integer                       :: cnt   
  integer                       :: idis  
  integer                       :: iend  
  integer                       :: ii
  integer                       :: irad  
  integer                       :: is
  integer                       :: istrt 
  integer                       :: jend  
  integer                       :: jj
  integer                       :: js
  integer                       :: jstrt 
  integer                       :: maxi  
  integer                       :: maxj  
  integer                       :: mini  
  integer                       :: minj  
  integer                       :: numval
  integer                       :: ztest 
  real                          :: a
  real                          :: b
  real                          :: meanv 
  real                          :: rad   
  real                          :: ri
  real                          :: rj
  real                          :: rrad  
  real                          :: rsqd  
  real                          :: sigmav
  real                          :: srcdis
  real                          :: srcwt 
  real                          :: sumsqr
  real                          :: sumval
  real                          :: vdiff 
  real                          :: wt(0:20)  
  real                          :: x 
  real                          :: y   
  data wt / 1.00, 0.98, 0.96, 0.96, 0.94,   &
       0.92, 0.90, 0.88, 0.85, 0.80, 0.70, &
       0.60, 0.50, 0.40, 0.30, 0.15,  & 
       0.10, 0.08, 0.06, 0.04, 0.02 /  

!     ------------------------------------------------------------------
!     define an internal function, sumsqr (the sum of the squares)
!     ------------------------------------------------------------------

  sumsqr(x, a, y, b) = ((x - a) * (x - a)) + ((y - b) * (y - b))   

!     ------------------------------------------------------------------
!     executable code begins here ... loop thru points in the hemisphere
!     ------------------------------------------------------------------

  cnt = 0
  
  do js = 1, LIS_rc%lnr(n)
     do is = 1, LIS_rc%lnc(n)

! ----------------------------------------------------------------------
!       Don't spread obs if 'background' source is CMORPH.
! ----------------------------------------------------------------------
        if (cmorphpixel(is,js) .eqv. .false.) then

           if( mrgp(is,js) .gt. -9998.0 )then

!     ------------------------------------------------------------------
!           precip amount is valid.
!           set the proper spread weight. if the value being spread is  
!           a real zero amount or a bogus zero amount, it gets a
!           different srcwts value. 
!     ------------------------------------------------------------------

              ztest = nint( mrgp(is,js) * 1000 )

              if( ztest .eq. 0 )then
                 if( src(is,js) .eq. 1 )then
                    srcwt = srcwts(7)
                 elseif( src(is,js) .eq. 3 )then
                    srcwt = srcwts(8)
                 else
                    srcwt = srcwts(src(is,js))
                    if( src(is,js) .le. 0 ) srcwt = 0.0
                    if( src(is,js) .gt. 6 ) srcwt = 0.0
                 endif
              else
                 srcwt = srcwts(src(is,js))
                 if( src(is,js) .lt. 1 ) srcwt = 0.0
                 if( src(is,js) .gt. 6 ) srcwt = 0.0
              endif
           
              if( srcwt .gt. 0.0 )then

!     ------------------------------------------------------------------
!             proceed with spreading this data point.  first determine
!             a new radius for the gridpoint.
!     ------------------------------------------------------------------

                 if (ztest .eq. 0 .and. src(is,js) .eq. 1) then
                 
                    irad = radius(is,js) + addrad(7)
                 
                 elseif (ztest .eq. 0 .and. src(is,js) .eq. 3) then
                 
                    irad = radius(is,js) + addrad(8)

                 else

                    irad = radius(is,js) + addrad(src(is,js))
                 
                 end if
              
                 if( src(is,js) .eq. 1 )then

!     ------------------------------------------------------------------
!               find the difference between the real grid point value   
!               and the mean value of the obs points around it, out to  
!               a distance of the standard radius.  
!     ------------------------------------------------------------------

                    jstrt = max( js - irad, 1 )
                    jend  = min( js + irad, LIS_rc%lnr(n) )
                    istrt = max( is - irad, 1 )
                    iend  = min( is + irad, LIS_rc%lnc(n) )
                    numval = 0
                    sumval = 0.0

                    do jj = jstrt, jend
                       do ii = istrt, iend
                          if( src(ii,jj) .eq. 1 )then
                             sumval = sumval + mrgp(ii,jj)
                             numval = numval + 1
                          endif
                       enddo
                    enddo

                    if( numval .gt. 1 )then

!     ------------------------------------------------------------------
!                 use an exponential distribution to describe the   
!                 frequency of occurrence of real precip amounts.   
!                 evaluate amt difference between local mean parsed 
!                 real precip value and pt's value. adjust spread   
!                 radius accordingly.   
!     ------------------------------------------------------------------

                       meanv = sumval / real( numval )
                    
                       if( nint(meanv) .eq. 0 )then
                          irad = irad + varrad(1)
                       else
                          vdiff = abs( mrgp(is,js) - meanv )
                          sigmav = meanv
                          if( vdiff .le. (0.5 * sigmav) )then
                             irad = irad + varrad(2)
                          elseif( vdiff .le. (2.0 * sigmav) )then
                             irad = irad + varrad(3)
                          else
                             irad = irad + varrad(4)
                          endif
                       endif
                    endif
                 endif
              
!     ------------------------------------------------------------------
!             now set the max and min radius values based upon source   
!             type. 
!     ------------------------------------------------------------------

                 irad = min( irad, maxrad(src(is,js)) )
                 irad = max( irad, minrad(src(is,js)) )

!     ------------------------------------------------------------------
!             make a real version of the pt's spread radius and calc the
!             max and min i/j coords of the grid pts to be spread to.   
!             only spread to box of +/- (irad-1) to save cpu time.  
!     ------------------------------------------------------------------

                 rad  = real( irad )
                 mini = max( (is - irad + 1), 1 )
                 minj = max( (js - irad + 1), 1 )
                 maxi = min( (is + irad - 1), LIS_rc%lnc(n) )
                 maxj = min( (js + irad - 1), LIS_rc%lnr(n) )
                 ri   = real( is )
                 rj   = real( js )
              
!     ------------------------------------------------------------------
!             loop thru the major landmass grid pts to which this   
!             precip value will be spread.  
!     ------------------------------------------------------------------

                 do jj = minj, maxj
                    do ii = mini, maxi

                       if( LIS_LMLC(n)%landmask(ii,jj) .gt. 0 )then
                         
                     
!     ------------------------------------------------------------------
!                   calculate the weights.  each estimate type uses 
!                   its own radii and weights.  
!     ------------------------------------------------------------------

                          rsqd = sumsqr( ri, real(ii), rj, real(jj) )
                          rrad = sqrt( rsqd )
                       
                          if( rrad .le. rad )then
                             idis = int( 20.0 * (rrad / rad) )
                             srcdis = srcwt * wt(idis)
                             wgt(ii,jj) = wgt(ii,jj) + srcdis
                             tmp(ii,jj) = tmp(ii,jj) + (srcdis * mrgp(is,js))
                          endif
                       
                       endif
                    
                    enddo
                 enddo

              endif
           endif
        endif
     enddo
  enddo

!     ------------------------------------------------------------------
!     now, loop through all points in the hemisphere and divide out the 
!     sum of the weights from the total values.  this will give a final 
!     weighted amount.  if the point contains a real (src = 1) then 
!     ensure the real value is replaced over the merged analysis.   
!     ------------------------------------------------------------------
  

! ----------------------------------------------------------------------
! Don't spread obs if 'background' source is CMORPH.
! ----------------------------------------------------------------------
  do js = 1, LIS_rc%lnr(n)
     do is = 1, LIS_rc%lnc(n)
        if (cmorphpixel(is,js) .eqv. .false.) then

!     ------------------------------------------------------------------
!         if merged values were spread to this grid pt, calculate a new 
!         merged amount.
!     ------------------------------------------------------------------

           if( wgt(is,js) .gt. 0.0 )then
              if( src(is,js) .ne. 1 )then
                 mrgp(is,js) = tmp(is,js) / wgt(is,js)
              endif
              cnt = cnt + 1
              tmp(is,js) = 0.0
              wgt(is,js) = 0.0
           else
              mrgp(is,js) = -9999.0
           endif
        endif
     enddo
  enddo

  write(LIS_logunit,6000) cnt

  return

!     ------------------------------------------------------------------
!     format statement
!     ------------------------------------------------------------------

6000 format (/, 1x, 55('-'), &
          /, 3x, 'routine barnes:',&
          /, 5x, '# of pts to which values were spread = ', i6,&
          /, 1x, 55('-'))
  
end subroutine AGRMET_pcp_barnes
