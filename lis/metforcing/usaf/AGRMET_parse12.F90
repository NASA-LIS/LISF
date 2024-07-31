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
! !ROUTINE: AGRMET_parse12
!  \label{AGRMET_parse12}
!
! !REVISION HISTORY:
!
!
!     15 jun 96  initial version................mr moore sysm(agromet)  
!      8 oct 99  ported to ibm sp-2, updated prolog, incorporated
!                FORTRAN 90 features.................capt hidalgo/agrmet
!     21 feb 01  renamed routine parse12 (formally parse).  removed
!                calculation of "ratio" as this field was not being
!                used......................................mr gayno/dnxm
!     3 nov  05 Sujay Kumar, Initial Code
!    30 aug  10  Modified to include observations that, due to land mask
!                resolution limits, were not included as they fell over 
!                grid cells flagged as "water" in coarser resolutions
!                ......................................Michael Shaw/WXE
!      5 nov 10  Modified to prevent array index errors.Michael Shaw/WXE
!
! !INTERFACE:    
subroutine AGRMET_parse12( n, p12, e, p, imax, jmax )
! !USES: 
  use LIS_LMLCMod, only : LIS_LMLC
  implicit none
! !ARGUMENTS: 
  integer, intent(in)               :: n
  integer, intent(in)               :: imax
  integer, intent(in)               :: jmax
  real,    intent(in)               :: e(imax,jmax,4)
  real,    intent(out)              :: p(imax,jmax,4)
  real,    intent(in)               :: p12(imax,jmax)
!
! !DESCRIPTION:
!
!   to break-up 'observed' 12-hrly real amounts into
!   3-hrly real amounts. \newline
!
!   loop thru the grid points. for each
!   major landmass point with a valid 12-hrly observed
!   real precip amount... \newline
!   a. build a 12-hrly estimated precip amount from the
!      3-hrly estimated amounts in the four 3-hrly periods. \newline
!    b. if the 12-hrly estimated precip amount = 0.0
!       break the 12-hrly real amount into four equal
!       parts. \newline
!    c. if the 12-hrly estimated precip amount > 0.0
!       break the 12-hrly real amount into four amounts,
!       each determined by the ratio of the 3-hrly
!       estimated amount to the 12-hrly estimated amount. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[e]        3-hrly estimated precip amounts (mm)
!   \item[e12]      sum of 3-hrly precip estimates for the 12 hrs 
!   \item[hemi]     hemisphere
!   \item[i]        loop index
!   \item[imax]     number of gridpoints in east/west direction
!   \item[j]        loop index
!   \item[jmax]     number of gridpoints in north/south direction
!   \item[k]        loop index   
!   \item[land]     point processing switches
!   \item[p]        3-hrly parsed real precip amounts (mm)
!   \item[p12]      12-hrly real precip amounts (mm)
!  \end{description}
!EOP  

  real                              :: e12
  integer                           :: k
  integer                           :: i
  integer                           :: j
  logical                           :: bounds_check

!     ------------------------------------------------------------------
!     executable code starts here...for each major landmass point,  
!     with a valid 12-hrly real precip amount,...   
!     ------------------------------------------------------------------

      do j = 1, jmax  
        do i = 1, imax

           bounds_check = .false. 
           if(LIS_LMLC(n)%landmask(i,j)     .gt. 0) then 
              bounds_check = .true. 
           endif
           if(i.lt.imax) then 
              if(LIS_LMLC(n)%landmask(i+1,j)   .gt. 0) then 
                 bounds_check = .true.
              endif
           endif
           if(j.lt.jmax) then 
              if(LIS_LMLC(n)%landmask(i,j+1)   .gt. 0) then 
                 bounds_check = .true.
              endif
           endif
           if(i.lt.imax.and.j.lt.jmax) then 
              if(LIS_LMLC(n)%landmask(i+1,j+1)   .gt. 0) then 
                 bounds_check = .true.
              endif
           endif
           if(i.gt.1) then 
              if(LIS_LMLC(n)%landmask(i-1,j)   .gt. 0) then 
                 bounds_check = .true.
              endif
           endif
           if(j.gt.1) then 
              if(LIS_LMLC(n)%landmask(i,j-1)   .gt. 0) then 
                 bounds_check = .true.
              endif
           endif
           if(i.gt.1.and.j.gt.1) then 
              if(LIS_LMLC(n)%landmask(i-1,j-1)   .gt. 0) then 
                 bounds_check = .true.
              endif
           endif
!           LAND_WATER : if( LIS_LMLC(n)%landmask(i,j)     .gt. 0 .or. &
!                         (i .gt. 1 .and. i .lt. imax .and. j .gt. 1 .and. j .lt. jmax) .and. &
!                         (LIS_LMLC(n)%landmask(i+1,j)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i,j+1)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i+1,j+1) .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i-1,j)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i,j-1)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i-1,j-1) .gt. 0) )then
 
           LAND_WATER : if( bounds_check) then 
            if( p12(i,j) .gt. -9990.0 )then  

!     ------------------------------------------------------------------
!             sum the four estimates for the 1st thru 4th 3-hrly periods
!             if all of the estimates are invalid (9999.0) then set 
!             the sum to be invalid (9999.0), else, sum only the valid
!             estimates into the 12 hour sum
!     ------------------------------------------------------------------

              e12 = 0.0 
              if ( (e(i,j,1) .lt. -9990.0) .and. (e(i,j,2) .lt. -9990.0) &
                   .and. (e(i,j,3) .lt. -9990.0) &
                   .and. (e(i,j,4) .lt. -9990.0)) then
                e12 = -9999.0
             else
                do k = 1, 4   
                   if ( e(i,j,k) .gt. -9990.0) then
                      e12 = e12 + e(i,j,k)
                   endif
                enddo
             endif
             
!     ------------------------------------------------------------------
!             if there were no valid 3-hrly estimates, parse the 12-
!             hrly precip amount into four equal parts. otherwise, use  
!             the 3-hrly estimates to parse the 12-hrly precip amount.  
!     ------------------------------------------------------------------
      
             if( (nint(e12).eq.0) .or. (e12.lt.-9998.0) )then
                p(i,j,1) = 0.25 * p12(i,j)
                p(i,j,2) = 0.25 * p12(i,j)
                p(i,j,3) = 0.25 * p12(i,j)
                p(i,j,4) = 0.25 * p12(i,j)
              else
                 if ( e(i,j,1) .gt. -9990.0 ) &
                      p(i,j,1) = ( e(i,j,1)/e12 ) * p12(i,j)
                if ( e(i,j,2) .gt. -9990.0 ) &
                     p(i,j,2) = ( e(i,j,2)/e12 ) * p12(i,j)
                if ( e(i,j,3) .gt. -9990.0 ) &
                     p(i,j,3) = ( e(i,j,3)/e12 ) * p12(i,j)
                if ( e(i,j,4) .gt. -9990.0 ) &
                     p(i,j,4) = ( e(i,j,4)/e12 ) * p12(i,j)
             endif
          endif
       endif LAND_WATER
    enddo
 enddo

 return
      
end subroutine AGRMET_parse12
