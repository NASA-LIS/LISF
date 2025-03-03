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
! !ROUTINE: AGRMET_parse6
!  \label{AGRMET_parse6}
!
! !REVISION HISTORY:
!
!     21 feb 01  initial version based on original parse routine........
!                ..........................................mr gayno/dnxm
!      4 nov 05  incorporated into LIS, sujay kumar
!      5 nov 10  Modified to include observations that, due to land mask
!                resolution limits, were not included as they fell over 
!                grid cells flagged as "water" in coarser resolutions
!                 ......................................Michael Shaw/WXE
!      5 nov 10  Modified to prevent array index errors.Michael Shaw/WXE
! 
! !INTERFACE:    
subroutine AGRMET_parse6( n, p6, estpcp, p, imax, jmax )
! !USES: 
  use LIS_LMLCMod, only : LIS_LMLC

  implicit none
  
  integer                           :: i
  integer, intent(in)               :: n
  integer, intent(in)               :: imax
  integer                           :: j
  integer, intent(in)               :: jmax
  integer                           :: k
  real,    intent(in)               :: estpcp(imax,jmax,4)
  real                              :: e6
  real,    intent(out)              :: p(imax,jmax,4)
  real,    intent(in)               :: p6(imax,jmax)

! !DESCRIPTION:
!
!    to break-up 'observed' 6-hrly real rain gauge amounts
!    into 3-hrly amounts. \newline
!
!    loop thru the points for a hemisphere. for each
!    land point with a valid 6-hrly observed
!    real precip amount... \newline
!    a. build a 6-hrly estimated precip amount from the
!       3-hrly estimated amounts in the two 3-hrly periods. \newline
!    b. if the 6-hrly estimated precip amount = 0.0
!       break the 6-hrly real amount into two equal
!       parts. \newline
!    c. if the 6-hrly estimated precip amount > 0.0
!       break the 6-hrly real amount into two amounts,
!       each determined by the ratio of the 3-hrly
!       estimated amount to the 6-hrly estimated amount. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[estpcp]    3-hrly estimated precip amounts (mm)
!   \item[e6]        sum of 3-hrly precip estimates for the 6-hrly period
!   \item[hemi]      hemisphere
!   \item[i]         loop index
!   \item[imax]      number of gridpoints in east/west direction
!   \item[j]         loop index
!   \item[jmax]      number of gridpoints in north/south direction
!   \item[k]         loop index   
!   \item[land]      point processing switches
!   \item[p]         3-hrly parsed real (rain gauge) precip amounts (mm)
!   \item[p6]        6-hrly real (rain gauge) precip amounts (mm)
!  \end{description}
!EOP
  
!     ------------------------------------------------------------------
!     executable code starts here...process land points with a
!     6-hrly real (rain gauge) precip amount...   
!     ------------------------------------------------------------------
  logical     :: bounds_check 

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
!        LAND_WATER : if( LIS_LMLC(n)%landmask(i,j)     .gt. 0 .or. &
!                         (i .gt. 1 .and. i .lt. imax .and. j .gt. 1 .and. j .lt. jmax) .and. &
!                        (LIS_LMLC(n)%landmask(i+1,j)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i,j+1)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i+1,j+1) .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i-1,j)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i,j-1)   .gt. 0 .or. &
!                         LIS_LMLC(n)%landmask(i-1,j-1) .gt. 0) )then 
        LAND_WATER : if(bounds_check) then 

           if( p6(i,j) .gt. -9990.0 )then  

!     ------------------------------------------------------------------
!             sum the two estimates for the six-hourly period.
!             if all of the estimates are invalid (9999.0) then set 
!             the sum to be invalid (9999.0), else, sum only the valid
!             estimates into the six hour sum.
!     ------------------------------------------------------------------

              e6 = 0.0 
              if ( (estpcp(i,j,1) .lt. -9990.0) .and. &
                   (estpcp(i,j,2) .lt. -9990.0) ) then
                 e6 = -9999.0
              else
                 do k = 1, 2   
                    if ( estpcp(i,j,k) .gt. -9990.0) &
                         e6 = e6 + estpcp(i,j,k)
                 enddo
              endif
              
!     ------------------------------------------------------------------
!             if there were no valid 3-hrly estimates, parse the 6-
!             hrly real precip amount into two equal parts. otherwise,
!             parse real amounts based on the ratio of the 3-hrly
!             estimate to the 6-hrly total estimate.
!     ------------------------------------------------------------------
      
              if( (nint(e6) .eq. 0) .or. (e6 .lt. -9998.0) )then
                p(i,j,1) = 0.5 * p6(i,j)
                p(i,j,2) = 0.5 * p6(i,j)
                
              else
                 
                 if ( estpcp(i,j,1) .gt. -9990.0 )  &
                      p(i,j,1) = ( estpcp(i,j,1)/e6 ) * p6(i,j)
                 if ( estpcp(i,j,2) .gt. -9990.0 )  &
                      p(i,j,2) = ( estpcp(i,j,2)/e6 ) * p6(i,j)

              endif
              
           endif
           
        endif LAND_WATER
        
     enddo
  enddo
  
  return
  
end subroutine AGRMET_parse6
