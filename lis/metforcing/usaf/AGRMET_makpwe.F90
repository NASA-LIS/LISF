!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
!BOP
! 
! !ROUTINE: AGRMET_makpwe
! \label{AGRMET_makpwe}
!
! !REVISION HISTORY: 
!    04 dec 97  initial version.....mr moore, capt andrus/dnxm(agromet) 
!     7 oct 99  ported to ibm sp-2, updated prolog, incorporated
!                FORTRAN 90 features................capt hidalgo/agrmet
!    28 oct 10  changed to accept 2-digit wmo block number in place
!               of afwa 6-digit bsn..........Chris Franks/16WS/WXE/SEMS
! !INTERFACE: 
subroutine AGRMET_makpwe( preswx, pastwx, lat, month, wmoblk, pwprc )

! !DESCRIPTION:
!
!    to create an estimate of precipitation from the present
!    and past weather codes
!
!    \textbf{Method}
!    
!    1. calculate a 1st guess precip rate by combining   
!       the array values for the reported present and past wx. \newline
!    2. determine a wmo region adjustment. \newline
!    3. determine a gross latitude adjustment. \newline
!    4. calculate a seasonal latitude adjustment.  \newline
!    5. adjust the 1st guess precip rate with the three  
!       adjustment factors and convert the result to mm/hr. \newline
! 
!  The arguments and variables are: 
!  \begin{description}
!   \item[blkadj]    array of wmo block adjustments 
!   \item[diff]      obs site's displacement from latitude of max precip
!   \item[fguess]    first guess total present weather estimate value
!   \item[gross]     gross latitudinal adjustment
!   \item[lat]       latitude of the observation site   
!   \item[maxlat]    monthly varying latitude of max precipitation  
!   \item[month]     current month of the observation   
!   \item[pastfg]    selected first guess past wx estimate value   
!   \item[pastamt]    array of possible first guess past wx values 
!   \item[pastwx]    past weather code reported in the observation  
!   \item[presfg]    selected first guess present wx value
!   \item[presamt]    array of possible first guess present wx values  
!   \item[preswx]    present weather code reported in the observation   
!   \item[pwprc]     final precipitation rate (mm/hr)
!   \item[season]    calculated seasonal latitude adjustment
!   \item[wmoadj]    selected wmo block adjustment  
!   \item[wmoblk]    wmo region number for the report   
!  \end{description}
!EOP
  implicit none

  integer                               :: month 
  integer                               :: pastwx  
  integer                               :: preswx  
  integer                               :: pwprc  
  integer                               :: wmoblk  
  
  real                          :: blkadj(100) 
  real                          :: diff
  real                          :: fguess  
  real                          :: gross   
  real                          :: lat  
  real                          :: maxlat(12)  
  real                          :: pastamt(10)  
  real                          :: pastfg  
  real                          :: presamt(100) 
  real                          :: presfg  
  real                          :: season  
  real                          :: wmoadj  
  
  data blkadj   &
       / 1.00, 1.0, 0.60, 0.80, 1.00, 0.88, 0.90, 1.0, 0.80, 0.90, &
       1.00, 1.0, 1.00, 1.00, 0.75, 0.90, 0.60, 1.0, 1.00, 1.00, &
       0.75, 1.0, 1.00, 1.00, 0.70, 0.80, 1.00, 1.0, 1.00, 1.00, &
       1.00, 0.9, 0.95, 1.00, 0.90, 1.00, 0.95, 0.9, 1.00, 1.00, &
       1.00, 1.0, 1.00, 0.80, 1.00, 1.00, 1.00, 1.0, 1.00, 0.93, &
       0.60, 0.6, 0.80, 0.93, 0.50, 0.70, 0.80, 0.9, 1.00, 0.50, &
       1.00, 1.0, 0.85, 1.00, 1.00, 1.00, 1.00, 1.0, 1.00, 1.00, &
       1.00, 1.0, 1.00, 1.00, 1.00, 1.00, 1.00, 1.0, 1.00, 1.30, &
       1.00, 1.0, 1.00, 1.50, 1.00, 1.00, 1.00, 1.0, 1.00, 1.00, &
       1.00, 1.0, 1.00, 1.00, 1.00, 1.00, 0.80, 1.0, 1.00, 1.00 /

  data maxlat  & 
       / -10., -5., 0., 5., 10., 15., 15., 10., 5., 0., -5., -10. /

  data pastamt &  
       / 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 2.0, 1.8, 2.0, 6.0 /

  data presamt   &
       / 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, &
       0.0, 0.0, 0.0, 0.7, 0.1,  0.5, 0.7, 0.9, 0.3, 1.5, &
       0.4, 0.9, 0.1, 0.6, 0.9,  0.5, 0.1, 0.6, 0.0, 1.5, &
       0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, &
       0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, &
       0.3, 0.4, 0.5, 0.6, 2.0,  4.0, 0.1, 0.5, 1.0, 5.0, &
       0.8, 2.4, 2.6, 6.6, 5.6, 13.0, 0.1, 4.0, 2.0, 7.2, &
       0.2, 1.0, 1.6, 3.0, 2.0,  3.6, 0.1, 0.7, 0.1, 3.0, &
       1.3, 3.9, 7.0, 0.5, 3.1,  0.1, 2.2, 0.2, 0.7, 1.3, &
       2.1, 2.2, 9.0, 0.2, 2.4,  3.0, 3.3, 8.5, 1.8, 6.0 /   

!     ----------------------------------------------------------------  
!     executable code begins here... initialize 1st guess amts  
!     for present and past weather, and the wmo blk adjustment.
!     ----------------------------------------------------------------  
    
  if( (preswx .ge. 0) .or. (pastwx .ge. 0) )then
     presfg = 0  
     pastfg = 0  
     wmoadj = 1.0

!     ----------------------------------------------------------------  
!       adjust wx codes (0-99 & 0-9) to work with arrays 
!       (1-100 & 1-10).
!     ----------------------------------------------------------------  

     preswx = preswx + 1
     pastwx = pastwx + 1
     
!     ----------------------------------------------------------------  
!       calculate a first guess values.
!     ----------------------------------------------------------------  
     
     if( (preswx .ge. 1) .and. (preswx .le. 100) )then
        presfg = presamt(preswx)
     endif
     
     if( (pastwx .ge. 1) .and. (pastwx .le. 10) )then
        pastfg = pastamt(pastwx)
     endif
     
     fguess = presfg + pastfg
    
!     ----------------------------------------------------------------  
!       determine the proper wmo block adjustment.
!     ----------------------------------------------------------------  

     if( (wmoblk .ge. 1) .and. (wmoblk .le. 99) )then
        wmoadj = blkadj(wmoblk)
     else
        wmoadj = 1.0
     endif

!     ----------------------------------------------------------------  
!       determine the proper gross latitude adjustment.
!     ----------------------------------------------------------------  

     if( (lat .gt. 70.0) .or. (lat .lt. -70.0) )then
        gross = 0.14
     elseif( lat .gt. 10.0 )then
        gross = 0.48
     elseif( lat .ge. -10.0 )then
        gross = 0.93
     else
        gross = 0.48
     endif

!     ----------------------------------------------------------------  
!       determine the proper seasonal latitude adjustment.  adjustment  
!       is set to 0.0 for latitudes further than 50 degrees away from   
!       the latitude of maximum precipitation which shifts monthly.
!     ----------------------------------------------------------------  

     diff = abs( real( maxlat(month) ) - lat )
     season = max( 0.05 * (50.0 - diff), 0.0 )
     
!     ----------------------------------------------------------------  
!       calculate the final precipiation estimate amount and convert   
!       its value from mm/3hr to mm/hr.   
!     ----------------------------------------------------------------  

     pwprc = int(gross * wmoadj * (fguess + (fguess * season)) &
          + 0.5)
     
  endif
  
  return
end subroutine AGRMET_makpwe

