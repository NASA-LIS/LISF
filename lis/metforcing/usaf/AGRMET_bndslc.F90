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
! !ROUTINE: AGRMET_bndslc
!  \label{AGRMET_bndslc}
!
! !REVISION HISTORY:
!
!     09 jul 94  initial version..............capt bertone, sysm(agromet)   
!     16 sep 99  ported to ibm sp2. removed box logic. added intent
!                attributes to arguments.  added call to utility 
!                routine pstoll, which replaces an obsolete utility
!                routine grdtll............................mr gayno/dnxm
!     01 aug 05  Sujay Kumar, Adopted in LIS with modifications
!
! !INTERFACE:    
subroutine AGRMET_bndslc(i,j, hemi, latbnd, lonslc ) 

  implicit none
! !ARGUMENTS: 
  integer,   intent(in)         :: hemi  
  integer,   intent(out)        :: latbnd
  integer,   intent(out)        :: lonslc
  integer, intent(in)           :: i
  integer, intent(in)           :: j  
!
! !DESCRIPTION:
!     to determine and return to the parent routine the   
!     latitude band and longitude slice into which this point falls.
!  
!     \textbf{Method} \newline
!    
!     - find the band and slice into which this point's 
!       lat and lon fall. \newline
!
!  The arguments and variables are: 
!  \begin{description}  
!   \item[hemi]      hemisphere (1 = north, 2 = south )  
!   \item[i]         i-coordinate of agrmet grid point  
!   \item[j]         j-coordinate of agrmet grid point    
!   \item[latbnd]    latitude band   
!   \item[lon]       longitude of agrmet point
!   \item[lonslc]    longitude slice 
!   \item[ri]        real i-coordinate of agrmet grid point  
!   \item[rj]        real j-coordinate of agrmet grid point  
!  \end{description}
!   
!  The routines invoked are: 
!  \begin{description}
!  \item[pstoll](\ref{pstoll}) \newline
!   computes the lat lon values of a PS grid point
!  \end{description}
!EOP

  real                             :: ri
  real                             :: rj
  real                             :: lat
  real                             :: lon  
!     ------------------------------------------------------------------
!     executable code begins here... convert the point's grid   
!     coord's to a lat and lon by invoking utility routine.
!     ------------------------------------------------------------------
    
  ri = float ( i )   
  rj = float ( j )  
 
  call pstoll( hemi, 1, ri, rj, 8, lat, lon )

  if ( hemi .eq. 2 ) then   
     if ( lat .le. -60. ) then   
        latbnd = 1
     elseif ( lat .le. -42.5 ) then  
        latbnd = 2
     elseif ( lat .le. -25.) then
        latbnd = 3
     elseif ( lat .le. -10.) then
        latbnd = 4
     elseif ( lat .lt. 0. ) then 
        latbnd = 5
     endif
  else  
     if ( lat .lt. 10. ) then
        latbnd = 5
     elseif ( lat .lt. 25. ) then
        latbnd = 4
     elseif ( lat .lt. 42.5 ) then   
        latbnd = 3
     elseif ( lat .lt. 60. ) then
        latbnd = 2
     else
        latbnd = 1
     endif
  endif
  
!     ------------------------------------------------------------------ 
!     determine the point's longitude slice from its longitude.
!     ------------------------------------------------------------------
    
  if ( lon .lt. 45. ) then  
     lonslc = 1  
  elseif ( lon .lt. 90. ) then  
     lonslc = 2  
  elseif ( lon .lt. 135. ) then 
     lonslc = 3  
  elseif ( lon .lt. 180. ) then 
     lonslc = 4  
  elseif ( lon .lt. 225. ) then 
     lonslc = 5  
  elseif ( lon .lt. 270. ) then 
     lonslc = 6  
  elseif ( lon .lt. 315. ) then 
     lonslc = 7  
  else  
     lonslc = 8  
  endif
  
  return

end subroutine AGRMET_bndslc
