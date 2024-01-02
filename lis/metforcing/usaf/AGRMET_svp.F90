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
! !ROUTINE: AGRMET_svp
! \label{AGRMET_svp}
!
! !REVISION HISTORY:
!
!     30 jul 99  initial version..........lt col (ima) ken mitchell/dnxm
!     07 sep 99  ported to ibm sp-2.  updated prolog.  added intent
!                attributes to arguments...................mr gayno/dnxm
!     01 aug 05 Sujay Kumar, Adopted in LIS
!
! !INTERFACE:    
subroutine AGRMET_svp( qs, es, p, t )   

  implicit none
! !ARGUMENTS: 
  real,  parameter              :: cpv = 1870.0
  real,  parameter              :: cw  = 4187.0
  real,  parameter              :: eso = 611.2
  real,  parameter              :: rv  = 461.5
  real,  parameter              :: to  = 273.15      
  
  real,     intent(out)        :: es   
  real,     intent(in)         :: p
  real,     intent(out)        :: qs   
  real,     intent(in)         :: t
!
! !DESCRIPTION:
!     calculates saturation vapor pressure and saturation 
!     specific humidity
!
!     \textbf{Method} \newline
!
!     - call function e to calculate saturation vapor pressure \newline
!     - calculate saturation specific humidity based on saturation
!       vapor pressure \newline
! 
!   The arguments are:
!   \begin{description}
!     \item[e]
!        saturation vapor pressure function
!     \item[es]
!       output saturation vapor pressure (pascals)
!     \item[p]
!        input air pressure (pascals) 
!     \item[qs]
!       output saturation mixing ratio (kg kg-1)  
!     \item[t]
!        input air temperature (k)
!   \end{description}
!EOP
  real                         :: lw
!     ------------------------------------------------------------------
!     execution starts here....calculate the saturation vapor pressure.
!     ------------------------------------------------------------------
  lw = 2.501e6 - ( cw - cpv ) * ( t - to )
  es = eso * exp(lw *(1/to - 1/t)/rv)    

!     ------------------------------------------------------------------
!     calculate a saturation mixing ratio. 
!     ------------------------------------------------------------------
  
  qs = 0.622 * es / p 
  
  return

end subroutine AGRMET_svp
    
