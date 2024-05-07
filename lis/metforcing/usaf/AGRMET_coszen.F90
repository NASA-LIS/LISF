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
! !ROUTINE: AGRMET_coszen
!  \label{AGRMET_coszen}
!
! !REVISION HISTORY:
! 27 Apr 2016: James Geiger; Refactored from AGRMET_tr_coeffs
!
! !INTERFACE:    
subroutine AGRMET_coszen(yr1, mo1, da1, hr1, hemi, i, j, mesh, coszen)
! !USES:
   use LIS_constantsMod, only : LIS_CONST_PI, LIS_CONST_SOLAR
   use LIS_timeMgrMod,   only : LIS_get_julhr

   implicit none
! !ARGUMENTS: 
   integer, intent(in)  :: yr1, mo1, da1, hr1
   integer, intent(in)  :: hemi
   integer, intent(in)  :: i, j
   integer, intent(in)  :: mesh
   real,    intent(out) :: coszen

! !DESCRIPTION:
! Compute the cosine of the zenith angle.
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[yr1]
!   current year.
!  \item[mo1]
!   current month.
!  \item[da1]
!   current day.
!  \item[hemi]
!   hemisphere.
!  \item[i]
!   row index.
!  \item[j]
!   column index.
!  \item[mesh]
!   mesh size of the input data.
!  \item[coszen]
!   solar zenith angle.
!  \item[julday]
!   julian day of the calendar year. 
!  \item[frac]
!   intermediate step to fracsq. 
!  \item[fracsq]
!   fraction of sunlight based on time of year.  
!  \item[asolcn]
!   amount of solar radiation entering the top of
!   the atmosphere for a particular time of the year.
!  \item[alat]
!    latitude of point.  
!  \item[alon]
!    longitude of the point.  
!  \item[slat]
!   sine of sun's elevation angle.   
!  \item[pid180]
!   coefficient for sun' angle.
!  \item[clat]
!   cosine of the sun's elevation angle.
!  \item[deltim]
!   longitudnal zulu time difference.
!  \item[sdec]
!   solar declination angle
!  \item[sdel]
!   solar zenith angle
!  \item[cdel]
!   cosine of the sun's declination angle. 
!  \item[ztime]
!   current zulu time.   
!  \item[hrangl]
!   hour angle of the sun.   
!  \item[pid12]
!   coefficient for hour angle.  
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_get\_julhr](\ref{LIS_get_julhr}) \newline
!    computes the current julian hour
!  \item[pstoll](\ref{pstoll}) \newline
!    converts polar stereographic grid i/j points to lat/lon
!  \end{description}
!EOP

   real                         :: frac 
   real                         :: fracsq   
   real                         :: asolcn   
   real                         :: alat
   real                         :: alon
   real                         :: slat 
   real                         :: clat 
   real                         :: deltim   
   real                         :: sdec 
   real                         :: sdel 
   real                         :: cdel 
   real                         :: ztime
   real                         :: hrangl   
   real                         :: pid180   
   real                         :: pid12
   integer                      :: julday
   integer                      :: julhr
   integer                      :: days   ( 12 ) 

   data days   / 0,31,59,90,120,151,181,212,243,273,304,334 /
   data pid180 / 0.0174532925 /  
   data pid12  / 0.2617993878 /  

!     ------------------------------------------------------------------
!     compute a julian day using month and day fm array datime.
!     ------------------------------------------------------------------

   julday = days(mo1) + da1

!     ------------------------------------------------------------------
!     adjust the julian day for leap years.  
!     ------------------------------------------------------------------
!
   if ( (mo1 .gt. 2) .and. (mod(yr1,4) .eq. 0) ) then  
      julday = julday + 1
   endif

!     ------------------------------------------------------------------
!     adjust the solar constant based on julian day.
!     ------------------------------------------------------------------

   frac   = (1.00014 + &
      (0.016726 * cos(2.0 * LIS_CONST_PI * (julday - 2) / 365.0)))
   fracsq = frac * frac  
   asolcn = LIS_CONST_SOLAR * fracsq  

!     ------------------------------------------------------------------
!     Michael Shaw - calc the sun's incident angle based on the location's 
!     latitude and the location's sun time based on longitude. 
!     need to tell pstoll what mesh we're working with based on the 
!     native grid bounds.   e.g., 1024/64 = 16 (th mesh), 512/64 = 
!     8 (th mesh). 
!     ------------------------------------------------------------------
   call pstoll( hemi, 0, float(i), float(j), mesh, alat, alon )

   slat   = sin( alat * pid180 )   
   clat   = cos( alat * pid180 )   
   deltim = alon / 15.0  

!     ------------------------------------------------------------------
!     compute the sun's declination angle and the sun's zenith angle.   
!     ------------------------------------------------------------------

   sdec   = 23.5 * (LIS_CONST_PI / 180.0) * &
      sin(2.0 * LIS_CONST_PI * (julday - 80) / 365.0) 
   sdel   = sin(sdec)  
   cdel   = cos(sdec)  

   call LIS_get_julhr(yr1, mo1, da1, hr1, 0, 0,julhr)
   ztime  = float( mod(julhr,24) ) 
   hrangl = ( 12.0 - ztime + deltim ) * pid12
   coszen = ( slat * sdel ) + ( clat * cdel * cos( hrangl ) )

end subroutine AGRMET_coszen
