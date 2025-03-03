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
! !ROUTINE: AGRMET_make03
!  \label{AGRMET_make03}
!
! !REVISION HISTORY:
!
!     15 jun 96  initial version...........................mr moore/sysm
!     10 apr 97  brought up to software standards.  modified to more
!                efficiently convert supergrid merged precip array into
!                box arrays, write to file, and sum into 12-hr merged
!                amount.................................ssgt miller/sysm
!     02 may 97  corrected error in output caused by dividing 12-hour
!                merged precipitation values by 10800. updated prolog
!                and brought up to standards............capt andrus/sysm
!      7 oct 99  ported to ibm sp-2, updated prolog, incorporated 
!                FORTRAN 90 features.  modified the part of code that 
!                divides the merged values by 10800 so that only the 
!                values needed by FLUX3 (the 3hrly merged values) are
!                divided by 10800 (to get the units in terms of a flux)
!                ....................................capt hidalgo/agrmet
!     05 jan 01  modified variable names to reflect 6-hourly cycles.....
!                ..........................................mr gayno/dnxm
!     3 nov 05 Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine AGRMET_make03( n, mrgp, mrgp6)
! !USES: 
  use LIS_coreMod,  only  : LIS_rc

  implicit none
  
! !ARGUMENTS: 
  integer,       intent(in)    :: n
  real,          intent(inout) :: mrgp(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real,          intent(inout) :: mrgp6(LIS_rc%lnc(n), LIS_rc%lnr(n))
!
! !DESCRIPTION:
!
!      to write the 3-hrly merged precip amounts to file and
!      accumulate the 3-hrly merged precip amounts into a
!      6-hrly merged precip total.
!
!      
!      \textbf{Method} \newline
!      1. loop thru grid points in the hemisphere
!        - sum the 3-hrly amts into 6-hrly amts \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[i]           loop index
!   \item[imax]        number of gridpoints in east/west direction 
!   \item[j]           loop counter
!   \item[jmax]        number of gridpoints in north/south direction
!   \item[mrgp]        supergrid 3-hrly merged precip amounts (mm)
!   \item[mrgp6]       supergrid 6-hrly merged precip amounts (mm)
!  \end{description}
!EOP
  integer                      :: i
  integer                      :: j

!     ------------------------------------------------------------------
!     executable code begins here...loop thru points in hemisphere
!     ------------------------------------------------------------------

  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)

!         --------------------------------------------------------------
!         sum 3-hrly amounts into 6-hrly amounts.
!         --------------------------------------------------------------

        if( mrgp(i,j) .gt.-9990.0)then
           if( mrgp6(i,j) .gt.-9990.0 )then
              mrgp6(i,j) = mrgp6(i,j) + mrgp(i,j)
           else
              mrgp6(i,j) = mrgp(i,j)
           endif
           
!         --------------------------------------------------------------
!         convert the 3-hr merged precip values from units of mm/3hr
!         to units of mm/s.  these are the required units for flux3.
!         --------------------------------------------------------------

           mrgp(i,j) = mrgp(i,j) / 10800.0
        else
           mrgp(i,j) = LIS_rc%udef
        endif
        
     enddo
  enddo

  return
  
end subroutine AGRMET_make03
