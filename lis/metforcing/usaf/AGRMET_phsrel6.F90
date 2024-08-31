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
! !ROUTINE: AGRMET_phsrel6
!  \label{AGRMET_phsrel6}
!
!
! !REVISION HISTORY:
!
!     21 feb 01  initial version based on original routine phsrel.......
!                ..........................................mr gayno/dnxm
!     10 jun 02  modified to reflect change from rtneph to cdfs2 data...
!                ..........................................mr gayno/dnxm
!
!     3 nov 05 Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine AGRMET_phsrel6 ( n, estpcp, j6hr, p6,&
     quad9r, cdfs2est, source, &
     obswch, p3, pathpcp)
! !USES: 
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use LIS_historyMod, only : LIS_writevar_bin
  use LIS_timeMgrMod, only : LIS_julhr_date
  use LIS_mpiMod

  implicit none
! !ARGUMENTS:     
  integer, intent(in)               :: n
  integer                           :: j3hr
  integer, intent(in)               :: j6hr
  integer                           :: k
  integer, intent(in)               :: obswch
  integer, intent(inout)            :: source(LIS_rc%lnc(n),LIS_rc%lnr(n),4)  
  real,    intent(inout)            :: cdfs2est(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  real,    intent(inout)            :: estpcp(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  real,    intent(out)              :: p3(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  real,    intent(in)               :: p6(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,    intent(in)               :: quad9r
  character(len=*), intent(in)      :: pathpcp
  
! !DESCRIPTION:
!
!     to parse 6-hourly real observed precipitation amounts into 3-hour
!     amounts, and to write all output arrays to file.  
!
!     \textbf{Method} \newline
!     
!     - fill the 3 hrly observed output array with 9999.0.  \newline
!     - call parse routine to divide up the 6-hrly real precip  
!       (rain gauge) amounts into 3-hour precip amounts. \newline
!     - write out the 3-hour parsed real precip amounts to file. \newline
!     - write out the 3-hour estimated precip amounts to file.  \newline
!     - write out the 3-hour estimate sources to file.  \newline
!     - write out the 3-hour cdfs2-based precip estimates to file. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[cdfs2est]   cdfs2-based precip estimate (mm/3hr)
!   \item[estpcp]     3-hrly estimated precip amounts (mm)
!   \item[hemi]       hemisphere
!   \item[j6hr]       current 6-hrly start time (julian hour)
!   \item[j3hr]       current 3-hrly time (julian hour)
!   \item[k]          loop index, 'i'th 3-hrly time period
!   \item[obswch]     observation processing control switch (1 - use obs)
!   \item[p3]         3-hrly parsed real precip amounts (mm)
!   \item[p6]         6-hrly real (rain gauge) precip amounts (mm)
!   \item[quad9r]     scalar value of 9999.0
!   \item[source]     array that contains source-of-the-estimate
!                     information with specified source flags
!  \item[pathpcp]     Path to agrmet precip files
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_parse6] (\ref{AGRMET_parse6}) \newline
!  \end{description}
!EOP
  character*100                     :: ofil
  character*10                      :: date10_03
  integer                           :: yr1,mo1,da1,hr1
  character*4                       :: fyr
  character*2                       :: fmo,fda
  integer                           :: ftn

!
!     ------------------------------------------------------------------
!     executable code begins here ... initialized parsed precip array
!     with "missing" values
!     ------------------------------------------------------------------

  p3 = quad9r

!     ------------------------------------------------------------------
!     parse the 6-hrly real (rain gauge) precip obs (p6) into
!     3-hrly amounts (p3) for this 6 hour period
!     ------------------------------------------------------------------

  if (obswch == 1) then

     write(LIS_logunit,*)" "
     write(LIS_logunit,*)"- PARSE 6 HOURLY RAIN GAUGE DATA INTO 3 HOURLY AMTS"
     
     call AGRMET_parse6( n, p6, estpcp, p3, LIS_rc%lnc(n), LIS_rc%lnr(n))

  end if

end subroutine AGRMET_phsrel6

