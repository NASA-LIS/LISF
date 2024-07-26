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
! !ROUTINE: AGRMET_phsrel12
!  \label{AGRMET_phsrel12}
!
! !REVISION HISTORY:
!
!     21 feb 01  initial version based on original routine phsrel.......
!                ..........................................mr gayno/dnxm
!     10 jun 02  modified to reflect change from rtneph to cdfs2 data...
!                ..........................................mr gayno/dnxm
!
!     3 nov  05 Sujay Kumar, incorporated into LIS
!
! !INTERFACE:    
subroutine AGRMET_phsrel12( n, p12, j6hr, &
     estpcp, source, cdfs2est,&
     quad9r, obswch, p3x, pathpcp)
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_writevar_bin
  use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use LIS_fileIOMod,  only : LIS_readData
  use LIS_timeMgrMod, only : LIS_julhr_date
  use LIS_mpiMod

  implicit none

 ! !ARGUMENTS:    
  integer,       intent(in)    :: n
  integer,       intent(in)    :: j6hr
  integer,       intent(in)    :: obswch    
  real,          intent(inout) :: cdfs2est(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  real,          intent(inout) :: estpcp(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  real,          intent(out)   :: p3x(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  real,          intent(in)    :: p12(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real,          intent(in)    :: quad9r
  integer,       intent(inout) :: source(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  character(len=*), intent(in) :: pathpcp

! !DESCRIPTION:
!
!     to parse 12-hourly real observed precipitation amounts into 3-hour
!     amounts, and to write all precip arrays to file.  
!
!     \textbf{Method} \newline
!   
!     - fill the 3 hrly observed output array with 9999.0.  \newline
!     - read in estimates from 6 and 9 hours ago. \newline
!     - call parse routine to divide up the 12-hrly real precip  
!       (rain gauge) amounts into 3-hour precip amounts. \newline
!     - write out the 3-hour parsed real precip data (over the 12
!       hour period) to file. \newline
!     - write out the 3-hour estimated precip amounts to file. \newline
!     - write out the 3-hour estimate sources to file.  \newline
!     - write out the 3-hour cdfs2-based precip estimates to file. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[cdfs2est]   cdfs2-based precip estimate (mm/3hr)
!   \item[e]          array of precip estimates (4 - 3hrly periods)
!   \item[estpcp]     array of precip estimates (last 2 3hrly periods)
!   \item[hemi]       hemisphere
!   \item[imax]       number of gridpoints in east/west direction
!   \item[jmax]       number of gridpoints in north/south direction
!   \item[j6hr]       current 6-hrly start time (julian hour)
!   \item[obswch]     observation processing control switch (1 - use obs)
!   \item[p3]         3-hrly parsed real precip amounts (mm)
!   \item[p12]        12-hrly real (rain gauge) precip amounts (mm)
!   \item[quad9r]     scalar value of -9999.0
!   \item[source]     array that contains source-of-the-estimate
!                     information with specified source flags
!  \item[pathpcp]     Path to agrmet precip files
!  \item[ifil]        input file name including directory path
!  \item[ofil]        output file name including directory path
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_parse12] (\ref{AGRMET_parse12}) \newline
!  \end{description}
!EOP
  integer                      :: k 
  real                         :: e(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  real                         :: p3(LIS_rc%lnc(n),LIS_rc%lnr(n),4)
  character*100                :: ofil
  character*100                :: ifil
  character*10                 :: date10_03
  logical                      :: exists
  integer                      :: j3hr
  integer                      :: yr1,mo1,da1,hr1
  character*4                  :: fyr
  character*2                  :: fmo,fda
  real                         :: gridDesc(6)
  integer                      :: ftn

!     ------------------------------------------------------------------
!     in order to parse the 12 hourly rain gauge amounts into 
!     3 hourly amounts, need 12 hours worth of estimates.  the most
!     recent 6 hours of estimates are stored in estpcp.  the remaining
!     estimates were calculated the last time the precip program ran,
!     so read them in.  store the 4 time periods of estimates in
!     array "e" for passing to routine parse12.
!     ------------------------------------------------------------------

  p3 = quad9r
  e  = quad9r 
  USING_OBS : if (obswch == 1) then
     write(LIS_logunit,*)" "
     write(LIS_logunit,*)"- PARSE 12 HOURLY RAIN GAUGE DATA INTO 3 HOURLY AMTS"
     write(LIS_logunit,*)"- READ ESTIMATES FROM 6 AND 9 HOURS AGO"
     write(LIS_logunit,*)" "


     e(:,:,1) = estpcp(:,:,1)
     e(:,:,2) = estpcp(:,:,2)
!     ------------------------------------------------------------------
!       fill the most recent three hourly estimates with data
!       calculated during this precip run.
!     ------------------------------------------------------------------
     
     e(:,:,3) = estpcp(:,:,3)
     e(:,:,4) = estpcp(:,:,4)

!     ------------------------------------------------------------------
!       parse 12 hourly rain gauge amounts into 3 hourly amounts, then
!       write the 3 hourly amounts to file.
!     ------------------------------------------------------------------
     
     call AGRMET_parse12( n, p12, e, p3, LIS_rc%lnc(n), LIS_rc%lnr(n))
  end if USING_OBS
  
  p3x = p3 
  
!     ------------------------------------------------------------------
!     write out estimates for last half of 12 hour period (the
!     estimates for the first half of the period were calculated
!     the last time the precip program ran.
!     ------------------------------------------------------------------


end subroutine AGRMET_phsrel12
   
