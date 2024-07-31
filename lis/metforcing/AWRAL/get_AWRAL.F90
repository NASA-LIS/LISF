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
! !ROUTINE: get_AWRAL
! \label{get_AWRAL}
!
! !REVISION HISTORY:
! 30 Jan 2017: Sujay Kumar, Initial version
!
! !INTERFACE:
subroutine get_AWRAL(n, findex)

! !USES:
  use LIS_coreMod, only     : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only  : LIS_tick, LIS_get_nstep
  use LIS_logMod,      only : LIS_logunit, LIS_endrun
  use AWRAL_forcingMod, only : AWRAL_struc
  use LIS_constantsMod, only : LIS_CONST_CDAY

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n 
  integer, intent(in) :: findex

! !DESCRIPTION:
!  Opens, reads, and interpolates AWRAL forcing. 
!  At the beginning of a simulation, the code reads the most
!  recent past data (nearest the hour interval), and the nearest
!  future data. These two datasets are used to temporally 
!  interpolate the data to the current model timestep. 
!
!  The arguments are: 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the AWRAL data times
!  \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    Computes the neighbor, weights for bilinear interpolation
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    Computes the neighbor, weights for conservative interpolation
!  \item[AWRALfile](\ref{AWRALfile}) \newline
!    Puts together appropriate file name for 1-hour intervals
!  \item[read\_AWRAL](\ref{read_AWRAL}) \newline
!    Interpolates AWRAL data to LIS grid
!  \end{description}
!EOP
   
!== Local Variables =======================
    integer :: ferror_AWRAL           ! Error flags for precip data sources
    integer :: order

    real*8  :: timenext
    real*8  :: AWRAL_file_timep       ! End boundary time for STAGEIV file
    real*8  :: AWRAL_file_timec       ! End boundary time for STAGEIV file

    integer :: doyp, yrp, mop, dap, hrp, mnp, ssp
    integer :: doyc, yrc, moc, dac, hrc, mnc, ssc
    integer :: doyn, yrn, mon, dan, hrn, mnn, ssn
    real    :: gmtp, gmtc, gmtn, tsp, tsc, tsn                    

    integer :: index1

!=== End Variable Definition =======================

    yrn = LIS_rc%yr  !get time next day
    mon = LIS_rc%mo
    dan = LIS_rc%da
    hrn = 0
    mnn = 0
    ssn = 0
    tsn = LIS_CONST_CDAY
    call LIS_tick( timenext, doyn, gmtn, yrn, mon, dan, hrn, mnn, ssn, tsn )

!-- Determine LIS's current time and the time of the AWRAL file:
    yrp = LIS_rc%yr !get time previous day
    mop = LIS_rc%mo
    dap = LIS_rc%da
    hrp = 0
    mnp = 0
    ssp = 0
    tsp = -LIS_CONST_CDAY
    call LIS_tick( AWRAL_file_timep, doyp, gmtp, yrp, mop, dap, hrp, mnp, ssp, tsp )

!-- AWRAL product time; end accumulation time data
    yrc = LIS_rc%yr !get current time
    moc = LIS_rc%mo
    dac = LIS_rc%da
    hrc = 0
    mnc = 0
    ssc = 0
    tsc = 0
    call LIS_tick( AWRAL_file_timec, doyc, gmtc, yrc, moc, dac, hrc, mnc, ssc, tsc )

!-- Ensure that data is found during first time step
    if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then 
         LIS_rc%rstflag(n) = 0
    endif

!-- Check for and get AWRAL data
    ferror_AWRAL = 0
    order = 2

  ! LIS timestep < AWRAL time interval (1da)
    if( LIS_rc%ts < AWRAL_struc(n)%ts ) then
      write(LIS_logunit,*) "[ERR] AWRAL READER CANNOT HANDLE TIME INTERPOLATION "
      write(LIS_logunit,*) "[ERR] TIMESTEP < 86400 secs -- AT THIS TIME.  STOPPING ..."
      call LIS_endrun
    endif


  ! LIS timestep == AWRAL time interval (1da)
    if( LIS_rc%ts == AWRAL_struc(n)%ts ) then

     ! Determine and return filename of AWRAL file 
       write(LIS_logunit,*) '[INFO] Getting new AWRAL data, no time interpolation necessary'
     ! Open, read, and reinterpolate AWRAL field to LIS-defined grid
       call read_AWRAL ( order, n, findex, yrc, doyc, ferror_AWRAL )
     ! Assign latest AWRAL file time to stored AWRAL time variable
       AWRAL_struc(n)%AWRALtime = AWRAL_file_timec
       index1 = LIS_domain(n)%gindex(1,1)
    else
      write(LIS_logunit,*) "[ERR] AWRAL READER CANNOT HANDLE LIS "
      write(LIS_logunit,*) "[ERR] TIMESTEP > 86400 secs -- AT THIS TIME.  STOPPING ..."
      call LIS_endrun
    endif

end subroutine get_AWRAL

