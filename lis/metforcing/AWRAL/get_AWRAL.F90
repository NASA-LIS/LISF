!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
!    integer :: endtime_AWRAL         ! 1=get a new file 
    integer :: order

    real*8  :: timenext
    real*8  :: AWRAL_file_time1       ! End boundary time for STAGEIV file
    real*8  :: AWRAL_file_time2       ! End boundary time for STAGEIV file
    character(80) :: file_name       ! Filename variables for precip data sources

    integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
    integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
    real    :: gmt1, gmt2,ts1,ts2                    
    real    :: gridDesci(LIS_rc%nnest,50)

!=== End Variable Definition =======================

    yr1 = LIS_rc%yr  !current time
    mo1 = LIS_rc%mo
    da1 = LIS_rc%da
    hr1 = 0
    mn1 = 0
    ss1 = 0
    ts1 = 86400
    call LIS_tick( timenext, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

!-- Determine LIS's current time and the time of the AWRAL file:
    yr1 = LIS_rc%yr
    mo1 = LIS_rc%mo
    da1 = LIS_rc%da
    hr1 = 0
    mn1 = 0
    ss1 = 0
    ts1 = -86400
    call LIS_tick( AWRAL_file_time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
!-- AWRAL product time; end accumulation time data
    yr2 = LIS_rc%yr     
    mo2 = LIS_rc%mo
    da2 = LIS_rc%da
    hr2 = 0
    mn2 = 0
    ss2 = 0
    ts2 = 0
    call LIS_tick( AWRAL_file_time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!-- Ensure that data is found during first time step
    if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then 
!         endtime_AWRAL = 1
         LIS_rc%rstflag(n) = 0
    endif

!-- Check for and get AWRAL Precipitation data
    ferror_AWRAL = 0
    order = 2

  ! LIS timestep < AWRAL time interval (1hr)
    if( LIS_rc%ts < AWRAL_struc(n)%ts ) then
      if ( LIS_rc%time > AWRAL_struc(n)%AWRALtime ) then
      ! Determine and return filename of AWRAL file 
        call AWRALfile( file_name, AWRAL_struc(n)%AWRALdir, yr2, doy2)
        write(LIS_logunit,*) '[INFO] Getting new AWRAL precip data:: ', file_name
      ! Open, read, and reinterpolate AWRAL field to LIS-defined grid
        call read_AWRAL ( n, file_name, findex, order, ferror_AWRAL )
      ! Assign latest AWRAL file time to stored AWRAL time variable
        AWRAL_struc(n)%AWRALtime = timenext
      endif

  ! LIS timestep == AWRAL time interval (1hr)
    elseif( LIS_rc%ts == AWRAL_struc(n)%ts ) then

     ! Determine and return filename of AWRAL file 
       call AWRALfile( file_name, AWRAL_struc(n)%AWRALdir, yr1, doy1 )
       write(LIS_logunit,*) '[INFO] Getting new AWRAL precip data:: ', file_name
     ! Open, read, and reinterpolate AWRAL field to LIS-defined grid
       call read_AWRAL ( n, file_name, findex, order, ferror_AWRAL )
     ! Assign latest AWRAL file time to stored AWRAL time variable
       AWRAL_struc(n)%AWRALtime = AWRAL_file_time1

    else
      write(LIS_logunit,*) "[ERR] AWRAL READER CANNOT HANDLE LIS "
      write(LIS_logunit,*) "[ERR] TIMESTEP > 3600 secs -- AT THIS TIME.  STOPPING ..."
      call LIS_endrun
    endif

end subroutine get_AWRAL

