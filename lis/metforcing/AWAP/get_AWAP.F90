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
! !ROUTINE: get_AWAP
! \label{get_AWAP}
!
! !REVISION HISTORY:
! 30 Jan 2017: Sujay Kumar, Initial version
!
! !INTERFACE:
subroutine get_AWAP(n, findex)

! !USES:
  use LIS_coreMod, only     : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only  : LIS_tick, LIS_get_nstep
  use LIS_logMod,      only : LIS_logunit, LIS_endrun
  use AWAP_forcingMod, only : AWAP_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n 
  integer, intent(in) :: findex

! !DESCRIPTION:
!  Opens, reads, and interpolates AWAP forcing. 
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
!    determines the AWAP data times
!  \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    Computes the neighbor, weights for bilinear interpolation
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    Computes the neighbor, weights for conservative interpolation
!  \item[AWAPfile](\ref{AWAPfile}) \newline
!    Puts together appropriate file name for 1-hour intervals
!  \item[read\_AWAP](\ref{read_AWAP}) \newline
!    Interpolates AWAP data to LIS grid
!  \end{description}
!EOP
   
!== Local Variables =======================
    integer :: ferror_AWAP           ! Error flags for precip data sources
!    integer :: endtime_AWAP         ! 1=get a new file 
    integer :: order

    real*8  :: timenext
    real*8  :: AWAP_file_time1       ! End boundary time for STAGEIV file
    real*8  :: AWAP_file_time2       ! End boundary time for STAGEIV file
    character(len=LIS_CONST_PATH_LEN) :: file_name ! Filename variables for precip data sources

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

!-- Determine LIS's current time and the time of the AWAP file:
    yr1 = LIS_rc%yr
    mo1 = LIS_rc%mo
    da1 = LIS_rc%da
    hr1 = 0
    mn1 = 0
    ss1 = 0
    ts1 = -86400
    call LIS_tick( AWAP_file_time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )
!-- AWAP product time; end accumulation time data
    yr2 = LIS_rc%yr     
    mo2 = LIS_rc%mo
    da2 = LIS_rc%da
    hr2 = 0
    mn2 = 0
    ss2 = 0
    ts2 = 0
    call LIS_tick( AWAP_file_time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!-- Ensure that data is found during first time step
    if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then 
!         endtime_AWAP = 1
         LIS_rc%rstflag(n) = 0
    endif

!-- Check for and get AWAP Precipitation data
    ferror_AWAP = 0
    order = 2

  ! LIS timestep < AWAP time interval (1hr)
    if( LIS_rc%ts < AWAP_struc(n)%ts ) then
      if ( LIS_rc%time > AWAP_struc(n)%AWAPtime ) then
      ! Determine and return filename of AWAP file 
        call AWAPfile( file_name, AWAP_struc(n)%AWAPdir, yr2, doy2)
        write(LIS_logunit,*) '[INFO] Getting new AWAP precip data:: ', trim(file_name)
      ! Open, read, and reinterpolate AWAP field to LIS-defined grid
        call read_AWAP ( n, file_name, findex, order, ferror_AWAP )
      ! Assign latest AWAP file time to stored AWAP time variable
        AWAP_struc(n)%AWAPtime = timenext
      endif

  ! LIS timestep == AWAP time interval (1hr)
    elseif( LIS_rc%ts == AWAP_struc(n)%ts ) then

     ! Determine and return filename of AWAP file 
       call AWAPfile( file_name, AWAP_struc(n)%AWAPdir, yr1, doy1 )
       write(LIS_logunit,*) '[INFO] Getting new AWAP precip data:: ', trim(file_name)
     ! Open, read, and reinterpolate AWAP field to LIS-defined grid
       call read_AWAP ( n, file_name, findex, order, ferror_AWAP )
     ! Assign latest AWAP file time to stored AWAP time variable
       AWAP_struc(n)%AWAPtime = AWAP_file_time1

    else
      write(LIS_logunit,*) "[ERR] AWAP READER CANNOT HANDLE LIS "
      write(LIS_logunit,*) "[ERR] TIMESTEP > 3600 secs -- AT THIS TIME.  STOPPING ..."
      call LIS_endrun
    endif

end subroutine get_AWAP

