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
! !ROUTINE: get_stg4
! \label{get_stg4}
!
! !REVISION HISTORY:
! 25 May 2006: Kristi Arsenault;  Data and Code implementation
!
! !INTERFACE:
subroutine get_stg4(n, findex)

! !USES:
  use LIS_coreMod, only     : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only  : LIS_tick, LIS_get_nstep
  use LIS_logMod,      only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,only : LIS_CONST_PATH_LEN
  use stg4_forcingMod, only : stg4_struc

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n 
  integer, intent(in) :: findex

! !DESCRIPTION:
!  Opens, reads, and interpolates hourly STAGE4 forcing. 
!  At the beginning of a simulation, the code reads the most
!  recent past data (nearest the hour interval), and the nearest
!  future data. These two datasets are used to temporally 
!  interpolate the data to the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[lis\_time]
!    Current LIS Time 
!  \item[stg4\_file\_time]
!    End boundary time of Stage IV file 
!  \item[file\_name]
!    Stage IV filename - passed back to getstg 
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the STAGE4 data times
!  \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    Computes the neighbor, weights for bilinear interpolation
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    Computes the neighbor, weights for conservative interpolation
!  \item[stg4file](\ref{stg4file}) \newline
!    Puts together appropriate file name for 1-hour intervals
!  \item[read\_stg4](\ref{read_stg4}) \newline
!    Interpolates STAGE4 data to LIS grid
!  \end{description}
!EOP
   
!== Local Variables =======================
    integer :: ferror_stg4           ! Error flags for precip data sources
!    integer :: endtime_stg4         ! 1=get a new file 
    integer :: order

    real*8  :: stg4_file_time1       ! End boundary time for STAGEIV file
    real*8  :: stg4_file_time2       ! End boundary time for STAGEIV file
    character(len=LIS_CONST_PATH_LEN) :: file_name       ! Filename variables for precip data sources

    integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
    integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
    real    :: gmt1, gmt2,ts1,ts2                    
    real    :: gridDesci(LIS_rc%nnest,50)

!=== End Variable Definition =======================

!    endtime_stg4 = 0

!-- Determine LIS's current time and the time of the STAGE 4 file:
    yr1 = LIS_rc%yr
    mo1 = LIS_rc%mo
    da1 = LIS_rc%da
    hr1 = LIS_rc%hr    
    mn1 = 0
    ss1 = 0
    ts1 = 0
    call LIS_tick( stg4_file_time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

!-- STAGE4 product time; end accumulation time data
    yr2 = LIS_rc%yr     
    mo2 = LIS_rc%mo
    da2 = LIS_rc%da
    hr2 = LIS_rc%hr+1  ! Advance STAGEIV by one hour increment
    mn2 = 0
    ss2 = 0
    ts2 = 0
    call LIS_tick( stg4_file_time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!-- Determine if STAGE IV Grid LIS_domain has changed (before opening file)   
    if( LIS_rc%time > stg4_struc(n)%griduptime1 .and. &
       stg4_struc(n)%gridchange1 ) then 

       write(LIS_logunit,*) "** Changing STAGE4 Grid to 2004 -" 

!   -- Reinitialize the weights and neighbors
       gridDesci(n,:) = 0.0
       gridDesci(n,1) = 5.0                 ! Projection type (UPS)
       gridDesci(n,2) = stg4_struc(n)%ncol  ! X-dir amount of points
       gridDesci(n,3) = stg4_struc(n)%nrow  ! y-dir amount of points
       gridDesci(n,4) = 23.11700            ! Starting latitude point
       gridDesci(n,5) = -119.02300          ! Starting longitude point
       gridDesci(n,6) = 8.0
       gridDesci(n,7) = 0.0                 ! Orientation
       gridDesci(n,8) = 4.7625              ! X-spacing length (meters) 
       gridDesci(n,9) = 4.7625              ! Y-spacing length (meters)
       gridDesci(n,10) = 60.0               ! True lat 
       gridDesci(n,11) = -105.0             ! Standard longitude (orientation)
       gridDesci(n,20) = 64.0               ! ??
!       gridDesci(n,20) = 0.

       select case( LIS_rc%met_interp(findex) )
        case( "bilinear" )
          call bilinear_interp_input( n, gridDesci(n,:),&
               stg4_struc(n)%n111, stg4_struc(n)%n121,  &
               stg4_struc(n)%n211, stg4_struc(n)%n221,  &
               stg4_struc(n)%w111, stg4_struc(n)%w121,  &
               stg4_struc(n)%w211, stg4_struc(n)%w221 )

        case( "budget-bilinear" )
          call conserv_interp_input( n, gridDesci(n,:),  &
                 stg4_struc(n)%n112, stg4_struc(n)%n122, &
                 stg4_struc(n)%n212, stg4_struc(n)%n222, &
                 stg4_struc(n)%w112, stg4_struc(n)%w122, &
                 stg4_struc(n)%w212, stg4_struc(n)%w222 )
       end select
       stg4_struc(n)%gridchange1 = .false.   ! Grid change has been applied
    endif

!-- Ensure that data is found during first time step
    if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then 
!         endtime_stg4 = 1
         LIS_rc%rstflag(n) = 0
    endif

!-- Check for and get STAGE4 Precipitation data
    ferror_stg4 = 0
    order = 2

  ! LIS timestep < STAGEIV time interval (1hr)
    if( LIS_rc%ts < stg4_struc(n)%ts ) then
      if ( LIS_rc%time > stg4_struc(n)%stg4time ) then

      ! Determine and return filename of STAGE IV file 
        call stg4file( file_name, stg4_struc(n)%stg4dir, yr2, mo2, da2, hr2 )
        write(LIS_logunit,*) 'Getting new STAGE4 precip data:: ', trim(file_name)
      ! Open, read, and reinterpolate STAGE IV field to LIS-defined grid
        call read_stg4 ( n, file_name, findex, order, ferror_stg4 )
      ! Assign latest STAGE IV file time to stored STAGE IV time variable
        stg4_struc(n)%stg4time = stg4_file_time2
      endif

  ! LIS timestep == STAGEIV time interval (1hr)
    elseif( LIS_rc%ts == stg4_struc(n)%ts ) then

     ! Determine and return filename of STAGE IV file 
       call stg4file( file_name, stg4_struc(n)%stg4dir, yr1, mo1, da1, hr1 )
       write(LIS_logunit,*) 'Getting new STAGE4 precip data:: ', trim(file_name)
     ! Open, read, and reinterpolate STAGE IV field to LIS-defined grid
       call read_stg4 ( n, file_name, findex, order, ferror_stg4 )
     ! Assign latest STAGE IV file time to stored STAGE IV time variable
       stg4_struc(n)%stg4time = stg4_file_time1

    else
      write(LIS_logunit,*) "ERR MSG: STAGE IV READER CANNOT HANDLE LIS "
      write(LIS_logunit,*) "  TIMESTEP > 3600 secs -- AT THIS TIME.  STOPPING ..."
      call LIS_endrun
    endif

end subroutine get_stg4

