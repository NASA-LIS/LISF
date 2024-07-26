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
  use LDT_coreMod,    only  : LDT_rc, LDT_domain
  use LDT_logMod,     only  : LDT_logunit, LDT_endrun
  use LDT_timeMgrMod, only  : LDT_tick, LDT_get_nstep, LDT_date2time
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
!  \item[ldt\_time]
!    Current LDT Time 
!  \item[stg4\_file\_time]
!    End boundary time of Stage IV file 
!  \item[file\_name]
!    Stage IV filename - passed back to getstg 
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    determines the STAGE4 data times
!  \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    Computes the neighbor, weights for bilinear interpolation
!  \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    Computes the neighbor, weights for conservative interpolation
!  \item[stg4file](\ref{stg4file}) \newline
!    Puts together appropriate file name for 1-hour intervals
!  \item[read\_stg4](\ref{read_stg4}) \newline
!    Interpolates STAGE4 data to LDT grid
!  \end{description}
!EOP
   
!== Local Variables =======================
    integer :: ferror_stg4               ! Error flags for precip data sources
    integer :: order
    real*8  :: stg4_file_time            ! Current end boundary time for STAGEIV file
    character(80) :: file_name           ! Filename variables for precip data sources

    integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
    integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
    real    :: gmt1, gmt2,ts1,ts2                    
    integer :: updoy
    real*8  :: start_time

    real    :: gridDesci(LDT_rc%nnest,20)

!=== End Variable Definition =======================

  ! Check for when first file date/time occurs:
    yr1 = 2002
    mo1 = 01
    da1 = 01
    hr1 = 00; mn1 = 0; ss1 = 0
    call LDT_date2time( start_time, updoy, gmt1,&
                        yr1, mo1, da1, hr1, mn1, ss1 )
    if( start_time > LDT_rc%time ) then
       write(LDT_logunit,*)" MSG: The Current LDT time is before the "
       write(LDT_logunit,*)"   STAGE IV start date and time. Stopping ..."
       call LDT_endrun
    endif

!-- STAGE4 start of product-accumulation time
    yr1 = LDT_rc%yr
    mo1 = LDT_rc%mo
    da1 = LDT_rc%da
    hr1 = LDT_rc%hr
    mn1 = 0
    ss1 = 0
    ts1 = 0
    call LDT_tick( stg4_struc(n)%stg4time, doy1, gmt1, &
                   yr1, mo1, da1, hr1, mn1, ss1, ts1 )

!-- STAGE4 product time; end accumulation time data
    yr2 = LDT_rc%yr
    mo2 = LDT_rc%mo
    da2 = LDT_rc%da
!    hr2 = LDT_rc%hr+1  ! Advance STAGEIV by one hour increment
    hr2 = LDT_rc%hr    ! Advance STAGEIV by one hour increment
    mn2 = 3600*(LDT_rc%mn/3600)
    mn2 = 0
    ss2 = 0
    ts2 = 0
    call LDT_tick( stg4_file_time, doy2, gmt2, &
                   yr2, mo2, da2, hr2, mn2, ss2, ts2 )

!-- Determine if STAGE IV Grid LDT_domain has changed (before opening file)   
    if( LDT_rc%time > stg4_struc(n)%griduptime1 .and. &
       stg4_struc(n)%gridchange1 ) then 

       write(LDT_logunit,*) "** Changing STAGE4 Grid to 2004 -" 

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

       select case( LDT_rc%met_gridtransform(findex) )
         case( "bilinear" )
         ! = Bilinear interpolation =
          call bilinear_interp_input( n, gridDesci(n,:),&
               stg4_struc(n)%n111, stg4_struc(n)%n121,  &
               stg4_struc(n)%n211, stg4_struc(n)%n221,  &
               stg4_struc(n)%w111, stg4_struc(n)%w121,  &
               stg4_struc(n)%w211, stg4_struc(n)%w221 )

         case( "budget-bilinear" )
         ! = Budget-bilinear interpolation =
          call conserv_interp_input( n, gridDesci(n,:),  &
                 stg4_struc(n)%n112, stg4_struc(n)%n122, &
                 stg4_struc(n)%n212, stg4_struc(n)%n222, &
                 stg4_struc(n)%w112, stg4_struc(n)%w122, &
                 stg4_struc(n)%w212, stg4_struc(n)%w222 )
       end select

       stg4_struc(n)%gridchange1 = .false.   ! Grid change has been applied
    endif

!-- Ensure that data is found during first time step
    if( LDT_get_nstep(LDT_rc,n) == 1 .or. LDT_rc%rstflag(n) == 1 ) then 
       LDT_rc%rstflag(n) = 0
    endif

!-- Check for and get STAGE4 Precipitation data
    ferror_stg4 = 0
    order = 2
    if( LDT_rc%ts >= stg4_struc(n)%ts ) then
       if ( LDT_rc%time == stg4_struc(n)%stg4time ) then
       ! Determine and return filename of STAGE IV file 
         call stg4file( file_name, stg4_struc(n)%stg4dir, yr2, mo2, da2, hr2 )
         write(LDT_logunit,*) "Getting new STAGE4 precip data:: ",trim(file_name)
       ! Open, read, and reinterpolate STAGE IV field to LDT-defined grid
         call read_stg4 ( n, file_name, findex, order, ferror_stg4 )
      end if
    endif  

end subroutine get_stg4

