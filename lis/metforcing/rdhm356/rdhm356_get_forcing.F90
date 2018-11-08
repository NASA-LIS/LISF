!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: getrdhm356
! \label{getrdhm356}
!
! !REVISION HISTORY:
! 25 May 2006: Kristi Arsenault;  Data and Code implementation
! 18 Dec 2013: Shugong Wang; implementation for RDHM356
! 
! !INTERFACE:
subroutine get_rdhm356(n,findex)

! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick, LIS_get_nstep
  use LIS_logMod, only : LIS_logunit
  use rdhm356_forcingMod, only : rdhm356_struc_precip, &
                                 rdhm356_struc_temper


  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n 
  integer, intent(in) :: findex

  
! !DESCRIPTION:
!  Opens, reads, and interpolates hourly rdhm356 forcing. 
!  At the beginning of a simulation, the code reads the most
!  recent past data (nearest the hour interval), and the nearest
!  future data. These two datasets are used to temporally 
!  interpolate the data to the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[timenow]
!    Current LIS Time 
!  \item[rdhm356\_file\_time2]
!    End boundary time of dmip II file 
!  \item[file\_name]
!    forcing file filename - passed back to getstg 
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the RDHM data times
!  \item[rdhm356\_precip\_file](\ref{rdhm356_precip_file}) \newline
!    Puts together appropriate file name for 1-hour intervals
!  \item[rhdm356\_read\_precip](\ref{rdhm356_read_precip}) \newline
!    Interpolates rdhm356 data to LIS grid
!  \end{description}
!EOP
   
    !== Local Variables =======================
    integer :: ferror_rdhm356                   ! Error flags for precip data sources
    integer :: endtime_rdhm356                  ! 1=get a new file 

    real*8  :: timenow, timenowback, rdhm356_file_time1, rdhm356_file_time2      ! Current LIS Time end boundary time and file1 file2 times for dmipII file
    real*8 :: rdhm356_file_time3
    character(128) :: file_name               ! Filename variables for precip data sources
    integer :: order

    integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
    integer :: doy0, yr0, mo0, da0, hr0, mn0, ss0
    integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
    integer :: doy3, yr3, mo3, da3, hr3, mn3, ss3
    real    :: gmt1, gmt2,gmt0,gmt3      
    real    :: ts1, ts0, ts2, ts3
                                             !  for precip data sources
    integer c,f
    !=== End Variable Definition =======================

    endtime_rdhm356 = 0
    !-- Determine LIS's current time and the time of the dmip 2 file:
    ! - LIS's current date and time
    ! put current time in working variables
    yr1 = LIS_rc%yr  
    mo1 = LIS_rc%mo
    da1 = LIS_rc%da
    hr1 = LIS_rc%hr
    mn1 = LIS_rc%mn
    ss1 = 0
    ts1 = 0      ! Time increment in seconds
    call LIS_tick( timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )


    ! store current time for use as is
    yr0 = LIS_rc%yr
    mo0 = LIS_rc%mo
    da0 = LIS_rc%da
    hr0 = LIS_rc%hr
    mn0 = LIS_rc%mn
    ss0 = 0
    ts0 = -3600      ! Time increment in seconds

    call LIS_tick( timenowback, doy0, gmt0, yr0, mo0, da0, hr0, mn0, ss0, ts0 )


    !-- Determine LIS's current time and the time of the dmip 2 file:
    ! - LIS's current date and time
    yr1 = LIS_rc%yr
    mo1 = LIS_rc%mo
    da1 = LIS_rc%da
    hr1 = 1*(int(real(LIS_rc%hr)/1.0))
    mn1 = 0
    ss1 = 0
    ts1 = 0      ! Time increment in seconds

    call LIS_tick( rdhm356_file_time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )

    !-- rdhm356 product time; end accumulation time data
    yr2 = LIS_rc%yr     
    mo2 = LIS_rc%mo
    da2 = LIS_rc%da
    !    hr2 = LIS_rc%hr+1  ! Advance dmipII by one hour increment
    hr2 = LIS_rc%hr
    mn2 = 0
    ss2 = 0
    ts2= 3600

    call LIS_tick( rdhm356_file_time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2 )

    yr3 = LIS_rc%yr
    mo3 = LIS_rc%mo
    da3 = LIS_rc%da
    hr3 = LIS_rc%hr
    mn3 = 0
    ss3 = 0
    ts3= 7200

    call LIS_tick( rdhm356_file_time3, doy3, gmt3, yr3, mo3, da3, hr3, mn3, ss3, ts3 )

    write(LIS_logunit,*) 'LIS_get_nstep is ', LIS_get_nstep(LIS_rc,n)


    !-- Check for and get precipitation data
    if ( timenow >= rdhm356_struc_precip(n)%rdhm356time2 ) then
        endtime_rdhm356 = 1
    endif 

    !-- Get new second file for both temperature and precip
    ! transfer time 2 temperature data into time 1 if not on first timestep
    if ( endtime_rdhm356 == 1 ) then  
       rdhm356_struc_precip(n)%rdhm356time1 = rdhm356_file_time1 
       rdhm356_struc_temper(n)%rdhm356time1 = rdhm356_file_time1
    endif
    
    ! reset restart flag
    LIS_rc%rstflag=0
    if ( endtime_rdhm356 == 1 ) then
       !   -- Determine and return filename of second dmip II precip file 
       call rdhm356_precip_file(file_name, &
                                rdhm356_struc_precip(n)%rdhm356dir, &
                                yr1, mo1, da1, hr1 )


       write(LIS_logunit,*) 'Getting new rdhm356 precip data: ', trim(file_name)
       
       ferror_rdhm356 = 0
       order = 2
       call rdhm356_read_precip ( n, trim(file_name), findex, order, ferror_rdhm356 )


       call rdhm356_temper_file(file_name, &
                           rdhm356_struc_temper(n)%rdhm356dir, &
                           yr1, mo1, da1, hr1 ) ! SY: Changed for Comp. W gorms

       write(LIS_logunit,*) 'Getting new rdhm356 temp data: ', trim(file_name)
      
       ferror_rdhm356 = 0
       order = 2
       call rdhm356_read_temper ( n, findex, order, trim(file_name), ferror_rdhm356 )


!   -- Assign latest dmip II file time to stored dmip II time variable
       rdhm356_struc_precip(n)%rdhm356time2 = rdhm356_file_time2
       rdhm356_struc_temper(n)%rdhm356time2 = rdhm356_file_time2

     endif  
  return

end subroutine get_rdhm356
