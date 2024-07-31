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
! !ROUTINE: get_TRMM3B42RT
! \label{get_TRMM3B42RT}
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
!
! !INTERFACE:
subroutine get_TRMM3B42RT(n,findex)
! !USES:
  use LIS_coreMod, only           : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod,  only       : LIS_time2date, LIS_tick, LIS_get_nstep, &
                                    LIS_isAlarmRinging ! SY
  use LIS_logMod, only            : LIS_logunit, LIS_endrun
  use TRMM3B42RT_forcingMod, only : TRMM3B42RT_struc
  use LIS_metforcingMod,  only    : LIS_forc ! SY
  use LIS_constantsMod,      only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, TRMM 3B42RT forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the TRMM 3B42RT data times
!  \item[TRMM3B42RTfile](\ref{TRMM3B42RTfile}) \newline
!    Puts together appropriate file name for 3 hour intervals
!  \item[read\_TRMM3B42RT](\ref{read_TRMM3B42RT}) \newline
!      Interpolates TRMM 3B42RT data to LIS grid
!  \end{description}
!EOP

   
!==== Local Variables=======================
  integer :: ferror_TRMM3B42RT            ! Error flag for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1,ts1     ! SY: Time parameters for start and end of current LDAS time step
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2    ! SY: Time parameters for TRMM data time nearest to start of model time step
  integer :: doy3, yr3, mo3, da3, hr3, mn3, ss3, ts3    ! SY: Time parameters for TRMM data time nearest to end of model time step
  real    :: gmt1, gmt2, gmt3
  character(len=LIS_CONST_PATH_LEN) :: name                    ! Filename variables for precip data sources
  real*8  :: LIS_timeAtTStepStart_add90min 
  real*8  :: LIS_timeAtTStepEnd_add90min   
  integer :: order
  logical :: alarmCheck 

!=== End Variable Definition =======================

!  Disable alarm.  It throws off the timing of reading the
!  next data set, and it throws off copying the bookends.
!  alarmCheck = LIS_isAlarmRinging(LIS_rc,"TRMM 3B42RT alarm") ! SY
!  if(alarmCheck) then 

 !------------------------------------------------------------------------
 ! SY : Determine start and end time of model time step
 !------------------------------------------------------------------------
  yr1 = LIS_rc%yr  !current time, ! SY: i.e., end of model time step
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = LIS_rc%ss
  ts1 = 0
  call LIS_tick( TRMM3B42RT_struc(n)%LIS_timeAtTStepEnd, doy1, gmt1, &
       yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
  ts1 = (-1)*LIS_rc%nts(n) ! SY: for start of model time step
  call LIS_tick( TRMM3B42RT_struc(n)%LIS_timeAtTStepStart, doy1, gmt1, &
       yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
 
 !------------------------------------------------------------------------
 ! SY: 3B42RT product times.
 !     I get these times based on a very simple principle: modify the
 !     time step start/end by adding 90 min. Whichever TRMM data stamp the
 !     modified start/end coincides with or immediately comes after, will be
 !     the TRMM data nearest to the original time step start/end. Hence,
 !     values of that TRMM data are assigned to that original start/end as
 !     preparation for weighting data in the time_interp* subroutine.
 !------------------------------------------------------------------------

   ! SY: add 90 min to start of model time step
   yr2 = LIS_rc%yr
   mo2 = LIS_rc%mo
   da2 = LIS_rc%da
   hr2 = LIS_rc%hr
   mn2 = LIS_rc%mn
   ss2 = LIS_rc%ss
   ts2 = (-1)*LIS_rc%nts(n) + 90*60
   call LIS_tick( LIS_timeAtTStepStart_add90min, doy2, gmt2, &
        yr2, mo2, da2, hr2, mn2, ss2, real(ts2))
   ! SY: Now start calculations for TRMM data time nearest to start of model time step
   TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart = yr2
   TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart = mo2
   TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart = da2
   TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart = 3*(hr2/3)
   mn2 = 0
   ss2 = 0
   ts2 = 0
   call LIS_tick( TRMM3B42RT_struc(n)%TRMM3B42RTtime_TStepStart, doy2, gmt2, &
        TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart, &
        TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart, &
        TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart, &
        TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart, mn2, ss2, real(ts2))
 
   ! SY: add 90 min to end of model time step
   yr3 = LIS_rc%yr
   mo3 = LIS_rc%mo
   da3 = LIS_rc%da
   hr3 = LIS_rc%hr
   mn3 = LIS_rc%mn
   ss3 = LIS_rc%ss
   ts3 = 90*60
   call LIS_tick( LIS_timeAtTStepEnd_add90min, doy3, gmt3, &
        yr3, mo3, da3, hr3, mn3, ss3, real(ts3))
   ! SY: Now start calculations for TRMM data time nearest to end of model time step
   TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd = yr3
   TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd = mo3
   TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd = da3
   TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd = 3*(hr3/3)
   mn3 = 0
   ss3 = 0
   ts3 = 0
   call LIS_tick( TRMM3B42RT_struc(n)%TRMM3B42RTtime_TStepEnd, doy3, gmt3, &
        TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd, &
        TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd, &
        TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd, &
        TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd, mn3, ss3, real(ts3))
 
 !------------------------------------------------------------------------
 ! Check for and get 3B42RT observed Precipitation data
 !------------------------------------------------------------------------
 
   ! SY: Now update TRMM data corresponding to start of model time step if required
   if  (LIS_get_nstep(LIS_rc, n).eq. 1 .or. LIS_rc%rstflag(n) .eq. 1) then

      if (LIS_rc%nts(n) .ge. (3*60-1)*60) then ! almost 3 hrs
         write(LIS_logunit,*) 'LIS time step should be < 3 hrs!'
         write(LIS_logunit,*) 'TRMM reader functionality for >= 3 hrs not written yet'
         write(LIS_logunit,*) 'will involve weights combining > 2 TRMM data files'
         write(LIS_logunit,*) 'Program stopping ... '
         call LIS_endrun()
      endif

     TRMM3B42RT_struc(n)%metdata1 = LIS_rc%udef
     TRMM3B42RT_struc(n)%metdata2 = LIS_rc%udef
     ferror_TRMM3B42RT = 0
     write(LIS_logunit, *) 'TRMM yr,mo,da,hr for time step start:', &
                           TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart, &
                           TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart, &
                           TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart, &
                           TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart
     call TRMM3B42RTfile( name, TRMM3B42RT_struc(n)%TRMM3B42RTdir,     &
                          TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart, &
                          TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart, &
                          TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart, &
                          TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart )
     write(LIS_logunit, *)'Getting new TRMM 3B42RT satellite precip data:', trim(name)
     order = 1
     call read_TRMM3B42RT(n, name, findex, order, ferror_TRMM3B42RT)
   elseif (.NOT. ((TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart_Previous) .AND. &
              (TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart_Previous) .AND. &
              (TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart_Previous) .AND. &
              (TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart_Previous))) then
     write(LIS_logunit, *) 'TRMM yr,mo,da,hr for time step start:', &
                           TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart, &
                           TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart, &
                           TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart, &
                           TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart
     write(LIS_logunit, *)'Values from time step end assigned'
     TRMM3B42RT_struc(n)%metdata1 = TRMM3B42RT_struc(n)%metdata2 
     ! SY: i.e., equal to TRMM3B42RT_struc(n)%metdata2 before possible change in TRMM3B42RT_struc(n)%metdata2 below
   endif

   ! SY: Now update TRMM data corresponding to end of model time step if required.
   if (.NOT. ((TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd_Previous) .AND. &
              (TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd_Previous) .AND. &
              (TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd_Previous) .AND. &
              (TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd .EQ. &
               TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd_Previous))) then
     ferror_TRMM3B42RT = 0
     write(LIS_logunit, *) 'TRMM yr,mo,da,hr for time step end:', &
                           TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd, &
                           TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd, &
                           TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd, &
                           TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd
     call TRMM3B42RTfile( name, TRMM3B42RT_struc(n)%TRMM3B42RTdir,   &
                          TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd, &
                          TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd, &
                          TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd, &
                          TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd )
     write(LIS_logunit, *)'Getting new TRMM 3B42RT satellite precip data:', trim(name)
     order = 2
     call read_TRMM3B42RT(n, name, findex, order, ferror_TRMM3B42RT)
   endif

   ! SY: Begin reassigning *_Previous values for use in next get* subroutine call
   TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepStart
   TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepStart
   TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepStart
   TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepStart
 
   TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RTyr_TStepEnd
   TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RTmo_TStepEnd
   TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RTda_TStepEnd
   TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd_Previous = &
                              TRMM3B42RT_struc(n)%TRMM3B42RThr_TStepEnd
   ! SY: End reassigning *_Previous values for use in next get* subroutine call
 
 ! J.Case (4/19/2013) -- test print of the suppdata2 array.
 ! write (96,*) suppdata2(1,:)
!  end if ! SY: if(alarmCheck) then

end subroutine get_TRMM3B42RT

!BOP
! !ROUTINE: TRMM3B42RTfile
! \label{TRMM3B42RTfile}
!
!
! !INTERFACE:
subroutine TRMM3B42RTfile( name, TRMM3B42RTdir, yr, mo, da, hr)
! !USES: 
  use TRMM3B42RT_forcingMod, only : TRMM3B42RT_struc
  use LIS_timeMgrMod, only : LIS_date2time
  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name
  character(len=*)  :: TRMM3B42RTdir
  integer           :: yr, mo, da, hr
!
! !DESCRIPTION:
!   This subroutine puts together TRMM 3B42RT file name for 
!   3 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[TRMM3B42RTdir]
!    Name of the TRMM 3B42RT data directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[name]
!   name of the timestamped TRMM 3B42RT file
!  \end{description}
!
!EOP

  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1, udoy
  real    :: ugmt
  real*8  :: utime
  character*2 :: m2, d2, h2 
  character*4 :: y4

  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0

  write(y4, '(I4.4)') uyr 
  write(m2, '(I2.2)') umo
  write(d2, '(I2.2)') uda 
  write(h2, '(I2.2)') uhr

  call LIS_date2time(utime, udoy, ugmt, uyr, umo, uda, uhr, umn, uss ) 
  if (utime .LT. TRMM3B42RT_struc(1)%griduptime1) then   ! RT V5
    name = trim(TRMM3B42RTdir)//"/"// y4 // m2 // "/3B42RT."//y4//m2//d2//h2 
  else
    name = trim(TRMM3B42RTdir)//"/"// y4 // m2 // "/3B42RT."//y4//m2//d2//h2//".6" 
  end if 

end subroutine TRMM3B42RTfile
