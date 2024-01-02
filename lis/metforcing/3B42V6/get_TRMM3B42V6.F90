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
! !ROUTINE: get_TRMM3B42V6
! \label{get_TRMM3B42V6}
!
! !REVISION HISTORY:
! 25Aug2006: Yudong Tian; Modified for LIS 4.2 release 
! 21 Jun 2013: Soni Yatheendradas; changes from earlier code to avoid
!              (a) alternate file skip,
!              (b) jump to previous day TRMM, and
!              (c) absence of rain rate weighting
!
! !INTERFACE:
subroutine get_TRMM3B42V6(n, findex)
! !USES:
  use LIS_coreMod, only           : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod, only        : LIS_time2date, LIS_tick, LIS_get_nstep, &
                                    LIS_isAlarmRinging ! SY
  use LIS_logMod, only            : LIS_logunit, LIS_endrun
  use TRMM3B42V6_forcingMod, only : TRMM3B42V6_struc
  use LIS_metforcingMod, only     : LIS_forc ! SY
  use LIS_constantsMod,      only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly TRMM 3B42V6 forcing.
!  At the beginning of a simulation, the code
!  reads the most recent past data (nearest 6 hour interval), and
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
!    determines the TRMM 3B42V6 data times
!  \item[TRMM3B42V6file](\ref{TRMM3B42V6file}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[read\_TRMM3B42V6](\ref{read_TRMM3B42V6}) \newline
!      Interpolates TRMM 3B42V6 data to LIS grid
!  \end{description}
!EOP

!==== Local Variables=======================
  integer :: ferror_TRMM3B42V6      ! Error flag for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1,ts1                ! SY: Time parameters for start and end of current LDAS time step
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2               ! SY: Time parameters for TRMM data time nearest to start of model time step
  integer :: doy3, yr3, mo3, da3, hr3, mn3, ss3, ts3               ! SY: Time parameters for TRMM data time nearest to end of model time step
  real    :: gmt1, gmt2, gmt3 ! SY ,kgmt3, mgmt3
  character(len=LIS_CONST_PATH_LEN) :: name ! Filename variables for precip data sources
  real*8 :: LIS_timeAtTStepStart_add90min ! SY
  real*8 :: LIS_timeAtTStepEnd_add90min ! SY
  integer :: order
  logical             :: alarmCheck ! SY

!=== End Variable Definition =======================

!  Disable alarm.  It throws off the timing of reading the
!  next data set, and it throws off copying the bookends.
!  alarmCheck = LIS_isAlarmRinging(LIS_rc,"TRMM 3B42V6 alarm") ! SY
!  if(alarmCheck) then ! SY

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
   call LIS_tick( TRMM3B42V6_struc(n)%LIS_timeAtTStepEnd, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
   ts1 = (-1)*LIS_rc%nts(n) ! SY: for start of model time step
   call LIS_tick( TRMM3B42V6_struc(n)%LIS_timeAtTStepStart, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
 
 !------------------------------------------------------------------------
 ! SY: 3B42V6 product times.
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
   TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart = yr2
   TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart = mo2
   TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart = da2
   TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart = 3*(hr2/3)
   mn2 = 0
   ss2 = 0
   ts2 = 0
   call LIS_tick( TRMM3B42V6_struc(n)%TRMM3B42V6time_TStepStart, doy2, gmt2, &
        TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart, &
        TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart, &
        TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart, &
        TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart, mn2, ss2, real(ts2))
 
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
   TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd = yr3
   TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd = mo3
   TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd = da3
   TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd = 3*(hr3/3)
   mn3 = 0
   ss3 = 0
   ts3 = 0
   call LIS_tick( TRMM3B42V6_struc(n)%TRMM3B42V6time_TStepEnd, doy3, gmt3, &
        TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd, &
        TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd, &
        TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd, &
        TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd, mn3, ss3, real(ts3))
 
 !------------------------------------------------------------------------
 ! Check for and get 3B42V6 observed Precipitation data
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

     TRMM3B42V6_struc(n)%metdata1 = LIS_rc%udef
     TRMM3B42V6_struc(n)%metdata2 = LIS_rc%udef
     ferror_TRMM3B42V6 = 0
     write(LIS_logunit, *) 'TRMM yr,mo,da,hr for time step start:', &
                           TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart
     call TRMM3B42V6file( name, n, TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart, &
                          TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart, &
                          TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart, &
                          TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart )
     write(LIS_logunit, *)'Getting new TRMM 3B42V6 satellite precip data:', trim(name)
     order = 1
     call read_TRMM3B42V6(n, name, findex, order, ferror_TRMM3B42V6)
   elseif (.NOT. ((TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart_Previous) .AND. &
              (TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart_Previous) .AND. &
              (TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart_Previous) .AND. &
              (TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart_Previous))) then
     write(LIS_logunit, *) 'TRMM yr,mo,da,hr for time step start:', &
                           TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart
     write(LIS_logunit, *)'Values from time step end assigned'
     TRMM3B42V6_struc(n)%metdata1 = TRMM3B42V6_struc(n)%metdata2 ! SY: i.e., equal to TRMM3B42V6_struc(n)%metdata2 before possible change in TRMM3B42V6_struc(n)%metdata2 below
   endif
 
   ! SY: Now update TRMM data corresponding to end of model time step if required.
   if (.NOT. ((TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd_Previous) .AND. &
              (TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd_Previous) .AND. &
              (TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd_Previous) .AND. &
              (TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd .EQ. &
               TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd_Previous))) then
     ferror_TRMM3B42V6 = 0
     write(LIS_logunit, *) 'TRMM yr,mo,da,hr for time step end:', &
                           TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd, &
                           TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd
     call TRMM3B42V6file( name, n, TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd, &
                          TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd, &
                          TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd, &
                          TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd )
     write(LIS_logunit, *)'Getting new TRMM 3B42V6 satellite precip data:', trim(name)
     order = 2
     call read_TRMM3B42V6(n, name, findex, order, ferror_TRMM3B42V6)
   endif
 
   ! SY: Begin reassigning *_Previous values for use in next get* subroutine call
   TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepStart
   TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepStart
   TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepStart
   TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepStart
 
   TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6yr_TStepEnd
   TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6mo_TStepEnd
   TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6da_TStepEnd
   TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd_Previous = &
                              TRMM3B42V6_struc(n)%TRMM3B42V6hr_TStepEnd
   ! SY: End reassigning *_Previous values for use in next get* subroutine call
 
 ! J.Case (4/19/2013) -- test print of the suppdata2 array.
 ! write (96,*) suppdata2(1,:)

!  end if ! SY: if(alarmCheck) then
end subroutine get_TRMM3B42V6

!BOP
! !ROUTINE: TRMM3B42V6file
! \label{TRMM3B42V6file}
!
! !DESCRIPTION: This subroutine puts together 3B42V6 file name for
!               3 hour file intervals
! There are two filename formats: \newline
!  1. original: 3B42.980131.12.6.precipitation, 3B42.980130.3.6.precipitation
!     which was produced after data are converted from hdf to bin format \newline
!  2. renamed: TRMM3B42V6.2005110809, same data, just different file name \newline
!  The raw hdf data files are named like this: 3B42.060105.15.6.HDF
!  It is possible that we may implement HDF file reading in LIS code 
!  to handle raw input. For now, we do the binary format.  -- Yudong 8/25/06
!
! !INTERFACE:
subroutine TRMM3B42V6file( name, n, yr, mo, da, hr)

  use TRMM3B42V6_forcingMod, only : TRMM3B42V6_struc
  use LIS_constantsMod,      only : LIS_CONST_PATH_LEN

!EOP
  implicit none

  integer, intent(in) :: n

!==== Local Variables=======================

  character(len=*) :: name
  character(len=LIS_CONST_PATH_LEN) :: temp, TRMM3B42V6dir
  integer :: yr, mo, da, hr
  integer :: i, j
  integer :: uyr, umo, uda, uhr, umn, uss
  integer :: original

  uyr = mod(yr, 100)          ! 1980 -> 80
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0

   TRMM3B42V6dir = TRMM3B42V6_struc(n)%TRMM3B42V6dir

   name=''        !clean it
   temp=''

   original = 2
   if (original .eq. 1) then     !  1. original: /abc/3B42.980131.12.6.precipitation
     write(temp, '(a, a, I4, I2.2, a, 3I2.2, a, I2, a)') trim(TRMM3B42V6dir), '/', yr, mo, '/3B42.', uyr, umo, uda, '.', &
          uhr,  '.6.precipitation'
   else                          !  2. renamed: TRMM3B42V6.2005110809
     write(temp, '(a, a, I4, I2.2, a, I4, 3I2.2)') trim(TRMM3B42V6dir), '/', yr, mo, '/3B42V6.', yr, umo, uda, uhr
   end if

  !strip off the spaces
  name = trim(temp)

end subroutine TRMM3B42V6file


