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
! !ROUTINE: get_TRMM3B42V7
! \label{get_TRMM3B42V7}
!
! !REVISION HISTORY:
! 21 Jun 2013: Soni Yatheendradas; based on new TRMM 3B42V6 code that has
!              changes from earlier code to avoid (a) alternate file skip,
!              (b) jump to previous day TRMM, and
!              (c) absence of rain rate weighting
!
! !INTERFACE:
subroutine get_TRMM3B42V7(n, findex)
! !USES:
  use LDT_coreMod, only           : LDT_rc, LDT_masterproc
  use LDT_timeMgrMod, only        : LDT_time2date, LDT_tick, LDT_get_nstep, &
                                    LDT_isAlarmRinging ! SY
  use LDT_logMod, only            : LDT_logunit, LDT_endrun
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use TRMM3B42V7_forcingMod, only : TRMM3B42V7_struc
  use LDT_metforcingMod, only     : LDT_forc ! SY

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly TRMM 3B42V7 forcing.
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
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    determines the TRMM 3B42V7 data times
!  \item[TRMM3B42V7file](\ref{TRMM3B42V7file}) \newline
!    Puts together appropriate file filename for 6 hour intervals
!  \item[read\_TRMM3B42V7](\ref{read_TRMM3B42V7}) \newline
!      Interpolates TRMM 3B42V7 data to LDT grid
!  \end{description}
!EOP

!==== Local Variables=======================
  integer :: ferror_TRMM3B42V7      ! Error flag for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1,ts1                
           ! SY: Time parameters for start and end of current LDAS time step
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2               
           ! SY: Time parameters for TRMM data time nearest to start of model time step
  integer :: doy3, yr3, mo3, da3, hr3, mn3, ss3, ts3               
           ! SY: Time parameters for TRMM data time nearest to end of model time step
  real    :: gmt1, gmt2, gmt3     ! SY ,kgmt3, mgmt3
  character(len=LDT_CONST_PATH_LEN) :: filename       ! Filefilename variables for precip data sources
  real*8  :: LDT_timeAtTStepStart_add90min ! SY
  real*8  :: LDT_timeAtTStepEnd_add90min ! SY
  integer :: order
  logical :: alarmCheck ! SY

!=== End Variable Definition =======================

!  Disable alarm.  It throws off the timing of reading the
!  next data set, and it throws off copying the bookends.
!  alarmCheck = LDT_isAlarmRinging(LDT_rc,"TRMM 3B42V7 alarm") ! SY
!  if(alarmCheck) then ! SY

 !------------------------------------------------------------------------
 ! SY : Determine start and end time of model time step
 !------------------------------------------------------------------------
   yr1 = LDT_rc%yr  !current time, ! SY: i.e., end of model time step
   mo1 = LDT_rc%mo
   da1 = LDT_rc%da
   hr1 = LDT_rc%hr
   mn1 = LDT_rc%mn
   ss1 = LDT_rc%ss
   ts1 = 0
   call LDT_tick( TRMM3B42V7_struc(n)%LDT_timeAtTStepEnd, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
   ts1 = (-1)*LDT_rc%nts(n) ! SY: for start of model time step
   call LDT_tick( TRMM3B42V7_struc(n)%LDT_timeAtTStepStart, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
 
 !------------------------------------------------------------------------
 ! SY: 3B42V7 product times.
 !     I get these times based on a very simple principle: modify the
 !     time step start/end by adding 90 min. Whichever TRMM data stamp the
 !     modified start/end coincides with or immediately comes after, will be
 !     the TRMM data nearest to the original time step start/end. Hence,
 !     values of that TRMM data are assigned to that original start/end as
 !     preparation for weighting data in the time_interp* subroutine.
 !------------------------------------------------------------------------
 
   ! SY: add 90 min to start of model time step
   yr2 = LDT_rc%yr
   mo2 = LDT_rc%mo
   da2 = LDT_rc%da
   hr2 = LDT_rc%hr
   mn2 = LDT_rc%mn
   ss2 = LDT_rc%ss
   ts2 = (-1)*LDT_rc%nts(n) + 90*60
   call LDT_tick( LDT_timeAtTStepStart_add90min, doy2, gmt2, &
        yr2, mo2, da2, hr2, mn2, ss2, real(ts2))
   ! SY: Now start calculations for TRMM data time nearest to start of model time step
   TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart = yr2
   TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart = mo2
   TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart = da2
   TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart = 3*(hr2/3)
   mn2 = 0
   ss2 = 0
   ts2 = 0
   call LDT_tick( TRMM3B42V7_struc(n)%TRMM3B42V7time_TStepStart, doy2, gmt2, &
        TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart, &
        TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart, &
        TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart, &
        TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart, mn2, ss2, real(ts2))
 
   ! SY: add 90 min to end of model time step
   yr3 = LDT_rc%yr
   mo3 = LDT_rc%mo
   da3 = LDT_rc%da
   hr3 = LDT_rc%hr
   mn3 = LDT_rc%mn
   ss3 = LDT_rc%ss
   ts3 = 90*60
   call LDT_tick( LDT_timeAtTStepEnd_add90min, doy3, gmt3, &
        yr3, mo3, da3, hr3, mn3, ss3, real(ts3))
   ! SY: Now start calculations for TRMM data time nearest to end of model time step
   TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd = yr3
   TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd = mo3
   TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd = da3
   TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd = 3*(hr3/3)
   mn3 = 0
   ss3 = 0
   ts3 = 0
   call LDT_tick( TRMM3B42V7_struc(n)%TRMM3B42V7time_TStepEnd, doy3, gmt3, &
        TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd, &
        TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd, &
        TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd, &
        TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd, mn3, ss3, real(ts3))
 
 !------------------------------------------------------------------------
 ! Check for and get 3B42V7 observed Precipitation data
 !------------------------------------------------------------------------
 
   ! SY: Now update TRMM data corresponding to start of model time step if required
   if  (LDT_get_nstep(LDT_rc, n).eq. 1 .or. LDT_rc%rstflag(n) .eq. 1) then

      if (LDT_rc%nts(n) .ge. (3*60-1)*60) then ! almost 3 hrs
         write(LDT_logunit,*) 'LDT time step should be < 3 hrs!'
         write(LDT_logunit,*) 'TRMM reader functionality for >= 3 hrs not written yet'
         write(LDT_logunit,*) 'will involve weights combining > 2 TRMM data files'
         write(LDT_logunit,*) 'Program stopping ... '
         call LDT_endrun()
      endif

     LDT_forc(n,findex)%metdata1 = LDT_rc%udef
     LDT_forc(n,findex)%metdata2 = LDT_rc%udef
     ferror_TRMM3B42V7 = 0
     write(LDT_logunit, *) 'TRMM yr,mo,da,hr for time step start:', &
                           TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart
     call TRMM3B42V7file( filename, n, TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart, &
                          TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart, &
                          TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart, &
                          TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart )
     write(LDT_logunit, *)'Getting new TRMM 3B42V7 satellite precip data:', filename
     order = 1
     call read_TRMM3B42V7(n, filename, findex, order, ferror_TRMM3B42V7)
   elseif (.NOT. ((TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart_Previous) .AND. &
              (TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart_Previous) .AND. &
              (TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart_Previous) .AND. &
              (TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart_Previous))) then
     write(LDT_logunit, *) 'TRMM yr,mo,da,hr for time step start:', &
                           TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart
     write(LDT_logunit, *)'Values from time step end assigned'
     LDT_forc(n,findex)%metdata1 = LDT_forc(n,findex)%metdata2 
   ! SY: i.e., equal to LDT_forc(n,findex)%metdata2 before possible change in LDT_forc(n,findex)%metdata2 below
   endif
 
   ! SY: Now update TRMM data corresponding to end of model time step if required.
   if (.NOT. ((TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd_Previous) .AND. &
              (TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd_Previous) .AND. &
              (TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd_Previous) .AND. &
              (TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd .EQ. &
               TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd_Previous))) then
     ferror_TRMM3B42V7 = 0
     write(LDT_logunit, *) 'TRMM yr,mo,da,hr for time step end:', &
                           TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd, &
                           TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd
     call TRMM3B42V7file( filename, n, TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd, &
                          TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd, &
                          TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd, &
                          TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd )
     write(LDT_logunit, *)'Getting new TRMM 3B42V7 satellite precip data:', filename
     order = 2
     call read_TRMM3B42V7(n, filename, findex, order, ferror_TRMM3B42V7)
   endif
 
   ! SY: Begin reassigning *_Previous values for use in next get* subroutine call
   TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart
   TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart
   TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart
   TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart
 
   TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd
   TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd
   TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd
   TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd_Previous = &
                              TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd
   ! SY: End reassigning *_Previous values for use in next get* subroutine call
 
 ! J.Case (4/19/2013) -- test print of the suppdata2 array.
 ! write (96,*) suppdata2(1,:)

end subroutine get_TRMM3B42V7

!BOP
! !ROUTINE: TRMM3B42V7file
! \label{TRMM3B42V7file}
!
! !INTERFACE:
subroutine TRMM3B42V7file( filename, n, yr, mo, da, hr)

! !DESCRIPTION: This subroutine puts together 3B42V7 file filename for
!               3 hour file intervals
! There are two filename formats for the TRMM 3B42V7 data: \newline
!  1. original: 3B42.980131.12.6.precipitation, 3B42.980130.3.6.precipitation
!     which was produced after data are converted from hdf to bin format \newline
!  2. refilenamed: TRMM3B42V6.2005110809, same data, just different file filename \newline
!  The raw hdf data files are filenamed like this: 3B42.060105.15.6.HDF
!  It is possible that we may implement HDF file reading in LDT code 
!  to handle raw input. For now, we do the binary format.  -- Yudong 8/25/06
!

! !USES:
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use TRMM3B42V7_forcingMod, only : TRMM3B42V7_struc

!EOP
  implicit none

  integer, intent(in) :: n

!==== Local Variables=======================

  character(len=LDT_CONST_PATH_LEN) :: filename, TRMM3B42V7dir
  character(len=LDT_CONST_PATH_LEN) :: temp
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

   TRMM3B42V7dir = TRMM3B42V7_struc(n)%TRMM3B42V7dir

   filename=''        !clean it
   temp=''

   original = 2
   if (original .eq. 1) then     !  1. original: /abc/3B42.980131.12.6.precipitation
     write(temp, '(a, a, I4, I2.2, a, 3I2.2, a, I2, a)') &
           trim(TRMM3B42V7dir), '/', yr, mo, '/3B42.', uyr, umo, uda, '.', &
           uhr,  '.6.precipitation'
   else                          !  2. refilenamed: TRMM3B42V7.2005110809
     write(temp, '(a, a, I4, I2.2, a, I4, 3I2.2)') &
           trim(TRMM3B42V7dir), '/', yr, mo, '/3B42V7.', yr, umo, uda, uhr
   end if

  !strip off the spaces
  j = 1
  Do i=1, len(temp)
   if( temp(i:i) .ne. ' ' ) then
      filename(j:j)=temp(i:i)
      j=j+1
   end if
  End Do

end subroutine TRMM3B42V7file


