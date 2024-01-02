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
! !ROUTINE: get_TRMM3B42RTV7
! \label{get_TRMM3B42RTV7}
!
! !REVISION HISTORY:
! 17 Jul 2001: Jon Gottschalck; Initial code
! 10 Oct 2001: Jon Gottschalck; Modified to adjust convective precip
!               using a ratio of the model convective / total ratio
! 06 Jan 2015: KR Arsenault; Added support for latest V7 files
!
! !INTERFACE:
subroutine get_TRMM3B42RTV7(n,findex)

! !USES:
  use LDT_coreMod, only           : LDT_rc, LDT_masterproc
  use LDT_timeMgrMod,  only       : LDT_time2date, LDT_tick, LDT_get_nstep, &
                                    LDT_isAlarmRinging ! SY
  use LDT_logMod, only            : LDT_logunit, LDT_endrun
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use TRMM3B42RTV7_forcingMod, only : TRMM3B42RTV7_struc
  use LDT_metforcingMod,  only    : LDT_forc 

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, TRMM 3B42RT V7 forcing. 
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
!  \item[LDT\_tick](\ref{LDT_tick}) \newline
!    determines the TRMM 3B42RTV7 data times
!  \item[create_3B42RTV7_filename](\ref{TRMM3B42RTV7file}) \newline
!    Puts together appropriate file name for 3 hour intervals
!  \item[read\_TRMM3B42RTV7](\ref{read_TRMM3B42RTV7}) \newline
!    Interpolates TRMM 3B42RTV7 data to LDT grid
!  \end{description}
!EOP

   
!==== Local Variables=======================
  integer :: ferror_TRMM3B42RT      ! Error flag for precip data sources
! SY: Time parameters for start and end of current LDAS time step
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1, ts1     
! SY: Time parameters for TRMM data time nearest to start of model time step
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2, ts2   
! SY: Time parameters for TRMM data time nearest to end of model time step
  integer :: doy3, yr3, mo3, da3, hr3, mn3, ss3, ts3    
! SY: kgmt3, mgmt3
  real    :: gmt1, gmt2, gmt3             

  character(len=LDT_CONST_PATH_LEN) :: filename                ! Filename variables for precip data sources
  real*8  :: LDT_timeAtTStepStart_add90min ! SY
  real*8  :: LDT_timeAtTStepEnd_add90min   ! SY
  logical :: alarmCheck  ! SY
  integer :: order

!=== End Variable Definition =======================

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
  call LDT_tick( TRMM3B42RTV7_struc(n)%LDT_timeAtTStepEnd, doy1, gmt1, &
       yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
  ts1 = (-1)*LDT_rc%nts(n) ! SY: for start of model time step
  call LDT_tick( TRMM3B42RTV7_struc(n)%LDT_timeAtTStepStart, doy1, gmt1, &
       yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
 
 !------------------------------------------------------------------------
 ! SY: 3B42RTV7 product times.
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
   TRMM3B42RTV7_struc(n)%yr_TStepStart = yr2
   TRMM3B42RTV7_struc(n)%mo_TStepStart = mo2
   TRMM3B42RTV7_struc(n)%da_TStepStart = da2
   TRMM3B42RTV7_struc(n)%hr_TStepStart = 3*(hr2/3)
   mn2 = 0
   ss2 = 0
   ts2 = 0
   call LDT_tick( TRMM3B42RTV7_struc(n)%time_TStepStart, doy2, gmt2, &
        TRMM3B42RTV7_struc(n)%yr_TStepStart, &
        TRMM3B42RTV7_struc(n)%mo_TStepStart, &
        TRMM3B42RTV7_struc(n)%da_TStepStart, &
        TRMM3B42RTV7_struc(n)%hr_TStepStart, mn2, ss2, real(ts2))
 
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
   TRMM3B42RTV7_struc(n)%yr_TStepEnd = yr3
   TRMM3B42RTV7_struc(n)%mo_TStepEnd = mo3
   TRMM3B42RTV7_struc(n)%da_TStepEnd = da3
   TRMM3B42RTV7_struc(n)%hr_TStepEnd = 3*(hr3/3)
   mn3 = 0
   ss3 = 0
   ts3 = 0
   call LDT_tick( TRMM3B42RTV7_struc(n)%time_TStepEnd, doy3, gmt3, &
        TRMM3B42RTV7_struc(n)%yr_TStepEnd, &
        TRMM3B42RTV7_struc(n)%mo_TStepEnd, &
        TRMM3B42RTV7_struc(n)%da_TStepEnd, &
        TRMM3B42RTV7_struc(n)%hr_TStepEnd, mn3, ss3, real(ts3))
 
 !------------------------------------------------------------------------
 ! Check for and get 3B42RTV7 observed Precipitation data
 !------------------------------------------------------------------------
 
   ! SY: Now update TRMM data corresponding to start of model time step if required
   if (LDT_get_nstep(LDT_rc, n).eq. 1 .or. LDT_rc%rstflag(n) .eq. 1) then

      if (LDT_rc%nts(n) .ge. (3*60-1)*60) then ! almost 3 hrs
         write(LDT_logunit,*) 'LDT time step should be < 3 hrs!'
         write(LDT_logunit,*) 'TRMM reader functionality for >= 3 hrs not written yet'
         write(LDT_logunit,*) 'will involve weights combining > 2 TRMM data files'
         write(LDT_logunit,*) 'Program stopping ... '
         call LDT_endrun()
      endif

     LDT_forc(n,findex)%metdata1 = LDT_rc%udef
     LDT_forc(n,findex)%metdata2 = LDT_rc%udef
     ferror_TRMM3B42RT = 0
     write(LDT_logunit, *) 'TRMM 3B42RTV7 yr,mo,da,hr for time step start:'
     write(LDT_logunit, *) TRMM3B42RTV7_struc(n)%yr_TStepStart, &
                           TRMM3B42RTV7_struc(n)%mo_TStepStart, &
                           TRMM3B42RTV7_struc(n)%da_TStepStart, &
                           TRMM3B42RTV7_struc(n)%hr_TStepStart

     call create_3B42RTV7_filename( filename, TRMM3B42RTV7_struc(n)%directory, &
                          TRMM3B42RTV7_struc(n)%yr_TStepStart, &
                          TRMM3B42RTV7_struc(n)%mo_TStepStart, &
                          TRMM3B42RTV7_struc(n)%da_TStepStart, &
                          TRMM3B42RTV7_struc(n)%hr_TStepStart )
     write(LDT_logunit, *)'Getting new TRMM 3B42RTV7 satellite precip data(1):'!,&
!           trim(filename)

     order = 1
     call read_TRMM3B42RTV7(n, filename, findex, order, ferror_TRMM3B42RT)

   elseif (.NOT. ((TRMM3B42RTV7_struc(n)%yr_TStepStart .EQ. &
               TRMM3B42RTV7_struc(n)%yr_TStepStart_Previous) .AND. &
              (TRMM3B42RTV7_struc(n)%mo_TStepStart .EQ. &
               TRMM3B42RTV7_struc(n)%mo_TStepStart_Previous) .AND. &
              (TRMM3B42RTV7_struc(n)%da_TStepStart .EQ. &
               TRMM3B42RTV7_struc(n)%da_TStepStart_Previous) .AND. &
              (TRMM3B42RTV7_struc(n)%hr_TStepStart .EQ. &
               TRMM3B42RTV7_struc(n)%hr_TStepStart_Previous))) then
     write(LDT_logunit,*) 'TRMM 3B42RTV7 yr,mo,da,hr for time step start:'
     write(LDT_logunit,*)  TRMM3B42RTV7_struc(n)%yr_TStepStart, &
                           TRMM3B42RTV7_struc(n)%mo_TStepStart, &
                           TRMM3B42RTV7_struc(n)%da_TStepStart, &
                           TRMM3B42RTV7_struc(n)%hr_TStepStart
     write(LDT_logunit, *)'Values from time step end assigned'
     LDT_forc(n,findex)%metdata1 = LDT_forc(n,findex)%metdata2 
     ! SY: i.e., equal to LDT_forc(n,findex)%metdata2 before possible 
     !     change in LDT_forc(n,findex)%metdata2 below
   endif

   ! SY: Now update TRMM data corresponding to end of model time step if required.
   if (.NOT. ((TRMM3B42RTV7_struc(n)%yr_TStepEnd .EQ. &
               TRMM3B42RTV7_struc(n)%yr_TStepEnd_Previous) .AND. &
              (TRMM3B42RTV7_struc(n)%mo_TStepEnd .EQ. &
               TRMM3B42RTV7_struc(n)%mo_TStepEnd_Previous) .AND. &
              (TRMM3B42RTV7_struc(n)%da_TStepEnd .EQ. &
               TRMM3B42RTV7_struc(n)%da_TStepEnd_Previous) .AND. &
              (TRMM3B42RTV7_struc(n)%hr_TStepEnd .EQ. &
               TRMM3B42RTV7_struc(n)%hr_TStepEnd_Previous))) then
     ferror_TRMM3B42RT = 0
     write(LDT_logunit, *) 'TRMM 3B42RTV7 yr,mo,da,hr for time step end:'
     write(LDT_logunit, *) TRMM3B42RTV7_struc(n)%yr_TStepEnd, &
                           TRMM3B42RTV7_struc(n)%mo_TStepEnd, &
                           TRMM3B42RTV7_struc(n)%da_TStepEnd, &
                           TRMM3B42RTV7_struc(n)%hr_TStepEnd
     call create_3B42RTV7_filename( filename, TRMM3B42RTV7_struc(n)%directory,&
                          TRMM3B42RTV7_struc(n)%yr_TStepEnd, &
                          TRMM3B42RTV7_struc(n)%mo_TStepEnd, &
                          TRMM3B42RTV7_struc(n)%da_TStepEnd, &
                          TRMM3B42RTV7_struc(n)%hr_TStepEnd )
     write(LDT_logunit, *)'Getting new TRMM 3B42RTV7 satellite precip data(2):'!,&
!           trim(filename)
     order = 2
     call read_TRMM3B42RTV7(n, filename, findex, order, ferror_TRMM3B42RT)

   endif

   ! SY: Begin reassigning *_Previous values for use in next get* subroutine call
   TRMM3B42RTV7_struc(n)%yr_TStepStart_Previous = &
                              TRMM3B42RTV7_struc(n)%yr_TStepStart
   TRMM3B42RTV7_struc(n)%mo_TStepStart_Previous = &
                              TRMM3B42RTV7_struc(n)%mo_TStepStart
   TRMM3B42RTV7_struc(n)%da_TStepStart_Previous = &
                              TRMM3B42RTV7_struc(n)%da_TStepStart
   TRMM3B42RTV7_struc(n)%hr_TStepStart_Previous = &
                              TRMM3B42RTV7_struc(n)%hr_TStepStart
 
   TRMM3B42RTV7_struc(n)%yr_TStepEnd_Previous = &
                              TRMM3B42RTV7_struc(n)%yr_TStepEnd
   TRMM3B42RTV7_struc(n)%mo_TStepEnd_Previous = &
                              TRMM3B42RTV7_struc(n)%mo_TStepEnd
   TRMM3B42RTV7_struc(n)%da_TStepEnd_Previous = &
                              TRMM3B42RTV7_struc(n)%da_TStepEnd
   TRMM3B42RTV7_struc(n)%hr_TStepEnd_Previous = &
                              TRMM3B42RTV7_struc(n)%hr_TStepEnd
   ! SY: End reassigning *_Previous values for use in next get* subroutine call
 
end subroutine get_TRMM3B42RTV7

!BOP
! !ROUTINE: create_3B42RTV7_filename
! \label{create_3B42RTV7_filename}
!
!
! !INTERFACE:
subroutine create_3B42RTV7_filename( filename, TRMM3B42RTdir, yr, mo, da, hr)

! !USES: 
  use LDT_timeMgrMod, only : LDT_date2time
  use TRMM3B42RTV7_forcingMod, only : TRMM3B42RTV7_struc

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: TRMM3B42RTdir
  integer           :: yr, mo, da, hr
!
! !DESCRIPTION:
!   This subroutine puts together TRMM 3B42RTV7 file name for 
!   3 hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[TRMM3B42RTdir]
!    Name of the TRMM 3B42RTV7 data directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   filename of the timestamped TRMM 3B42RTV7 file
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

  call LDT_date2time(utime, udoy, ugmt, uyr, umo, uda, uhr, umn, uss ) 

  if (utime .LT. TRMM3B42RTV7_struc(1)%griduptime1) then   
    filename = trim(TRMM3B42RTdir)//"/"//y4//m2//"/3B42RT."//y4//m2//d2//h2//".7R2"
  else
    filename = trim(TRMM3B42RTdir)//"/"//y4//m2//"/3B42RT."//y4//m2//d2//h2//".7"
  end if 

end subroutine create_3B42RTV7_filename

