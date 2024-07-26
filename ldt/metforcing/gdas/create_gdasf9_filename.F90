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
! !ROUTINE: create_gdasf9_filename
! \label{create_gdasf9_filename}
!
! !INTERFACE:
subroutine create_gdasf9_filename(option, name00, name03, &
                                  gdasdir, yr, mo, da, hr, status )
! !USES: 
  use LDT_timeMgrMod, only : LDT_tick

  implicit none
! !ARGUMENTS: 
  integer,          intent(in)    :: option
  character(len=*), intent(in)    :: gdasdir
  integer, intent(in)             :: yr, mo, da, hr
  character(len=*), intent (out)  :: name00
  character(len=*), intent (out)  :: name03
  integer                         :: status
! !DESCRIPTION:
!   This subroutine puts together 9hr forecast GDAS file name
!
!   First read the notes in
!   create\_gdasfilename.F90 (\ref{create_gdasfilename}).
!
!   This routine is used to create the GDAS file name for the backup 9hr 
!   forecast file.  At analysis hours (hr = 00, 06, 12, and 18), LDT
!   reads the hr.f03 GDAS file for both the new instantaneous and the
!   new averaged forcing values (here, name00 == name03). Should this file be
!   missing, the 9-hour forecast of the previous analysis hour contains
!   comparable forcing values.
!
!   An hr.f09 file contains instantaneous forcing values valid at hour hr+9,
!   and it contains 3-hourly averaged forcing values valid for [hr+6, hr+9]
!
!   For example, when LDT is at hour 06, it reads ahead to get instantaneous
!   forcing values at hour 09 and to get averaged forcing data for the
!   period [06, 09].  The GDAS file 06.f03 contains these data.  Should
!   this file be missing, the GDAS file 00.f09 also contains instantaneous
!   forcing values at hour 09 and averged forcing values for the
!   period [06, 09].
!
!   Note that backup 9-hour forecast files exist only for anaylsis hours.
!   Should a name00 or name03 file be missing when LDT is at a non-analysis
!   hour (hr = 03, 09, 15, or 21), then there is no previous hour with a 9-hour
!   forecast that contains comparable forcing values.  Another strategy must
!   be developed to replace those missing files.
!
!   Note that the logic for reading and creating names for bookend1, which
!   is only done at initialization, is the same as described above except
!   imagine that you are at start hour - 3 reading ahead to the start hour.
!
!   Note that the case where LDT is reading name00, name03, and name06
!   data to read ahead is not handled by this routine.
!   Example, consider the case when reading ahead to hour 18.
!   Here LDT is at hour 15, and LDT needs 18f00 for instantaneous values
!   at hour 18, and LDT needs 12f03 and 12f06 for averaged values for the
!   period [15, 18].  12f03 contains averaged values for [12, 15].  If it
!   is missing, LDT could use 06f09 to obtain averaged values for [12, 15].
!   Again, this scenario is not handled by this routine.  LDT will simply
!   roll back to the the previous day to look for data.
! 
!  The arguments are:
!  \begin{description}
!  \item[option]
!    bookend option (1 or 2)
!   \item[name00]
!   name of the timestamped 9hr GDAS file
!   \item[name03]
!   name of the timestamped 9hr GDAS file
!  \item[gdasdir]
!    Name of the GDAS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!  \item[status]
!   flag indicating if the 9hr forecast is available. 
!  \end{description}
!
!EOP
  integer           :: i, c
  integer           :: uyr, umo, uda, uhr, umn, uss
  integer           :: uyr0, umo0, uda0, uhr0, umn0, uss0
  character(len=2)  :: initcode1,  fcstcode1
  character(len=8)  :: fdir
  character(len=10) :: ftime
  character(len=21) :: fsubs
  real*8            :: dumbtime
  integer           :: doy1, doy
  real              :: gmt
!=== End Variable Definition ===============

!=== formats for filename segments
!-----------------------------------------------------------------
!  Make variables for the time used to create the file
!  We don't want these variables being passed out
!-----------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0

  uyr0 = yr
  umo0 = mo
  uda0 = da
  uhr0 = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn0 = 0
  uss0 = 0

!-----------------------------------------------------------------
! Cheat sheet
!
! bookend 2: 
!  hour 00 look for 18F9
!  hour 06 look for 00F9
!  hour 12 look for 06F9
!  hour 18 look for 12F9
!
! bookend 1: 
!  hour 03 look for 18F9
!  hour 09 look for 00F9
!  hour 15 look for 06F9
!  hour 21 look for 12F9
!-----------------------------------------------------------------

  if ( option == 2 ) then !bookend 2
     ! backup f09 files exist only for analysis hours
     if ( uhr == 0 .or. uhr == 6 .or. uhr == 12 .or. uhr == 18 ) then
        status = 0

        ! read forcing from previous analysis_hour file
        write(unit=initcode1,fmt='(i2.2)') uhr-6

        !special case: need to go back to the previous day 18z
        if ( uhr0 == 0 ) then
           initcode1 = '18'
           call LDT_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, &
                uhr0, umn0, uss0, -1.0*6*60*60)
        endif
     else
        status = 1
     endif
  else !bookend 1
     if (uhr == 3 .or. uhr == 9 .or. uhr == 15 .or. uhr == 21 ) then 
        status = 0

        ! read forcing from previous analysis_hour file
        write(unit=initcode1,fmt='(i2.2)') uhr-9

        !special case: need to go back to the previous day 18z
        if ( uhr0 == 3 ) then
           initcode1 = '18'
           call LDT_tick(dumbtime, doy, gmt, uyr0, umo0, uda0, &
                         uhr0, umn0, uss0, -9.0*60*60)
        endif
     else
        status = 1
     endif
  endif

  if ( status == 0 ) then
     fcstcode1 = '09'
    
     !name00
     write(UNIT=fdir, fmt='(a1, i4.4, i2.2, a1)') '/', uyr0, umo0, '/'
     write(UNIT=ftime, fmt='(i4.4, i2.2, i2.2, a2)') uyr0, umo0, uda0, initcode1
     write(UNIT=fsubs, fmt='(a16, a2, a3)') '.gdas1.sfluxgrbf', fcstcode1, '.sg'
     name00 = trim(gdasdir)//fdir//ftime//fsubs

     !name03
     name03 = name00
  endif

end subroutine create_gdasf9_filename

