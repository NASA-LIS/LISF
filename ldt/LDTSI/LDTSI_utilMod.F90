!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LDTSI_utilMod
! 
! REVISION HISTORY:
! 08 Feb 2019  Eric Kemp  First ported to LDT.
! 09 May 2019  Eric Kemp  Renamed LDTSI.
!
! DESCRIPTION:
! Source code for util library for Air Force snow depth analysis.
!-------------------------------------------------------------------------

#include "LDT_misc.h"
   
module LDTSI_utilMod

   !Defaults
   implicit none
   private

   ! Public methods
   public :: abort_message
   public :: date10_julhr
   public :: error_message
   public :: julhr_date10
   public :: putget_int1
   public :: putget_int
   public :: putget_real
   public :: tmjul4

contains

   subroutine abort_message (program_name, routine_name, message)

      !*************************************************************************
      !*************************************************************************
      !**
      !**  NAME: ABORT MESSAGE
      !**
      !**  PURPOSE:  GENERATE STANDARD ABORT MESSAGE AND ABORT PROGRAM
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  ROUTINE_NAME, MESSAGE
      !**
      !**   OUTPUT: NONE
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILENAME/UNIT #                     W  DESCRIPTION
      !**  ---------------------------------- --- ----------------------------
      !**  abort.${ROUTINE_NAME}.${NN} / 99    R  ABORT MESSAGE
      !**
      !**  UPDATES
      !**  =======
      !**  11 MAY 11 INITIAL VERSION (FROM ERROR_MESSAGE)...MR LEWISTON/16WS/WXE
      !**  22 Mar 19 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**
      !*************************************************************************
      !*************************************************************************
      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun

      ! Defaults
      implicit none

      ! Local constants
      integer,      parameter     :: msglns = 20      ! MAXIMUM NUMBER OF LINES IN MESSAGE

      ! Arguments
      character*12, intent(in)    :: program_name     ! NAME OF CALLING ROUTINE
      character*12, intent(in)    :: routine_name     ! NAME OF CALLING ROUTINE
      character*90, intent(in)    :: message (msglns) ! ERROR MESSAGE FROM CALLER

      ! Local variables
      character*7                 :: access_type      ! FILE ACCESS TYPE
      character*100               :: errmsg  (msglns) ! ERROR MESSAGE TO OUTPUT
      character*40                :: message_file     ! MESSAGE FILE NAME
      integer                     :: i                ! DO LOOP COUNTER
      integer                     :: istat            ! I/O STATUS
      integer                     :: nlines           ! NUMBER OF LINES IN MESSAGE

      ! FIND OUT HOW MANY NON-BLANK LINES THERE ARE.
      errmsg = ' '
      do i = msglns, 1 , -1
         if (len_trim(message(i)) .gt. 0) then
            nlines = i
            exit
         endif
      enddo
      nlines = nlines + 2

      ! SET FIRST TWO LINES, THEN COPY CALLER'S LINES TO OUTPUT MESSAGE.
      errmsg(1) = 'PROGRAM: ' // trim (program_name) // ' - FATAL ERROR'
      errmsg(2) = 'ROUTINE: ' // trim (routine_name)
      do i = 3, nlines
         errmsg(i) = trim(message(i-2))
      enddo

      ! SET OUTPUT MESSAGE FILE NAME.
      message_file = 'abort.' // trim(program_name)

      ! OPEN MESSAGE FILE.
      ! WRITE MESSAGE TO LOG AND MESSAGE FILE.
      access_type = 'OPENING'

      open (unit=99, file=trim(message_file), err=5000, iostat=istat)

      access_type = 'WRITING'

      write (ldt_logunit, 6000)
      do i = 1, nlines
         write (ldt_logunit,  6200) trim(errmsg(i))
         write (99, 6400, err=5000) trim(errmsg(i))
      enddo
      write (ldt_logunit, 6600)

      close (99)
      !call abort
      call LDT_endrun()
      
      ! ERROR-HANDLING SECTION.
5000  write (ldt_logunit, 8000) routine_name, access_type, istat
      !call abort
      call LDT_endrun()

      ! FORMAT STATEMENTS.
6000  format (/, 1X, 75('*'))
6200  format ('[ERR] ', A)
6400  format (A)
6600  format (1X, 75('*'))
8000  format (/, 1X, 75('*'),                                          &
           /, 1X, '[ERR]  ABNORMAL ABORT FROM ABORT_MESSAGE ROUTINE',  &
           /, 1X, '[ERR]  CALLED BY ', A12,                            &
           /, 1X, '[ERR]  ERROR WHILE ', A7, ' MESSAGE FILE',          &
           /, 1X, '[ERR]  ISTAT = ', I6,                               &
           /, 1X, 75('*'))

   end subroutine abort_message

   subroutine date10_julhr(date10, julhr, program_name, routine_name)

      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      !
      !     Name: Date/time to Julian hour
      !     ====
      !
      !     Purpose:
      !     =======
      !     Convert a ten-digit date/time group (yyyymmddhh) to Julian hour.
      !
      !     Method:
      !     ======
      !     - Convert character date/time group to integer year, month,
      !       day and hour.
      !     - Check for invalid date/times. send back bad status -
      !       julhr of -1 if invalid.
      !     - If date/times good, calculate julian hour by summing the
      !       the number of hours since 1/1/68, accounting for leap years.
      !
      !     Interface:
      !     =========
      !     Inputs:   date10
      !
      !     Outputs:  julhr
      !
      !     Variables:
      !     =========
      !     Label           Description
      !     ---------       --------------------------------------------------
      !     date10          Ten-digit date/time group yyyymmddhh
      !     days            Total days
      !     dd              Day of the month
      !     flagly          Remainder flag to indicate current year is
      !                     a leap year
      !     hh              Time of day in hours and minutes
      !     hours           Total hours
      !     julhr           Total julian hours
      !     lyears          Number of leap years
      !     message         Error message
      !     mm              Month of the year
      !     month           Table of days accumulated at start of month
      !     temphr          Total hours in time
      !     total           Hours in time rounded to nearest hour
      !     years           Total years
      !     yyyy            Year
      !
      !     Updates:
      !     =======
      !     10 Aug 1999  Ported to IBM SP2.......................Mr Gayno/DNXM
      !     24 Feb 2010  Converted to FTN 90 for SNODEP...Mr Lewiston/16WS/WXE
      !     22 Mar 2019  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !     09 May 2019  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------

      ! Defaults
      implicit none

      ! Arguments
      character*10,  intent(in)   :: date10           ! 10-DIGIT DATE-TIME GROUP
      integer,       intent(out)  :: julhr            ! TOTAL JULIAN HOURS
      character*12,  intent(in)   :: program_name     ! NAME OF CALLING PROGRAM
      character*12,  intent(in)   :: routine_name     ! NAME OF CALLING ROUTINE

      ! Local variables
      character*90                :: message     (20) ! ERROR MESSAGE ! EMK
      character*12                :: librtne_name     ! NAME OF THIS ROUTINE
      integer                     :: days             ! TOTAL DAYS
      integer                     :: dd               ! DAY OF THE MONTH
      integer                     :: flagly           ! LEAP YEAR FLAG
      integer                     :: hh               ! HOUR OF THE DAY
      integer                     :: hours            ! TOTAL HOURS
      integer                     :: lyears           ! NUMBER OF LEAP YEARS
      integer                     :: mm               ! MONTH
      integer                     :: month       (12) ! TABLE OF DAYS ACCUMULATED AT START OF MONTH
      integer                     :: temphr           ! TOTAL HOURS IN TIME
      integer                     :: thours           ! HOURS IN TIME ROUNDED TO NEAREST HOUR
      integer                     :: years            ! TOTAL YEARS
      integer                     :: yyyy             ! YEAR

      ! Data initialization section.
      data month(1)  /0/
      data month(2)  /31/
      data month(3)  /59/
      data month(4)  /90/
      data month(5)  /120/
      data month(6)  /151/
      data month(7)  /181/
      data month(8)  /212/
      data month(9)  /243/
      data month(10) /273/
      data month(11) /304/
      data month(12) /334/
      data librtne_name / 'DATE10_JULHR' /

      ! Convert the 10 character date string into integer values.
      read(date10(1:4), '(i4)', err=9000) yyyy
      read(date10(5:6), '(i2)', err=9000) mm
      read(date10(7:8), '(i2)', err=9000) dd
      read(date10(9:10),'(i2)', err=9000) hh

      hh = hh * 100

      ! First check for valid hour, day, month and year.
      julhr = 0

      if ((hh .lt. 0 .or. hh .gt. 2400) .or.                             &
           (dd .lt. 1 .or. dd .gt. 31)   .or.                            &
           (mm .lt. 1 .or. mm .gt. 12)   .or.                            &
           (yyyy .lt. 1968 .or. yyyy .gt. 9999))  then

         julhr  = -1

      endif

      ! Flag invalid day/month combinations.
      ! If february, check for correct days in a leap year.
      if (mm .eq. 2) then
         if ((mod((yyyy - 68), 4) .ne. 0) .and. (dd .gt. 28)) then
            julhr  = -1
         endif
      elseif ((dd .gt. 30) .and.                                        &
           ((mm .eq. 4) .or.                                            &
           (mm .eq. 6) .or.                                             &
           (mm .eq. 9) .or.                                             &
           (mm .eq. 11))) then
         julhr  = -1
      endif

      ! If an input error was detected, return the -1 status,
      ! otherwise calculate the julian hour.
      if (julhr .ne. -1) then

         ! Calculate number of 'complete' years past 1968.
         years = yyyy - 1968

         ! Flag if this year is a leap year.
         flagly = mod((yyyy - 1968), 4)

         ! Calculate number of days (only 'complete' years) since 1968.
         days = years * 365

         ! Calculate number of leap years
         ! (not including 1968, if 68 is the passed yyyy)
         ! add 3 for rounding, since the operator '/' truncates
         years  = years + 3
         lyears = years / 4

         ! Add one day for every 'complete' leap year to total days.
         days = days + lyears

         ! Add in this year's number of days.
         days = days + month(mm) + dd

         ! Calculate the total number of hours
         hours = days * 24

         ! Calculate today's number of hours. rounded to nearest hour.
         ! note :
         !   since the operator '/' truncates for integer operations, we
         !   add 70 minutes to the passed hh, to facilitate being able to
         !   properly round a time.  examples:
         !           (1430 hrs  will 'round' properly to 1500 hrs if we add 70)
         !             1430 + 70 = 1500 / 100 = 15
         !           (1400 hrs will still 'round' properly to 1400 hrs)
         !             1400 + 70 = 1470 / 100 = 14
         temphr = hh + 70
         thours = temphr / 100

         ! Add today's number of hours to the total number of hours.
         thours = thours + hours

         ! If this a leap year and past feb, add 24 hours.
         if ((flagly .eq. 0) .and. (mm .gt. 2)) then
            thours = thours + 24
         endif
         julhr = thours
      endif

      return

      ! Error handling section.
9000  continue

      message(1) = '[ERR] LIBROUTINE: ' // librtne_name
      message(2) = '[ERR]  ERROR CONVERTING CHARACTER DATE/TIMES TO INTEGER'
      call abort_message (program_name, routine_name, message)
      return

   end subroutine date10_julhr

   subroutine error_message (program_name, routine_name, message)

      !*************************************************************************
      !*************************************************************************
      !**
      !**  NAME: ERROR MESSAGE
      !**
      !**  PURPOSE:  GENERATE STANDARD ERROR MESSAGE
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  ROUTINE_NAME, MESSAGE
      !**
      !**   OUTPUT: NONE
      !**
      !**  FILES ACCESSED
      !**  ==============
      !**  FILENAME/UNIT #                     W  DESCRIPTION
      !**  ---------------------------------- --- ----------------------------
      !**  alert.${ROUTINE_NAME}.${NN} / 99    R  ERROR MESSAGE
      !**
      !**  UPDATES
      !**  =======
      !**  22 FEB 01 INITIAL VERSION (REUSED FROM ADVCLD).......SSGT MILLER/DNXM
      !**  21 JUL 04 CONVERTED TO FORTRAN 90 FOR 16TH MESH......MR EYLANDER/DNXM
      !**  10 MAY 11 UPDATED FOR SNODEP UPGRADE.............MR LEWISTON/16WS/WXE
      !**  22 Mar 19  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 Ma  19  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*************************************************************************
      !*************************************************************************
      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun

      ! Defaults
      implicit none

      ! Local constants
      integer,      parameter     :: msglns = 20      ! MAXIMUM NUMBER OF LINES IN MESSAGE

      ! Arguments
      character*12, intent(in)    :: program_name     ! NAME OF CALLING ROUTINE
      character*12, intent(in)    :: routine_name     ! NAME OF CALLING ROUTINE
      character*90, intent(in)    :: message (msglns) ! ERROR MESSAGE FROM CALLER
      ! Local variables
      character*7                 :: access_type      ! FILE ACCESS TYPE
      character*2                 :: calert_number    ! ALERT NUMBER FOR FILE NAME
      character*100               :: errmsg  (msglns) ! ERROR MESSAGE TO OUTPUT
      character*40                :: message_file     ! MESSAGE FILE NAME
      integer                     :: alert_number     ! ALERT NUMBER
      integer                     :: i                ! DO LOOP COUNTER
      integer                     :: istat            ! I/O STATUS
      integer                     :: nlines           ! NUMBER OF LINES IN MESSAGE
      logical                     :: isfile           ! FLAG INDICATING WHETHER FILE EXISTS

      ! FIND OUT HOW MANY NON-BLANK LINES CALLER PROVIDED.
      ! INCREASE SIZE TO ALLOW NEW FIRST LINE.
      errmsg = ' '
      do i = msglns, 1 , -1
         if (len_trim(message(i)) .gt. 0) then
            nlines = i
            exit
         endif
      enddo
      nlines = nlines + 2

      ! SET FIRST TWO LINES, THEN COPY CALLER'S LINES TO OUTPUT MESSAGE.
      errmsg(1) = '[WARN] PROGRAM: ' // trim (program_name)
      errmsg(2) = '[WARN] ROUTINE: ' // trim (routine_name)
      do i = 3, nlines
         errmsg(i) = trim(message(i-2))
      enddo

      ! SET OUTPUT MESSAGE FILE NAME.
      ! USE ALERT_NUMBER TO AVOID OVERWRITING EXISTING FILES.
      alert_number = 1
      message_file = ' '
      isfile       = .true.

      do while (isfile)
         write (calert_number, '(i2.2)', iostat = istat) alert_number
         message_file = 'alert.' // trim(routine_name) // '.' //         &
              calert_number
         inquire (file=message_file, exist=isfile)
         alert_number = alert_number + 1
      end do

      ! OPEN MESSAGE FILE.
      ! WRITE MESSAGE TO LOG AND MESSAGE FILE.
      access_type = 'OPENING'
      open (unit=99, file=trim(message_file), err=5000, iostat=istat)
      access_type = 'WRITING'
      write (ldt_logunit, 6000)
      do i = 1, nlines
         write (ldt_logunit,  6200) trim(errmsg(i))
         write (99, 6400, err=5000) trim(errmsg(i))
      enddo
      write (ldt_logunit, 6600)

      close (99)

      return

      ! ERROR-HANDLING SECTION.
5000  write (ldt_logunit, 8000) routine_name, access_type, istat
      !call abort
      call LDT_endrun()

      ! FORMAT STATEMENTS.
6000  format (/, 1X, 75('*'))
6200  format ('[WARN] ', A)
6400  format (A)
6600  format (1X, 75('*'))
8000  format (/, 1X, 75('*'),                                          &
           /, 1X, '[ERR]  ABNORMAL ABORT FROM ERROR_MESSAGE ROUTINE',  &
           /, 1X, '[ERR]  CALLED BY ', A50,                            &
           /, 1X, '[ERR]  ERROR WHILE ', A7, ' MESSAGE FILE',          &
           /, 1X, '[ERR]  ISTAT = ', I6,                               &
           /, 1X, 70('*'))

   end subroutine error_message

   subroutine julhr_date10( julhr, date10, program_name, routine_call )

      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      !
      !     Name: Julhr to 10 character date/time group
      !     ====
      !
      !     Purpose:
      !     =======
      !     Convert from a Julian hour to a ten-digit date/time group
      !     (yyyymmddhh)
      !
      !     Method:
      !     ======
      !     - Call utility routine tmjul4 to convert from julian hours
      !       to year, month, day, hour.
      !     - Perform several checks to ensure date information passed
      !       back from tmjul4 is valid.
      !     - If date is good, convert from integer data to a ten-digit
      !       date/time group and pass back to calling routine.
      !
      !     Interface:
      !     =========
      !     Inputs:   JULHR
      !
      !     Outputs:  DATE10
      !
      !     Updates:
      !     =======
      !     15 Oct 1998  Initial version.........................Mr Moore/DNXM
      !     10 Aug 1999  Ported to IBM SP2.  Added intent attributes to
      !                  arguments...............................Mr Gayno/DNXM
      !     24 Feb 2010  Converted to FTN 90 for SNODEP...Mr Lewiston/16WS/WXE
      !     22 Mar 2019  Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !     09 May 2019  Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------

      ! Defaults
      implicit none

      ! Arguments
      integer,        intent(in)  :: julhr            ! AFWA JULIAN HOUR
      character*10,  intent(out)  :: date10           ! 10-DIGIT DATE-TIME GROUP
      character*12,  intent(in)   :: program_name     ! NAME OF CALLING PROGRAM
      character*12,  intent(in)   :: routine_call     ! NAME OF CALLING ROUTINE

      ! Local variables
      character*90                :: message     (20) ! ERROR MESSAGE ! EMK
      integer                     :: dd               ! DAY OF THE MONTH
      integer                     :: hh               ! HOUR OF THE DAY
      integer                     :: j                ! DO LOOP COUNTER
      integer                     :: mm               ! MONTH
      integer                     :: yyyy             ! YEAR

      ! Executable code begins here.
      ! Use tmjul4 to convert julhr to hour, day, month and year.
      message = ' '

      call tmjul4( hh, dd, mm, yyyy, julhr )
     
      ! Check for valid hour, day, month, and year
      if( (  hh .lt.    0 .or.   hh .gt.   23) .or.                     &
           (  dd .lt.    1 .or.   dd .gt.   31) .or.                    &
           (  mm .lt.    1 .or.   mm .gt.   12) .or.                    &
           (yyyy .lt. 1968 .or. yyyy .gt. 9999) ) then
         goto 9000
      endif

      ! Flag invalid day/month combinations, for example...
      ! if February, check for correct number of days in month.
      if( mm .eq. 2 )then

         if( (mod((yyyy - 68), 4) .ne. 0) .and. (dd .gt. 28) )then

            goto 9100

         endif

      elseif (  (dd .gt. 30) .and.                                     &
           ( (mm .eq. 4)  .or.                                         &
           (mm .eq. 6)  .or.                                           &
           (mm .eq. 9)  .or.                                           &
           (mm .eq. 11) ) ) then

         goto 9200

      endif

      ! Convert integer values into a 10 character date/time string
      write(date10(1:4) ,'(i4)', err=9300) yyyy
      write(date10(5:6) ,'(i2)', err=9300) mm
      write(date10(7:8) ,'(i2)', err=9300) dd
      write(date10(9:10),'(i2)', err=9300) hh

      do j = 1, 10
         if( date10(j:j) .eq. ' ' )then
            date10(j:j) = '0'
         endif
      enddo

      return

      ! Error handling section.
9000  continue
      message(1) = '[ERR] LIBROUTINE: JULHR_DATE10'
      message(2) = '[ERR] INVALID JULHR TO DATE/TIME CONVERSION'
      call abort_message (program_name, routine_call, message)
      return

9100  continue
      message(1) = '[ERR] LIBROUTINE: JULHR_DATE10'
      message(2) = '[ERR] CREATED 29 FEB IN NON-LEAP YEAR'
      call abort_message (program_name, routine_call, message)
      return

9200  continue
      message(1) = '[ERR] LIBROUTINE: JULHR_DATE10'
      message(2) = '[ERR] CREATED 31ST DAY FOR MONTH WITH ONLY 30 DAYS'
      call abort_message (program_name, routine_call, message)
      return

9300  continue
      message(1) = '[ERR] LIBROUTINE: JULHR_DATE10'
      message(2) = '[ERR] ERROR CONVERTING INTEGER DATE/TIMES TO CHARACTER'
      call abort_message (program_name, routine_call, message)
      return

   end subroutine julhr_date10

   subroutine putget_int1 ( buffer, iofunc, file_name, program_name, &
        routine_name, igrid, jgrid )

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: Put or get 1-byte integer array
      !**
      !**  PURPOSE:  Read 1-byte integer data from or write to direct access files.
      !**
      !**  METHOD
      !**  ======
      !**   - Open file
      !**   - If iofunc = 'r' or 'R', read buffer from file, abort on error
      !**   - If iofunc = 'w' or 'W', write buffer to file, abort on error
      !**   - If iofunc = anything else, abort
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  BUFFER, IOFUNC, FILE_NAME, PROGRAM_NAME, ROUTINE_NAME, IGRID, JGRID
      !**
      !**   OUTPUT: BUFFER
      !**
      !**  FILES ACCESSED:  Defined in FILE_NAME
      !**
      !**  UPDATES
      !**  =======
      !**  24 Feb 10 Initial version (adapted from putget_int)....Mr Lewiston/16WS/WXE
      !**  22 Mar 19 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Defaults
      implicit none

      ! Local constants
      integer,      parameter      :: msglns = 20                ! MAX NUMBER OF LINES IN MESSAGE

      ! Arguments
      integer*1,     intent(inout) :: buffer      (igrid, jgrid) ! DATA TO BE READ/WRITTEN
      character*1,   intent(in)    :: iofunc                     ! I/O FUNCTION ('r', 'w')
      character*100, intent(in)    :: file_name                  ! FILE PATH AND NAME
      character*12,  intent(in)    :: program_name               ! NAME OF CALLING PROGRAM
      character*12,  intent(in)    :: routine_name               ! NAME OF CALLING ROUTINE
      integer,       intent(in)    :: igrid                      ! SIZE OF GRID IN I-DIRECTION
      integer,       intent(in)    :: jgrid                      ! SIZE OF GRID IN I-DIRECTION

      ! Local variables
      character*7                  :: access_type                ! FILE ACCESS TYPE
      character*4                  :: cstat                      ! I/O STATUS FOR MESSAGE
      character*90                 :: message     (msglns)       ! ERROR MESSAGE
      integer                      :: istat                      ! I/O STATUS
      integer                      :: istat1                     ! CONVERSION STATUS
      integer                      :: reclen                     ! FILE RECORD LENGTH
      logical                      :: isopen                     ! INDICATES IF FILE IS OPEN

      ! Executable code starts here. Open file, abort on error.
      isopen  = .false.
      message = ' '
      reclen  = igrid * jgrid

      access_type = 'OPENING'
      open( 2, file=trim(file_name), form='unformatted',                &
           access='direct', recl=reclen, iostat=istat, err=1000 )
      isopen = .true.

      ! Read from file, abort on error
      if ( (iofunc .eq. 'r') .or. (iofunc .eq. 'R') ) then
         access_type = 'READING'
         read( 2, rec=1, iostat=istat, err=1000 ) buffer

      ! write to file, abort on error
      else if ( (iofunc .eq. 'w') .or. (iofunc .eq. 'W') ) then
         access_type = 'WRITING'
         write( 2, rec=1, iostat=istat, err=1000 ) buffer

      ! Else abort due to invalid IOFUNC value.
      else
         go to 4000
      endif

      close(2)

      return

      ! Error handling.

1000  continue
      if (isopen) close (2)
      message(1) = '[ERR] LIBROUTINE: PUTGET_INT1'
      message(2) = '[ERR] ERROR ' // access_type // ' ' // trim (file_name)
      write (cstat, '(i4)', iostat=istat1) istat
      if (istat1 .eq. 0) message(3) = '[ERR] ISTAT = ' // cstat
      call abort_message (program_name, routine_name, message)
      return

4000  continue
      close(2)
      message(1) = '[ERR] LIBROUTINE: PUTGET_INT1'
      message(2) = '[ERR] INVALID IOFUNC VALUE = '// trim(iofunc)
      call abort_message (program_name, routine_name, message)
      return

   end subroutine putget_int1

   subroutine putget_int ( buffer, iofunc, file_name, program_name,  &
        routine_name, igrid, jgrid )

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: Put or get integer array
      !**
      !**  PURPOSE:  Read integer data from or write to direct access files.
      !**
      !**  METHOD
      !**  ======
      !**   - Open file
      !**   - If iofunc = 'r' or 'R', read buffer from file, abort on error
      !**   - If iofunc = 'w' or 'W', write buffer to file, abort on error
      !**   - If iofunc = anything else, abort
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  BUFFER, IOFUNC, FILE_NAME, PROGRAM_NAME, ROUTINE_NAME, IGRID, JGRID
      !**
      !**   OUTPUT: BUFFER
      !**
      !**  FILES ACCESSED:  Defined in FILE_NAME
      !**
      !**  UPDATES
      !**  =======
      !**  04 Aug 99 Initial AGRMET version..............................Mr Gayno/DNXM
      !**  24 Feb 10 Converted to FORTRAN 90 for SNODEP...........Mr Lewiston/16WS/WXE
      !**  22 Mar 19 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19 Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Defaults
      implicit none

      ! Local constants
      integer,      parameter      :: msglns = 20                ! MAX NUMBER OF LINES IN MESSAGE

      ! Arguments
      integer,       intent(inout) :: buffer      (igrid, jgrid) ! DATA TO BE READ/WRITTEN
      character*1,   intent(in)    :: iofunc                     ! I/O FUNCTION ('r', 'w')
      character*100, intent(in)    :: file_name                  ! FILE PATH AND NAME
      character*12,  intent(in)    :: program_name               ! NAME OF CALLING PROGRAM
      character*12,  intent(in)    :: routine_name               ! NAME OF CALLING ROUTINE
      integer,       intent(in)    :: igrid                      ! SIZE OF GRID IN I-DIRECTION
      integer,       intent(in)    :: jgrid                      ! SIZE OF GRID IN I-DIRECTION

      ! Local variables
      character*7                  :: access_type                ! FILE ACCESS TYPE
      character*4                  :: cstat                      ! I/O STATUS FOR MESSAGE
      character*90                 :: message     (msglns)       ! ERROR MESSAGE
      integer                      :: istat                      ! I/O STATUS
      integer                      :: istat1                     ! CONVERSION STATUS
      integer                      :: reclen                     ! FILE RECORD LENGTH
      logical                      :: isopen                     ! INDICATES IF FILE IS OPEN

      ! Executable code starts here. Open file, abort on error.
      isopen  = .false.
      message = ' '
      reclen  = igrid * jgrid * 4

      access_type = 'OPENING'
      open( 2, file=trim(file_name), form='unformatted',                &
           access='direct', recl=reclen, iostat=istat, err=1000 )
      isopen = .true.

      ! Read from file, abort on error
      if ( (iofunc .eq. 'r') .or. (iofunc .eq. 'R') ) then
         access_type = 'READING'
         read( 2, rec=1, iostat=istat, err=1000 ) buffer

      ! Write to file, abort on error
      else if ( (iofunc .eq. 'w') .or. (iofunc .eq. 'W') ) then
         access_type = 'WRITING'
         write( 2, rec=1, iostat=istat, err=1000 ) buffer

      ! Else abort due to invalid IOFUNC value.         
      else
         go to 4000
      endif

      close(2)

      return

      ! Error handling.
1000  continue
      if (isopen) close (2)
      message(1) = '[ERR] LIBROUTINE: PUTGET_INT'
      message(2) = '[ERR] ERROR ' // access_type // ' ' // trim (file_name)
      write (cstat, '(i4)', iostat=istat1) istat
      if (istat1 .eq. 0) message(3) = '[ERR] ISTAT = ' // cstat
      call abort_message (program_name, routine_name, message)
      return

4000  continue
      close(2)
      message(1) = '[ERR] LIBROUTINE: PUTGET_INT'
      message(2) = '[ERR] INVALID IOFUNC VALUE = '// trim(iofunc)
      call abort_message (program_name, routine_name, message)
      return

   end subroutine putget_int

   subroutine putget_real ( buffer, iofunc, file_name, program_name, &
        routine_name, igrid, jgrid )

      !*******************************************************************************
      !*******************************************************************************
      !**
      !**  NAME: Put or get real array
      !**
      !**  PURPOSE:  Read real data from or write to direct access files.
      !**
      !**  METHOD
      !**  ======
      !**   - Open file
      !**   - If iofunc = 'r' or 'R', read buffer from file, abort on error
      !**   - If iofunc = 'w' or 'W', write buffer to file, abort on error
      !**   - If iofunc = anything else, abort
      !**
      !**  INTERFACE
      !**  =========
      !**   INPUT:  BUFFER, IOFUNC, FILE_NAME, PROGRAM_NAME, ROUTINE_NAME, IGRID, JGRID
      !**
      !**   OUTPUT: BUFFER
      !**
      !**  FILES ACCESSED:  Defined in FILE_NAME
      !**
      !**  UPDATES
      !**  =======
      !**  04 Aug 99 Initial AGRMET version..............................Mr Gayno/DNXM
      !**  08 Jun 04 Made the file_name length dynamic.......................AER, Inc.
      !**  24 Feb 10 Converted to FORTRAN 90 for SNODEP...........Mr Lewiston/16WS/WXE
      !**  22 Mar 19 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !**  09 May 19 Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !**
      !*******************************************************************************
      !*******************************************************************************

      ! Defaults
      implicit none

      ! Local constants
      integer,      parameter      :: msglns = 20                ! MAX NUMBER OF LINES IN MESSAGE

      ! Arguments
      real,          intent(inout) :: buffer      (igrid, jgrid) ! DATA TO BE READ/WRITTEN
      character*1,   intent(in)    :: iofunc                     ! I/O FUNCTION ('r', 'w')
      character(len=*), intent(in) :: file_name                  ! FILE PATH AND NAME
      character*12,  intent(in)    :: program_name               ! NAME OF CALLING PROGRAM
      character*12,  intent(in)    :: routine_name               ! NAME OF CALLING ROUTINE
      integer,       intent(in)    :: igrid                      ! SIZE OF GRID IN I-DIRECTION
      integer,       intent(in)    :: jgrid                      ! SIZE OF GRID IN I-DIRECTION

      ! Local variables
      character*7                  :: access_type                ! FILE ACCESS TYPE
      character*4                  :: cstat                      ! I/O STATUS FOR MESSAGE
      character*90                 :: message     (msglns)       ! ERROR MESSAGE
      integer                      :: istat                      ! I/O STATUS
      integer                      :: istat1                     ! CONVERSION STATUS
      integer                      :: reclen                     ! FILE RECORD LENGTH
      logical                      :: isopen                     ! INDICATES IF FILE IS OPEN

      ! Executable code starts here. Open file, abort on error.
      isopen  = .false.
      message = ' '
      reclen  = igrid * jgrid * 4

      access_type = 'OPENING'
      open( 2, file=trim(file_name), form='unformatted',                &
           access='direct', recl=reclen, iostat=istat, err=1000 )
      isopen = .true.

      ! Read from file, abort on error
      if ( (iofunc .eq. 'r') .or. (iofunc .eq. 'R') ) then
         access_type = 'READING'
         read( 2, rec=1, iostat=istat, err=1000 ) buffer

      ! Write to file, abort on error
      else if ( (iofunc .eq. 'w') .or. (iofunc .eq. 'W') ) then
         access_type = 'WRITING'
         write( 2, rec=1, iostat=istat, err=1000 ) buffer

     ! Else abort due to invalid IOFUNC value.
      else
         go to 4000
      endif

      close(2)

      return

      ! Error handling
1000  continue
      if (isopen) close (2)
      message(1) = '[ERR] LIBROUTINE: PUTGET_REAL'
      message(2) = '[ERR] ERROR ' // access_type // ' ' // trim (file_name)
      write (cstat, '(i4)', iostat=istat1) istat
      if (istat1 .eq. 0) message(3) = 'ISTAT = ' // cstat
      call abort_message (program_name, routine_name, message)
      return

4000  continue
      close(2)
      message(1) = '[ERR] LIBROUTINE: PUTGET_REAL'
      message(2) = '[ERR] INVALID IOFUNC VALUE = '// trim(iofunc)
      call abort_message (program_name, routine_name, message)
      return

   end subroutine putget_real

   subroutine tmjul4( hour, day, month, year, julhr ) 

      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      !
      !     Name:  Time/Julian hour conversion with 4-digit year output 
      !     ====
      !
      !     Purpose:
      !     =======
      !     Uses the Julian hour to determine the 2-digit Zulu time, day,  
      !     month, and 4-digit year.    
      !
      !     Method:
      !     ====== 
      !     - Determine the current Zulu hour  
      !     - Determine the total number of elapsed days   
      !     - Count forward to the current day/month/year  
      !  
      !     Interface:
      !     =========  
      !     Inputs:  julhr
      ! 
      !     Iutputs: day, hour, month, year
      !
      !     Variables:
      !     =========
      !     Label       Description
      !     ------      --------------------------------
      !     day         Day of the month (1..31)   
      !     done        Test if done counting days/year
      !     dypmon      Number of days in each month   
      !     elapdy      Number of elapsed days 
      !     hour        Zulu time (UTC) of the julian hour   
      !     i           Loop index 
      !     julhr       Julian hour being processed
      !     month       Month of the year (1..12)  
      !     year        Four-digit year
      !
      !     Files accessed:  None
      !     ==============
      !
      !     Updates:
      !     =======
      !     18 Mar 98 Initial AGRMET version..................SrA Milburn/DNXM
      !     08 Jul 99 Fixed error which initialzed done flag in a data
      !               statement.  Variables must be initialized using an
      !               assignment statement.  Ported to IBM SP2...Mr Gayno/DNXM
      !     15 Sep 09 Incorporated into SNODEP............Mr Lewiston/2WXG/WEA
      !     22 Mar 19 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
      !     09 May 19 Renamed LDTSI...Eric Kemp, NASA GSFC/SSAI
      !  
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------

      ! Defaults
      implicit none 

      ! Arguments
      integer,  intent(out)       :: hour 
      integer,  intent(out)       :: day  
      integer,  intent(out)       :: month
      integer,  intent(out)       :: year 
      integer,  intent(in)        :: julhr

      ! Local variables
      integer                     :: dypmon(12)   
      integer                     :: elapdy   
      integer                     :: i
      logical                     :: done 

      data dypmon   /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      ! Initialize done flag to false.
      done = .false.

      ! Extract the Zulu hour from the Julian hour.   
      hour = mod(julhr, 24) 

      ! Determine the number of days that have elapsed since Dec 31, 1967 
      ! (Julian hour 0).  
      elapdy = julhr / 24   

      ! Initialize the starting day, month, and year values.  
      if (elapdy .gt. 0) then   

         year   = 1968   
         day    = 1  
         month  = 1  
         elapdy = elapdy - 1 

      else  

         year   = 1967   
         day    = 31 
         month  = 12 

      endif

      ! Loop through the elapsed days to determine the year.  

100   continue  

      dypmon(2) = 28  

      ! Determine if year value is a leap year.  Leap years occur every 
      ! 4 years, with the exception of century years not evenly 
      ! divisible by 400.   
      if ((mod(year, 4)   .eq. 0) .and.                                  &   
           ((mod(year, 100) .ne. 0) .or.                                 &   
           (mod(year, 400) .eq. 0))) then  

         dypmon(2) = 29
      endif

      ! If the elapsed number of days is more than a year's worth,  
      ! subtract the appropriate number of days, and increment the year.  
      if (dypmon(2) .eq. 28) then 
         if (elapdy .ge. 365) then 
            year = year + 1 
            elapdy = elapdy - 365   
         else  
            done = .true.   
         endif
      else
         if (elapdy .ge. 366) then 
            year = year + 1 
            elapdy = elapdy - 366   
         else  
            done = .true.   
         endif
      endif

      ! If the elapsed number of days is less than a year's worth, then   
      ! exit loop.
      if (.not. done) goto 100  

      ! Count the days and months elapsed in the current year.
      do 200 i = 1, elapdy  
         day = day + 1   
         if (day .gt. dypmon(month)) then
            day = 1   
            month = month + 1 
         endif
200   end do

      return      
   end subroutine tmjul4

end module LDTSI_utilMod
