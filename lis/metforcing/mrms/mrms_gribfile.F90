!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!
! Jonathan Case; 13 Feb 2015: Modified for MRMS file acquisition.
! Jessica Erlingis; 5 September 2017: Modified for MRMS GRIB2 operational data
!
!BOP
! !ROUTINE: mrmsfile
! \label{mrmsfile}
!
! !INTERFACE:
subroutine mrms_gribfile( name, mrms_grib_dir, yr, mo, da, hr)

! USES:
   use LIS_logMod,       only: LIS_logunit

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a MRMS filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[mrms_grib_dir]
!    Name of the MRMS directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name]
!   name of the time-stamped MRMS file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, mo, da, hr

  character(150) :: name
  character(40) :: mrms_grib_dir
  character(4) :: cyear
  character(2) :: cmon, cday, chour

! Operational Convention at NCEP (14 September 2014 - present)
!   input filename (e.g.): 201708/MRMS_GaugeCorr_QPE_01H_00.00_20170831-140000.grib2.gz
!

  write ( cyear, '(i4)' ) yr
  if ( mo < 10 ) write ( cmon, '(a1,i1)' ) "0", mo
  if ( mo >= 10) write ( cmon, '(i2)' ) mo
  if ( da < 10 ) write ( cday, '(a1,i1)' ) "0", da
  if ( da >= 10) write ( cday, '(i2)' ) da
  if ( hr < 10 ) write ( chour, '(a1,i1)' ) "0", hr
  if ( hr >= 10) write ( chour, '(i2)' ) hr

  name = trim(mrms_grib_dir)//'/'//cyear//cmon// &
    '/GRIB2/MRMS_GaugeCorr_QPE_01H_00.00_' &
    //cyear//cmon//cday//'-'//chour//'0000.grib2' !From MSFC Archive

  return

end subroutine mrms_gribfile

