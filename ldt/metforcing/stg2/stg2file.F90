!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! !ROUTINE: stg2file
! \label{stg2file}
!
! !INTERFACE:
subroutine stg2file( filename, stg2dir, yr, mo, da, hr)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a STAGE2 filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[stg2dir]
!    Name of the STAGE II directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   name of the time-stamped STAGE II file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, mo, da, hr

  character(len=*) :: filename
  character(len=*) :: stg2dir
  character(4)  :: cyear
  character(2)  :: cmon, cday, chour

!- Build the filename to be opened
!   input filename (e.g.): 200501/ST2ml2005011912.Grb

   write ( cyear, '(i4)' ) yr
   if ( mo < 10 ) write ( cmon, '(a1,i1)' ) "0", mo
   if ( mo >= 10) write ( cmon, '(i2)' ) mo
   if ( da < 10 ) write ( cday, '(a1,i1)' ) "0", da
   if ( da >= 10) write ( cday, '(i2)' ) da
   if ( hr < 10 ) write ( chour, '(a1,i1)' ) "0", hr
   if ( hr >= 10) write ( chour, '(i2)' ) hr

   filename = trim(stg2dir)//'/'//cyear//cmon//'/ST2ml'//&
              cyear//cmon//cday//chour//'.Grb'

end subroutine stg2file

