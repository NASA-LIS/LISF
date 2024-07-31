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
! !ROUTINE: stg4file
! \label{stg4file}
!
! !INTERFACE:
subroutine stg4file( name, stg4dir, yr, mo, da, hr)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a STAGE4 filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[stg4dir]
!    Name of the STAGE IV directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[name]
!   name of the time-stamped STAGE IV file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, mo, da, hr

  character(80) :: name
  character(40) :: stg4dir
  character(4) :: cyear
  character(2) :: cmon, cday, chour

!- Build the filename to be opened
!   input filename (e.g.): 200501/ST4.2006050604.01h

   write ( cyear, '(i4)' ) yr
   if ( mo < 10 ) write ( cmon, '(a1,i1)' ) "0", mo
   if ( mo >= 10) write ( cmon, '(i2)' ) mo
   if ( da < 10 ) write ( cday, '(a1,i1)' ) "0", da
   if ( da >= 10) write ( cday, '(i2)' ) da
   if ( hr < 10 ) write ( chour, '(a1,i1)' ) "0", hr
   if ( hr >= 10) write ( chour, '(i2)' ) hr

   name = trim(stg4dir)//'/'//cyear//cmon//'/ST4.'//cyear//cmon//cday//chour//'.01h'

  return

end subroutine stg4file

