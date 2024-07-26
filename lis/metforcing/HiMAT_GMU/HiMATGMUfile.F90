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
! !ROUTINE: HiMATGMUfile
! \label{HiMATGMUfile}
!
! !INTERFACE:
subroutine HiMATGMUfile( name, HiMATGMUdir, yr, mo, da, hr)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a STAGE4 filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[HiMATGMUdir]
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

  character(len=*) :: name
  character(len=*) :: HiMATGMUdir
  character(4) :: cyear
  character(2) :: cmon, cday, chour

!- Build the filename to be opened
!   input filename (e.g.): 200501/ST4.2006050604.01h

   write ( cyear, '(i4)' ) yr
   write ( cmon,  '(i2.2)') mo
   write ( cday,  '(i2.2)') da
   write ( chour, '(i2.2)') hr

   name = trim(HiMATGMUdir)//'/'//cyear//cmon//'/gmu'//cyear//cmon//cday//chour//'.nc4'

  return

end subroutine HiMATGMUfile

