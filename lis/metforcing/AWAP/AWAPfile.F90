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
! !ROUTINE: AWAPfile
! \label{AWAPfile}
!
! !INTERFACE:
subroutine AWAPfile( name, AWAPdir, yr, doy)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a STAGE4 filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[AWAPdir]
!    Name of the STAGE IV directory
!  \item[yr]
!    year 
!  \item[doy]
!   day of year
!   \item[name]
!   name of the AWAP  file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, doy

  character(len=*) :: name
  character(len=*) :: AWAPdir
  character(4) :: cyear
  character(3) :: cdoy

   write ( cyear, '(i4)' ) yr
   write ( cdoy, '(i3)' ) doy

   name = trim(AWAPdir)//'/'//trim(cyear)//'/'//&
        trim(cdoy)//'/awap.grb'

end subroutine AWAPfile

