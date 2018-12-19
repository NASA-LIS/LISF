!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! !ROUTINE: AWRALfile
! \label{AWRALfile}
!
! !INTERFACE:
subroutine AWRALfile( name, AWRALdir, yr, doy)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a STAGE4 filename for 
!   one hour file intervals
! 
!  The arguments are:
!  \begin{description}
!  \item[AWRALdir]
!    Name of the STAGE IV directory
!  \item[yr]
!    year 
!  \item[doy]
!   day of year
!   \item[name]
!   name of the AWRAL  file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, doy

  character(80) :: name
  character(40) :: AWRALdir
  character(4) :: cyear
  character(3) :: cdoy

   write ( cyear, '(i4)' ) yr
   write ( cdoy, '(i3)' ) doy

   name = trim(AWRALdir)//'/AWRAL_'//trim(cyear)//'/'//&
        trim(cdoy)
!   name = trim(AWRALdir)//'/AWRAL_'//trim(cyear)//'/'//&
!        trim(cdoy)//'



end subroutine AWRALfile

