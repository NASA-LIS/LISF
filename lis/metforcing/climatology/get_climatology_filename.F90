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
! !ROUTINE: get_climatology_filename
! \label{get_climatology_filename}
!
! !REVISION HISTORY:
!  27 Sep 2016: K. Arsenault; Initial Implementation
!
! !INTERFACE:
 subroutine get_climatology_filename(doy,directory,filename)

   implicit none
! !ARGUMENTS: 
   integer, intent(in)        :: doy             ! File day of year (DOY)
   character(len=*), intent(in)  :: directory       ! File directory
   character(len=*), intent(out) :: filename  
!
! !DESCRIPTION:
!   This subroutine puts together LDT generated forcing
!    file names.
! 
!  The arguments are:
!  \begin{description}
!   \item[doy]
!     Integer-based day of year (DOY) input
!   \item[directory]
!     Name of the forcing directory
!   \item[filename]
!     Name of the time-stamped LDT-generated climatology forcing file
!  \end{description}
!
!EOP
   character*3  :: fdoy
 
  !=== end variable definition =============================================

   write(unit=fdoy, fmt='(i3.3)') doy

  !=== Assemble LDT-generated climatology filename:
  !   ./DIRECTORY/FORCING/LDT_FORC_CLIMO_088.nc

   filename = trim(directory)//"/FORCING/LDT_FORC_CLIMO_"&
              //fdoy//".nc"

 end subroutine get_climatology_filename

