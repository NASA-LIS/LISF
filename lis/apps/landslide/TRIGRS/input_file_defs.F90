!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module input_file_defs
  	character (len=255):: heading,slofil,zonfil,zfil
  	character (len=255):: depfil
  	character (len=255):: rizerofil,nxtfil,ndxfil
	character (len=255):: wffil,dscfil,init,title
  	character (len=255):: elevfil ! added 4/21/2010
  	character (len=255), allocatable:: rifil(:)
 	character (len=224):: folder,elfoldr
	character (len=8):: suffix

        integer:: RainReadSource ! SY: 1 for forcings coming through LIS interface, 2 for forcing files that native TRIGRS inputs (for testing purposes primarily)

        integer :: i_per ! SY: Counter of TRIGRS time periods for the case RainReadSource = 2 

        real :: TRIGRS_timestep ! SY

end module input_file_defs

