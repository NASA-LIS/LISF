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
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
        character (len=LIS_CONST_PATH_LEN):: heading,slofil,zonfil,zfil
  	character (len=LIS_CONST_PATH_LEN):: depfil
  	character (len=LIS_CONST_PATH_LEN):: rizerofil,nxtfil,ndxfil
	character (len=LIS_CONST_PATH_LEN):: wffil,dscfil,init,title
  	character (len=LIS_CONST_PATH_LEN):: elevfil ! added 4/21/2010
  	character (len=LIS_CONST_PATH_LEN), allocatable:: rifil(:)
 	character (len=LIS_CONST_PATH_LEN):: folder,elfoldr
	character (len=8):: suffix

        integer:: RainReadSource ! SY: 1 for forcings coming through LIS interface, 2 for forcing files that native TRIGRS inputs (for testing purposes primarily)

        integer :: i_per ! SY: Counter of TRIGRS time periods for the case RainReadSource = 2 

        real :: TRIGRS_timestep ! SY

end module input_file_defs

