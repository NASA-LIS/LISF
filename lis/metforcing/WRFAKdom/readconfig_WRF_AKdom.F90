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
!
! !ROUTINE: readconfig_WRF_AKdom
!  \label{readconfig_WRF_AKdom}
!
! !REVISION HISTORY:
!  21 Jun 2021: K.R. Arsenault; Updated for different WRF AK files
!
! !INTERFACE:    
subroutine readconfig_WRF_AKdom()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use WRF_AKdom_forcingMod, only : WRFAK_struc

  implicit none

! !DESCRIPTION:
!
!  This routine reads the options specific to WRF output forcing from 
!  the LIS configuration file. 
!  
!EOP

  integer :: n,rc
  
  call ESMF_ConfigFindLabel(LIS_config,"WRF AK forcing directory:",rc=rc)
  call LIS_verify(rc, 'WRF AK forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,WRFAK_struc(n)%WRFAKdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest    
     WRFAK_struc(n)%nest_id = 1
  enddo

  write(unit=LIS_logunit,fmt=*)'[INFO] Using WRF AK forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] WRF AK forcing directory :',trim(WRFAK_struc(n)%WRFAKdir)

     WRFAK_struc(n)%WRFouttime1 = 3000.0
     WRFAK_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine readconfig_WRF_AKdom
