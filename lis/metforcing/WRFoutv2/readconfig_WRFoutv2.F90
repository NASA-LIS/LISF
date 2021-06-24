!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readconfig_WRFoutv2
!  \label{readconfig_WRFoutv2}
!
! !REVISION HISTORY:
! 14 Mar 2013; Sujay Kumar, Initial Code
! 20 Nov 2020; K.R. Arsenault, Updated for different WRF output files
!
! !INTERFACE:    
subroutine readconfig_WRFoutv2()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use WRFoutv2_forcingMod, only : WRFoutv2_struc

  implicit none

! !DESCRIPTION:
!
!  This routine reads the options specific to WRF output forcing from 
!  the LIS configuration file. 
!  
!EOP

  integer :: n,rc
  
  call ESMF_ConfigFindLabel(LIS_config,"WRF output v2 forcing directory:",rc=rc)
  call LIS_verify(rc, 'WRF output v2 forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,WRFoutv2_struc(n)%WRFoutv2dir,rc=rc)
  enddo

!  call ESMF_ConfigFindLabel(LIS_config,"WRF nest id:",rc=rc)
  do n=1,LIS_rc%nnest    
     WRFoutv2_struc(n)%nest_id = 1
!     call ESMF_ConfigGetAttribute(LIS_config,WRFoutv2_struc(n)%nest_id,rc=rc)
  enddo


  write(unit=LIS_logunit,fmt=*)'[INFO] Using WRF output v2 forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] WRF output v2 forcing directory :',WRFoutv2_struc(n)%WRFoutv2dir

     WRFoutv2_struc(n)%WRFouttime1 = 3000.0
     WRFoutv2_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine readconfig_WRFoutv2
