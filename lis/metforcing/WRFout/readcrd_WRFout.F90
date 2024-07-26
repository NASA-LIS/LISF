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
! !ROUTINE: readcrd_WRFout
!  \label{readcrd_WRFout}
!
! !REVISION HISTORY:
! 14 Mar 2013; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_WRFout()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use WRFout_forcingMod, only : WRFout_struc

  implicit none

! !DESCRIPTION:
!
!  This routine reads the options specific to WRF output forcing from 
!  the LIS configuration file. 
!  
!EOP

  integer :: n,rc
  
  call ESMF_ConfigFindLabel(LIS_config,"WRF output forcing directory:",rc=rc)
  call LIS_verify(rc, 'WRF output forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,WRFout_struc(n)%WRFoutdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"WRF nest id:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,WRFout_struc(n)%nest_id,rc=rc)
  enddo


  write(unit=LIS_logunit,fmt=*)'[INFO] Using WRF output forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] WRF output forcing directory :',trim(WRFout_struc(n)%WRFoutdir)

     WRFout_struc(n)%WRFouttime1 = 3000.0
     WRFout_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine readcrd_WRFout
