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
! !ROUTINE: readcrd_gdas
! \label{readcrd_gdas}
!
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_gdas()
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_config
  use LIS_logMod
  use gdas_forcingMod, only: gdas_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GDAS forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"GDAS forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdas_struc(n)%gdasdir,rc=rc)
     call LIS_verify(rc,'GDAS forcing directory: not specified')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using GDAS forcing'
     write(LIS_logunit,*) '[INFO] GDAS forcing directory : ',trim(gdas_struc(n)%GDASDIR)
     gdas_struc(n)%GDASTIME1  = 3000.0
     gdas_struc(n)%GDASTIME2  = 0.0
  enddo

end subroutine readcrd_gdas
