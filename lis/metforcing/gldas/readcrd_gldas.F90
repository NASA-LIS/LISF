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
! !ROUTINE: readcrd_gldas
! \label{readcrd_gldas}
!
!
! !REVISION HISTORY:
!  19 Sept 2008: Sujay Kumar: Initial Implementation
!
! !INTERFACE:    
subroutine readcrd_gldas()
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_config
  use LIS_logMod,       only : LIS_logunit
  use gldas_forcingMod, only : gldas_struc
!
! !DESCRIPTION:
!
!  This routine reads the runtime options specific to GLDAS forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"GLDAS forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gldas_struc(n)%gldasdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using GLDAS forcing'
     write(LIS_logunit,*) 'GLDAS forcing directory :',trim(gldas_struc(n)%GLDASDIR)
     gldas_struc(n)%GLDASTIME1  = 3000.0
     gldas_struc(n)%GLDASTIME2  = 0.0
  enddo

end subroutine readcrd_gldas
