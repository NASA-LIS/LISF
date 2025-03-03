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
! !ROUTINE: readcrd_agrradps
! \label{readcrd_agrradps}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_agrradps()
! !USES:
  use agrradps_forcingMod, only : agrradps_struc
  use LIS_coreMod,         only : LIS_config,LIS_rc
  use ESMF 
  use LIS_logMod,          only : LIS_logunit

!
! !DESCRIPTION:
!
!  This routine reads the options specific to AGRRADPS forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  
  integer :: n,rc

  write(LIS_logunit,*)'Using AGRRADPS forcing'
  call ESMF_ConfigFindLabel(LIS_config, "AGRRADPS forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,agrradps_struc(n)%agrpsdir,rc=rc)
    write(LIS_logunit,*) 'AGRRADPS forcing directory :',&
                         trim(agrradps_struc(n)%agrpsdir)
  enddo


end subroutine readcrd_agrradps
