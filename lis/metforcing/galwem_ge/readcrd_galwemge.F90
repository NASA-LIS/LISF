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
! !ROUTINE: readcrd_galwemge
! \label{readcrd_galwemge}
!
! !REVISION HISTORY:
! 09 May 2022; Yeosang Yoon, Initial Code
!
! !INTERFACE:    
subroutine readcrd_galwemge()
! !USES:
  use LIS_logMod
  use LIS_coreMod
  use galwemge_forcingMod, only : galwemge_struc
  use ESMF
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GALWEM-GE forecast forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM-GE forecast forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwemge_struc(n)%odir,rc=rc)
     call LIS_verify(rc,'GALWEM-GE forecast forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM-GE forecast run mode:",rc=rc)
  call LIS_verify(rc, 'GALWEM-GE forecast run mode: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwemge_struc(n)%runmode,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM-GE forecast number of ensemble members:",rc=rc)
  call LIS_verify(rc, 'GALWEM-GE forecast number of ensemble members: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwemge_struc(n)%max_ens_members,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using GALWEM-GE forecast forcing'
     write(LIS_logunit,*) '[INFO] GALWEM-GE forecast forcing directory: ', trim(galwemge_struc(n)%odir)
     write(LIS_logunit,*) '[INFO] GALWEM-GE forecast run mode: ',galwemge_struc(n)%runmode
     write(LIS_logunit,*) '[INFO] GALWEM-GE forecast number of ensemble members:',&
           galwemge_struc(n)%max_ens_members
  enddo
end subroutine readcrd_galwemge
