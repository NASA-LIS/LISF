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
! !ROUTINE: readcrd_galwem
! \label{readcrd_galwem}
!
! !REVISION HISTORY:
! 11 Mar 2022; Yeosang Yoon, Initial Code
! 08 Sep 2022; Yeosang Yoon, Add codes to read GALWEM 25 DEG dataset
!
! !INTERFACE:    
subroutine readcrd_galwem()
! !USES:
  use LIS_logMod
  use LIS_coreMod
  use galwem_forcingMod, only : galwem_struc
  use ESMF
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GALWEM forecast forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM forecast forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwem_struc(n)%odir,rc=rc)
     call LIS_verify(rc,'GALWEM forecast forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM forecast run mode:",rc=rc)
  call LIS_verify(rc, 'GALWEM forecast run mode: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwem_struc(n)%runmode,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM forecast resolution:",rc=rc)
  call LIS_verify(rc, 'GALWEM forecast resolution: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwem_struc(n)%resol,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using GALWEM forecast forcing'
     write(LIS_logunit,*) '[INFO] GALWEM forecast forcing directory: ', trim(galwem_struc(n)%odir)
     write(LIS_logunit,*) '[INFO] GALWEM forecast run mode: ',galwem_struc(n)%runmode
     write(LIS_logunit,*) '[INFO] GALWEM forecast resolution: ',galwem_struc(n)%resol
  enddo
end subroutine readcrd_galwem
