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
! !ROUTINE : readclm2crd
! \label{readclm2crd}
!
! !REVISION HISTORY:
! 14 Oct 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readclm2crd()
! !USES:
  use ESMF 
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_verify
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use clm2_lsmMod, only : clm2_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to CLM from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc
  character*10 :: time

  call ESMF_ConfigFindLabel(LIS_config,"CLM model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'CLM model timestep: not defined')

     call LIS_parseTimeString(time,clm2_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CLM restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)

     call LIS_parseTimeString(time,clm2_struc(n)%rstInterval)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CLM restart file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clm2_struc(n)%clm_rfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CLM vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clm2_struc(n)%clm_vfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CLM canopy height table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clm2_struc(n)%clm_chtfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CLM initial soil moisture:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clm2_struc(n)%clm_ism,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CLM initial soil temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clm2_struc(n)%clm_it,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CLM initial snow mass:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clm2_struc(n)%clm_iscv,rc=rc)
     write(LIS_logunit,*)'Running CLM LSM:'
     write(LIS_logunit,*)'CLM Active Restart File: ', clm2_struc(n)%clm_rfile
     clm2_struc(n)%clmopen=0
  enddo

end subroutine readclm2crd
