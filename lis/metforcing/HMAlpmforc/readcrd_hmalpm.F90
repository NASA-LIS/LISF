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
! !ROUTINE: readcrd_hmalpm
! \label{readcrd_hmalpm}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:    
subroutine readcrd_hmalpm()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod
  use hmalpm_forcingMod, only : hmalpm_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to HMALPM forcing
!  from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,t,rc

  call ESMF_ConfigFindLabel(LIS_config,"HMA LPM forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hmalpm_struc(n)%hmalpmdir,&
          rc=rc)
     call LIS_verify(rc,&
          'HMA LPM forcing directory: not defined')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using HMALPM forcing'
     write(LIS_logunit,*) '[INFO] HMALPM forcing directory: ',&
          hmalpm_struc(n)%hmalpmDIR

     hmalpm_struc(n)%hmalpmtime1 = 3000.0
     hmalpm_struc(n)%hmalpmtime2 = 0.0

  enddo
end subroutine readcrd_hmalpm
