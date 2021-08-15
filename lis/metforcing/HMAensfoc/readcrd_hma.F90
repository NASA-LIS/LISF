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
! !ROUTINE: readcrd_hmaens
! \label{readcrd_hmaens}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:    
subroutine readcrd_hmaens()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod
  use hmaens_forcingMod, only : hmaens_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to HMAENS forcing
!  from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,t,rc

  call ESMF_ConfigFindLabel(LIS_config,"HMAENS forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hmaens_struc(n)%hmaensdir,&
          rc=rc)
     call LIS_verify(rc,&
          'HMAENS forcing directory: not defined')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using HMAENS forcing'
     write(LIS_logunit,*) '[INFO] HMAENS forcing directory: ',&
          hmaens_struc(n)%hmaensDIR

     hmaens_struc(n)%hmaenstime1 = 3000.0
     hmaens_struc(n)%hmaenstime2 = 0.0

  enddo
end subroutine readcrd_hmaens
