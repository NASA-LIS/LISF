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
! !ROUTINE: readcrd_COAMPSout
!  \label{readcrd_COAMPSout}
!
! !REVISION HISTORY:
! 14 Mar 2013; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_COAMPSout()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use COAMPSout_forcingMod, only : COAMPSout_struc

  implicit none

! !DESCRIPTION:
!
!  This routine reads the options specific to COAMPS output forcing from 
!  the LIS configuration file. 
!  
!EOP

  integer :: n,rc
  
  call ESMF_ConfigFindLabel(LIS_config,"COAMPS output forcing directory:",rc=rc)
  call LIS_verify(rc, 'COAMPS output forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,COAMPSout_struc(n)%COAMPSoutdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"COAMPS nest id:",rc=rc)
  call LIS_verify(rc, 'COAMPS nest id: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,COAMPSout_struc(n)%nest_id,rc=rc)
  enddo


  write(unit=LIS_logunit,fmt=*)'[INFO] Using COAMPS output forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] COAMPS output forcing directory :',trim(COAMPSout_struc(n)%COAMPSoutdir)

     COAMPSout_struc(n)%COAMPSouttime1 = 3000.0
     COAMPSout_struc(n)%COAMPSouttime2 = 0.0
  enddo

end subroutine readcrd_COAMPSout
