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
! !ROUTINE: readcrd_TRMM3B42RTV7
! \label{readcrd_TRMM3B42RTV7}
!
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 06 Jan 2015; KR Arsenault, Added support for latest V7 data
!
! !INTERFACE:    
subroutine readcrd_TRMM3B42RTV7()

! !USES:
  use ESMF 
  use TRMM3B42RTV7_forcingMod, only : TRMM3B42RTV7_struc
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit
!
! !DESCRIPTION:
!
!  This routine reads the options specific to TRMM 3B42RT V7 forcing 
!   from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"TRMM 3B42RTV7 forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,TRMM3B42RTV7_struc(n)%directory,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'Using TRMM 3B42RT V7 forcing'
     write(LIS_logunit,*) 'TRMM 3B42RT V7 forcing directory :',&
           trim(TRMM3B42RTV7_struc(n)%directory)
  enddo

end subroutine readcrd_TRMM3B42RTV7
