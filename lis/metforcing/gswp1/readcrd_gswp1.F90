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
! !ROUTINE: readcrd_gswp1
! \label{readcrd_gswp1}
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Code
!
! !INTERFACE:
subroutine readcrd_gswp1()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use gswp1_forcingMod, only : gswp1_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GSWP1 forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  
  integer   :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"GSWP1 forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gswp1_struc(n)%gswp1dir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'Using GSWP1 forcing'
     write(LIS_logunit,*)'GSWP1 forcing directory: ',trim(gswp1_struc(n)%GSWP1DIR)
  enddo

end subroutine readcrd_gswp1

