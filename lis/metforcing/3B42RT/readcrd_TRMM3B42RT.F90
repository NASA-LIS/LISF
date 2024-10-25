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
! !ROUTINE: readcrd_TRMM3B42RT
! \label{readcrd_TRMM3B42RT}
!
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_TRMM3B42RT()
! !USES:
  use ESMF 
  use TRMM3B42RT_forcingMod, only : TRMM3B42RT_struc
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod, only : LIS_logunit

!
! !DESCRIPTION:
!
!  This routine reads the options specific to TRMM 3B42RT forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"TRMM 3B42RT forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,TRMM3B42RT_struc(n)%TRMM3B42RTdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'Using TRMM 3B42RT forcing'
     write(LIS_logunit,*) 'TRMM 3B42RT forcing directory :',trim(TRMM3B42RT_struc(n)%TRMM3B42RTDIR)
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
     !TRMM3B42RT_struc(n)%TRMM3B42RTtime = 0.0
  enddo

end subroutine readcrd_TRMM3B42RT
