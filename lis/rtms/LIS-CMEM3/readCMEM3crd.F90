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
! !ROUTINE: readCMEM3crd
! \label{readCMEM3crd}
!
! !REVISION HISTORY:
! 21 Feb 2011; Yudong Tian: added support for freq-dependent zenith angles (e.g., WindSat)
!
! !INTERFACE:    
#include "LIS_misc.h"
subroutine readCMEM3crd()
! !USES:

#if (defined RTMS)
  use ESMF
  use CMEM3_Mod, only : cmem3_struc
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_logMod

!
! !DESCRIPTION:
!
!  This routine reads the options specific to CMEM3 from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: ftn 
  integer :: n, j
  real :: freqs(100), theta(100)
  character(len=LIS_CONST_PATH_LEN) :: fname(100)

  call ESMF_ConfigFindLabel(LIS_config,"CMEM3 sensor id:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,cmem3_struc(n)%sensor_id,rc=rc)
     call LIS_verify(rc,'CMEM3 sensor id: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CMEM3 number of frequencies:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,cmem3_struc(n)%nfs,rc=rc)
     call LIS_verify(rc,'CMEM3 number of frequencies: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CMEM3 frequencies file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, cmem3_struc(n)%freqfile,rc=rc)
     call LIS_verify(rc,'CMEM3 frequencies file: not defined')
  End Do
  
  ftn = LIS_getNextUnitNumber()
  open(ftn, file=cmem3_struc(1)%freqfile, form="formatted")  
  do j=1, cmem3_struc(1)%nfs
     read(ftn, *) freqs(j), theta(j), fname(j)
  end do
  call LIS_releaseUnitNumber(ftn)
    
  do n=1,LIS_rc%nnest
     do j=1, cmem3_struc(n)%nfs
        cmem3_struc(n)%fghz(j) = freqs(j) 
        cmem3_struc(n)%theta(j) = theta(j) 
        cmem3_struc(n)%chName(j) = fname(j)
     end do
  end do
  
#endif
end subroutine readCMEM3crd
