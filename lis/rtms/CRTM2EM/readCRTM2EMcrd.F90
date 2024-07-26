!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: readCRTM2EMcrd
! \label{readCRTM2EMcrd}
!
! !REVISION HISTORY:
! 26 Mar 2009; Sujay Kumar, Initial Code
!
! !INTERFACE:    
#include "LIS_misc.h"
subroutine readCRTM2EMcrd()
! !USES:
#if (defined RTMS)
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_verify
  use CRTM2_EMMod, only : crtm_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to CRTM2 from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n

  call ESMF_ConfigFindLabel(LIS_config,"CRTM number of sensors:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%nsensors,rc=rc)
     call LIS_verify(rc,'CRTM number of sensors: not defined')
  enddo

! This should be same as the number of levels in the ATM data. 
  call ESMF_ConfigFindLabel(LIS_config,"CRTM number of layers:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%nlayers,rc=rc)
     call LIS_verify(rc,'CRTM number of layers: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CRTM number of absorbers:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%nabsorbers,rc=rc)
     call LIS_verify(rc,'CRTM number of absorbers: not defined')
  enddo

! not sure if this is the correct qualifying description..
  call ESMF_ConfigFindLabel(LIS_config,"CRTM number of clouds:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%nclouds,rc=rc)
     call LIS_verify(rc,'CRTM number of clouds: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CRTM number of aerosols:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%naerosols,rc=rc)
     call LIS_verify(rc,'CRTM number of aerosols: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CRTM sensor id:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%sensor_id,rc=rc)
     call LIS_verify(rc,'CRTM sensor id: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CRTM coefficient data path:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%coeff_data,rc=rc)
     call LIS_verify(rc,'CRTM coefficient data path: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"CRTM zenith angle:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,crtm_struc(n)%zenith_angle,rc=rc)
     call LIS_verify(rc,'CRTM zenith angle: not defined')
  enddo

#endif
end subroutine readCRTM2EMcrd
