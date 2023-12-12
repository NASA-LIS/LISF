!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readtemplateOpenWatercrd
! \label{readtemplateOpenWatercrd}
!
! !REVISION HISTORY:
! 21 Jul 2004; Sujay Kumar, Initial Code
! 23 Oct 2007; Kristi Arsenault, Updated for use in LISv5.0
!
! !INTERFACE:    
subroutine readtemplateOpenWatercrd()

! !USES:
  use ESMF
  use LIS_coreMod,     only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,      only : LIS_logunit, LIS_verify
  use templateOpenWaterMod, only : templateOpenWater_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to templateOpenWater LSM 
!  option from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n
  character*10  :: time

  call ESMF_ConfigFindLabel(LIS_config,"Template open water timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Template open water timestep: not defined')

     call LIS_parseTimeString(time,templateOpenWater_struc(n)%ts)
  enddo

  write(LIS_logunit,*)'[INFO] Running Template Open water Option:'

end subroutine readtemplateOpenWatercrd
