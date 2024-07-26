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
! !ROUTINE: template_readcrd
! \label{template_readcrd}
!
! !REVISION HISTORY:
! 21 Jul 2004; Sujay Kumar, Initial Code
! 23 Oct 2007; Kristi Arsenault, Updated for use in LISv5.0
! 23 May 2014: David Mocko, Updated for LIS-7.0
!
! !INTERFACE:    
subroutine template_readcrd()

! !USES:
  use ESMF
  use LIS_coreMod,     only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,      only : LIS_logunit, LIS_verify
  use template_lsmMod, only : template_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to template LSM 
!  option from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n
  character*10  :: time

  call ESMF_ConfigFindLabel(LIS_config,"TEMPLATE model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'TEMPLATE model timestep: not defined')

     call LIS_parseTimeString(time,template_struc(n)%ts)
  enddo

  write(LIS_logunit,*)'Running Template LSM Option:'
  do n=1,LIS_rc%nnest
     template_struc(n)%templateopen=0
  enddo

end subroutine template_readcrd
