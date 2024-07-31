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
! !ROUTINE: hyssib_readcrd
! \label{hyssib_readcrd}
!
! !REVISION HISTORY:
!  14 Oct 2003: Sujay Kumar, Initial Code
!  21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
!   5 Sep 2007: Chuck Alonge, Updates for LIS 5.0
!  27 Oct 2010: David Mocko, changes for HY-SSiB in LIS6.1
!
! !INTERFACE:
subroutine hyssib_readcrd()
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,    only : LIS_logunit, LIS_verify
  use hyssib_lsmMod, only : hyssib_struc
!
! !DESCRIPTION:
!  Routine to read HY-SSiB specific parameters from the LIS
!  configuration file
!
!EOP
  implicit none
  integer :: rc
  integer :: n
  character*10 :: time

  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'HYSSIB model timestep: not defined')

     call LIS_parseTimeString(time, hyssib_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)

     call LIS_parseTimeString(time, hyssib_struc(n)%rstInterval)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB restart file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%rfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%vfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB albedo parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%afile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB topography stand dev file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%topostdfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB number of vegetation parameters:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%nvegp,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB number of monthly veg parameters:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%nvegip,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB reference height for forcing T and q:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%zh,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB reference height for forcing u and v:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%zm,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB initial soil moisture:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%initsm,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"HYSSIB initial soil temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,hyssib_struc(n)%initTemp,rc=rc)
  enddo

  write(LIS_logunit,*)'Running HY-SSiB LSM:'
  do n=1,LIS_rc%nnest
     hyssib_struc(n)%hyssibopen=0
  enddo

end subroutine hyssib_readcrd

