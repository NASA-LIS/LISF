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
! !ROUTINE: clsmf25_readcrd
! \label{clsmf25_readcrd}
!
! !DESCRIPTION:
!  Routine to read Catchment specific parameters from the card file. 
!
! !REVISION HISTORY:
! 12 Feb 2006; Sujay Kumar, Initial Code
! 23 Nov 2012: David Mocko, Additional configs for Catchment Fortuna-2.5
!
! !INTERFACE:    
subroutine clsmf25_readcrd()
! !USES:
  use clsmf25_lsmMod, only : clsmf25_struc
  use LIS_logMod,     only : LIS_logunit,LIS_verify
  use LIS_timeMgrMod
  use LIS_coreMod,    only : LIS_config,LIS_rc
  use ESMF
!EOP
  implicit none

  integer      :: n, rc
  character*20 :: time
  integer      :: flxfix

!  allocate(LIS_rc%tile_coord_file(LIS_rc%nnest))
!  allocate(LIS_rc%tile_veg_file(LIS_rc%nnest))

  call ESMF_ConfigFindLabel(LIS_config,"CLSM F2.5 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 model timestep: not defined')

     call LIS_parseTimeString(time,clsmf25_struc(n)%ts)
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 restart output interval: not defined')
     call LIS_parseTimeString(time,clsmf25_struc(n)%rstInterval)

  enddo
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 restart file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clsmf25_struc(n)%rfile,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 restart file: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 top soil layer depth:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clsmf25_struc(n)%dzsfcrd,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 top soil layer depth: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 initial soil moisture:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clsmf25_struc(n)%initSM,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 initial soil moisture: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 initial soil temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clsmf25_struc(n)%initST,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 initial soil temperature: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 fixed reference height:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clsmf25_struc(n)%RefHfixed,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 fixed reference height: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 turbulence scheme:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clsmf25_struc(n)%turbscheme,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 turbulence scheme: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CLSM F2.5 use MODIS albedo flag:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,clsmf25_struc(n)%usemodisalbflag,rc=rc)
     call LIS_verify(rc,'CLSM F2.5 use MODIS albedo flag: not defined')
  enddo

  write(LIS_logunit,*) '[INFO] Running Catchment Fortuna-2.5 LSM ...'
  
end subroutine clsmf25_readcrd
