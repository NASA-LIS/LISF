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
! !ROUTINE: vic411_readcard
! \label{vic411_readcard}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
!
! !INTERFACE:    
subroutine vic411_readcard()

! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_endrun
  use vic411_lsmMod,  only : vic411_struc

!
! !DESCRIPTION:
!  This routine reads the options specific to VIC 4.1.1 LSM 
!  option from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n
  character*10 :: time

  call ESMF_ConfigFindLabel(LIS_config,"VIC411 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'VIC411 model timestep: not defined')
     call LIS_parseTimeString(time,vic411_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC411 restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'VIC411 restart output interval: not defined')
     call LIS_parseTimeString(time,vic411_struc(n)%rstInterval)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC411 global parameter file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,vic411_struc(n)%global_param,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC411 veg tiling scheme:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vic411_struc(n)%veg_tiling_scheme,rc=rc)
    if ( vic411_struc(n)%veg_tiling_scheme == 0 .and. &
         LIS_rc%surface_maxt /= 1 ) then
       write(LIS_logunit,*) "ERR: When using VIC's vegetation tiling scheme, LIS must be configured to run at 1 tile per grid."
       call LIS_endrun()
    endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC411 total number of veg types:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,vic411_struc(n)%NT,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC411 convert units:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vic411_struc(n)%debugging_convert_units,rc=rc)
  enddo

  write(LIS_logunit,*)'Running VIC 4.1.1 LSM Option:'

  do n=1,LIS_rc%nnest
     vic411_struc(n)%vicopen=0
  enddo

end subroutine vic411_readcard
