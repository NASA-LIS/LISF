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
! !ROUTINE: vic412_readcard
! \label{vic412_readcard}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
!
! !INTERFACE:    
subroutine vic412_readcard()

! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_endrun
  use vic412_lsmMod,  only : vic412_struc

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

  call ESMF_ConfigFindLabel(LIS_config,"VIC412 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'VIC412 model timestep: not defined')
     call LIS_parseTimeString(time,vic412_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC412 restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'VIC412 restart output interval: not defined')
     call LIS_parseTimeString(time,vic412_struc(n)%rstInterval)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC412 restart file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,vic412_struc(n)%rfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC412 restart file format:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,vic412_struc(n)%rfile_format,rc=rc)
     call LIS_verify(rc,'VIC412 restart file format: not defined (binary or netcdf)')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"VIC412 veg tiling scheme:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vic412_struc(n)%veg_tiling_scheme,rc=rc)
    if ( vic412_struc(n)%veg_tiling_scheme == 0 .and. &
         LIS_rc%surface_maxt /= 1 ) then
       write(LIS_logunit,*) "ERR: When using VIC's vegetation tiling scheme, LIS must be configured to run at 1 tile per grid."
       call LIS_endrun()
    endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC412 total number of veg types:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,vic412_struc(n)%NT,rc=rc)
  enddo

!  call ESMF_ConfigFindLabel(LIS_config,"VIC412 size of state chunk:",rc=rc)
!  if(rc /= 0) then
!    write(LIS_logunit, *) "*******************************************************************"
!    write(LIS_logunit, *) "VIC412 size of state chunk is not defined in LIS configuration file"
!    write(LIS_logunit, *) "Default value (1024) is used. Enlarge the value if number of snow"
!    write(LIS_logunit, *) "band (elevation band) is greater than 10"
!    write(LIS_logunit, *) "*******************************************************************"
!    write(LIS_logunit, *)
!    do n=1,LIS_rc%nnest
!      vic412_struc(n)%state_chunk_size = 1024
!    enddo   
!  else
!    do n=1,LIS_rc%nnest
!       call ESMF_ConfigGetAttribute(LIS_config,vic412_struc(n)%state_chunk_size,rc=rc)
!       call LIS_verify(rc,"VIC412 size of state chunk: not defined")
!    enddo
!  endif

  call ESMF_ConfigFindLabel(LIS_config,"VIC412 convert units:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vic412_struc(n)%debugging_convert_units,rc=rc)
  enddo

  write(LIS_logunit,*)'Running VIC 4.1.2 LSM Option:'

  do n=1,LIS_rc%nnest
     vic412_struc(n)%vicopen=0
  enddo

end subroutine vic412_readcard
