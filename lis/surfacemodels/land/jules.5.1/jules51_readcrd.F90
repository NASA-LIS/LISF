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
! !ROUTINE: jules51_readcrd
! \label{jules51_readcrd}
!
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES.5.1 
!
! !INTERFACE:    
subroutine jules51_readcrd()

! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use jules51_lsmMod, only : jules51_struc
  use jules_lis_exchange
!
! !DESCRIPTION:
!
!  This routine reads the options specific to jules51 LSM 
!  option from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n
  character*10  :: time

  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.1 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'JULES.5.1 model timestep: not defined')

     call LIS_parseTimeString(time,jules51_struc(n)%ts)
  enddo
 
  if (LIS_rc%startcode== "restart" ) then
    call ESMF_ConfigFindLabel(LIS_config,"JULES.5.1 restart file:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, JULES51_struc(n)%rfile, rc=rc)
       call LIS_verify(rc,'JULES.5.1 restart file: not defined')
    enddo
  endif

  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.1 namelist directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, jules51_struc(n)%namelist_dir,rc=rc)
     call LIS_verify(rc,'JULES.5.1 namelist directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "JULES.5.1 restart output interval:", rc = rc)
  do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
      call LIS_verify(rc,"JULES.5.1 restart output interval: not defined")
      call LIS_parseTimeString(time, jules51_struc(n)%rstInterval)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.1 reference height for forcing T and q:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,jules51_struc(n)%z1_tq,rc=rc)
     call LIS_verify(rc,'JULES.5.1 reference height for forcing T and q: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.1 reference height for forcing u and v:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,jules51_struc(n)%z1_uv,rc=rc)
     call LIS_verify(rc,'JULES.5.1 reference height for forcing u and v: not defined')
  enddo

  write(LIS_logunit,*)'Running JULES.5.1 LSM Option:'
  do n=1,LIS_rc%nnest
     jules51_struc(n)%jules51open=0
     jules51_struc(n)%rformat="netcdf"
  enddo

  !!!! jules start time and end time 
  write(main_run_start, '(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') & 
        LIS_rc%syr, '-', LIS_rc%smo, '-', LIS_rc%sda, ' ', LIS_rc%shr, ':', LIS_rc%smn, ':', LIS_rc%sss
  write(main_run_end, '(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') & 
        LIS_rc%eyr, '-', LIS_rc%emo, '-', LIS_rc%eda, ' ', LIS_rc%ehr, ':', LIS_rc%emn, ':', LIS_rc%ess
end subroutine jules51_readcrd
