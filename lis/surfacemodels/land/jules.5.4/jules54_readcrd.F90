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
! !ROUTINE: jules54_readcrd
! \label{jules54_readcrd}
!
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES 5.0 
! 28 Nov 2018; Shugong Wang; updated for JULES 5.2 
! 12 Dec 2018; Shugong Wang; updated for JULES 5.4 
!
! !INTERFACE:    
subroutine jules54_readcrd()

! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use jules54_lsmMod, only : jules54_struc
  use jules_lis_exchange
!
! !DESCRIPTION:
!
!  This routine reads the options specific to jules54 LSM 
!  option from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n
  character*10  :: time

  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.4 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'JULES.5.4 model timestep: not defined')

     call LIS_parseTimeString(time,jules54_struc(n)%ts)
  enddo
 
  if (LIS_rc%startcode== "restart" ) then
    call ESMF_ConfigFindLabel(LIS_config,"JULES.5.4 restart file:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, JULES54_struc(n)%rfile, rc=rc)
       call LIS_verify(rc,'JULES.5.4 restart file: not defined')
    enddo
  endif

  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.4 namelist directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, jules54_struc(n)%namelist_dir,rc=rc)
     call LIS_verify(rc,'JULES.5.4 namelist directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "JULES.5.4 restart output interval:", rc = rc)
  do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
      call LIS_verify(rc,"JULES.5.4 restart output interval: not defined")
      call LIS_parseTimeString(time, jules54_struc(n)%rstInterval)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.4 reference height for forcing T and q:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,jules54_struc(n)%z1_tq,rc=rc)
     call LIS_verify(rc,'JULES.5.4 reference height for forcing T and q: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"JULES.5.4 reference height for forcing u and v:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,jules54_struc(n)%z1_uv,rc=rc)
     call LIS_verify(rc,'JULES.5.4 reference height for forcing u and v: not defined')
  enddo

  write(LIS_logunit,*)'Running JULES.5.4 LSM Option:'
  do n=1,LIS_rc%nnest
     jules54_struc(n)%jules54open=0
     jules54_struc(n)%rformat="netcdf"
  enddo

  !!!! jules start time and end time 
  write(main_run_start, '(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') & 
        LIS_rc%syr, '-', LIS_rc%smo, '-', LIS_rc%sda, ' ', LIS_rc%shr, ':', LIS_rc%smn, ':', LIS_rc%sss
  write(main_run_end, '(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') & 
        LIS_rc%eyr, '-', LIS_rc%emo, '-', LIS_rc%eda, ' ', LIS_rc%ehr, ':', LIS_rc%emn, ':', LIS_rc%ess
end subroutine jules54_readcrd
