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
! !ROUTINE: noah271_readcrd
! \label{noah271_readcrd}
!
! !REVISION HISTORY:
!  14 Oct 2003; Sujay Kumar, Initial Code
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
!
! !INTERFACE:    
subroutine noah271_readcrd()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_warning
  use noah271_lsmMod, only : noah271_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to Noah2.7.1 LSM from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n,i
  character*10 :: time

  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 model timestep: not defined')

     call LIS_parseTimeString(time,noah271_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 restart output interval: not defined')

     call LIS_parseTimeString(time,noah271_struc(n)%rstInterval)
  enddo
  if ( trim(LIS_rc%startcode) == "restart" ) then
     call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 restart file:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%rfile,rc=rc)
        call LIS_verify(rc,'Noah.2.7.1 restart file: not defined')
     enddo
  endif
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%vfile,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 vegetation parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 soil parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%sfile,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 soil parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 use PTF for mapping soil properties:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%useptf,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 use PTF for mapping soil properties: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 number of vegetation parameters:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%nvegp,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 number of vegetation parameters: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 soils scheme:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%soilscheme,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 soils scheme: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 number of soil classes:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%nstxts,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 number of soil classes: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 number of soil layers:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%nslay,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 number of soil layers: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 layer thicknesses: ",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah271_struc(n)%lyrthk(noah271_struc(n)%nslay))
     do i = 1,noah271_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%lyrthk(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.2.7.1 layer thicknesses: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 initial skin temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%initskintemp,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 initial skin temperature: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 initial soil temperatures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah271_struc(n)%inittemp(noah271_struc(n)%nslay))
     do i = 1,noah271_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%inittemp(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.2.7.1 initial soil temperatures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 initial total soil moistures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah271_struc(n)%initsm(noah271_struc(n)%nslay))
     do i = 1,noah271_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%initsm(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.2.7.1 initial total soil moistures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 initial liquid soil moistures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah271_struc(n)%initsmliq(noah271_struc(n)%nslay))
     do i = 1,noah271_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%initsmliq(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.2.7.1 initial liquid soil moistures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 initial canopy water:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%initcanopywater,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 initial canopy water: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 initial snow depth:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%initsnowdepth,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 initial snow depth: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 initial snow equivalent:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%initsnowequiv,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 initial snow equivalent: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 reference height for forcing T and q:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%zh,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 reference height for forcing T and q: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 reference height for forcing u and v:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%zm,rc=rc)
     call LIS_verify(rc,'Noah.2.7.1 reference height for forcing u and v: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"Noah.2.7.1 reinitialize parameters from OPTUE output:",rc=rc)
  do n=1,LIS_rc%nnest
     noah271_struc(n)%param_rst = 0 
     call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%param_rst,rc=rc)
     call LIS_warning(rc,&
          'Noah.2.7.1 reinitialize parameters from OPTUE output: not specified')

     if(noah271_struc(n)%param_rst.eq.1) then 
        call ESMF_ConfigGetAttribute(LIS_config,noah271_struc(n)%prstfile,&
             label="Noah.2.7.1 parameter restart file (from OPTUE):",rc=rc)
        call LIS_verify(rc,&
             "Noah.2.7.1 parameter restart file (from OPTUE): not defined")
     endif

  enddo

  
  write(LIS_logunit,*) 'Running Noah2.7.1 LSM:'
  do n=1,LIS_rc%nnest
     noah271_struc(n)%noah271open=0
     noah271_struc(n)%nsoilp   = 10
  enddo
end subroutine noah271_readcrd
