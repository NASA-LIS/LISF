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
! !ROUTINE: noah32_readcrd
! \label{noah32_readcrd}
!
! !REVISION HISTORY:
!  14 Oct 2003; Sujay Kumar, Initial Code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!
! !INTERFACE:    
subroutine noah32_readcrd()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,  only : LIS_logunit, LIS_verify
  use noah32_lsmMod, only : noah32_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to Noah3.2 LSM from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n,i
  character*10 :: time

  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Noah.3.2 model timestep: not defined')

     call LIS_parseTimeString(time,noah32_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Noah.3.2 restart output interval: not defined')

     call LIS_parseTimeString(time, noah32_struc(n)%rstInterval)
  enddo
  if (LIS_rc%startcode== "restart" ) then
     call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 restart file:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%rfile,rc=rc)
        call LIS_verify(rc,'Noah.3.2 restart file: not defined')
     enddo
  endif
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%vfile,rc=rc)
     call LIS_verify(rc,'Noah.3.2 vegetation parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 soil parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%sfile,rc=rc)
     call LIS_verify(rc,'Noah.3.2 soil parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 general parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%gfile,rc=rc)
     call LIS_verify(rc,'Noah.3.2 general parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 use PTF for mapping soil properties:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%useptf,rc=rc)
     call LIS_verify(rc,'Noah.3.2 use PTF for mapping soil properties: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 soils scheme:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%soilscheme,rc=rc)
     call LIS_verify(rc,'Noah.3.2 soils scheme: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 number of soil layers:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%nslay,rc=rc)
     call LIS_verify(rc,'Noah.3.2 number of soil layers: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 layer thicknesses: ",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah32_struc(n)%lyrthk(noah32_struc(n)%nslay))
     do i = 1,noah32_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%lyrthk(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.2 layer thicknesses: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 use distributed soil depth map: ",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%usedsoilmap,&
          default=0,rc=rc)
!     call LIS_verify(rc,'Noah.3.2 use distributed soil depth map: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 use distributed root depth map: ",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%usedrootmap,&
          default=0,rc=rc)
!     call LIS_verify(rc,'Noah.3.2 use distributed root depth map: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 initial skin temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%initskintemp,rc=rc)
     call LIS_verify(rc,'Noah.3.2 initial skin temperature: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 initial soil temperatures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah32_struc(n)%inittemp(noah32_struc(n)%nslay))
     do i = 1,noah32_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%inittemp(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.2 initial soil temperatures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 initial total soil moistures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah32_struc(n)%initsm(noah32_struc(n)%nslay))
     do i = 1,noah32_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%initsm(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.2 initial total soil moistures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 initial liquid soil moistures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah32_struc(n)%initsmliq(noah32_struc(n)%nslay))
     do i = 1,noah32_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%initsmliq(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.2 initial liquid soil moistures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 initial canopy water:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%initcanopywater,rc=rc)
     call LIS_verify(rc,'Noah.3.2 initial canopy water: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 initial snow depth:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%initsnowdepth,rc=rc)
     call LIS_verify(rc,'Noah.3.2 initial snow depth: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 initial snow equivalent:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%initsnowequiv,rc=rc)
     call LIS_verify(rc,'Noah.3.2 initial snow equivalent: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 fixed max snow albedo:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%fixedmxsnalb,rc=rc)
     call LIS_verify(rc,'Noah.3.2 fixed max snow albedo: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 fixed deep soil temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%fixedtbot,rc=rc)
     call LIS_verify(rc,'Noah.3.2 fixed deep soil temperature: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 fixed vegetation type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%fixedvegtype,rc=rc)
     call LIS_verify(rc,'Noah.3.2 fixed vegetation type: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 fixed soil type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%fixedsoiltype,rc=rc)
     call LIS_verify(rc,'Noah.3.2 fixed soil type: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 fixed slope type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%fixedslopetype,rc=rc)
     call LIS_verify(rc,'Noah.3.2 fixed slope type: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 sfcdif option:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%sfcdifoption,rc=rc)
     call LIS_verify(rc,'Noah.3.2 sfcdif option: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 z0 veg-type dependence option:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%iz0tlnd,rc=rc)
     call LIS_verify(rc,'Noah.3.2 z0 veg-type dependence option: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 greenness fraction:",rc=rc)
  do n=1,LIS_rc%nnest
     do i = 1,12
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%shdfac_monthly(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.2 greenness fraction: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 background albedo:",rc=rc)
  do n=1,LIS_rc%nnest
     do i = 1,12
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%albedo_monthly(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.2 background albedo: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 background roughness length:",rc=rc)
  do n=1,LIS_rc%nnest
     do i = 1,12
        call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%z0brd_monthly(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.2 background roughness length: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 reference height for forcing T and q:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%zh,rc=rc)
     call LIS_verify(rc,'Noah.3.2 reference height for forcing T and q: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.2 reference height for forcing u and v:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah32_struc(n)%zm,rc=rc)
     call LIS_verify(rc,'Noah.3.2 reference height for forcing u and v: not defined')
  enddo

  write(LIS_logunit,*) 'Running Noah3.2 LSM:'

end subroutine noah32_readcrd
