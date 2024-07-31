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
! !ROUTINE: noah33_readcrd
! \label{noah33_readcrd}
!
! !REVISION HISTORY:
!  14 Oct 2003; Sujay Kumar, Initial Code
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!
! !INTERFACE:    
subroutine noah33_readcrd()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,  only : LIS_logunit, LIS_verify
  use noah33_lsmMod, only : noah33_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to Noah3.3 LSM from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n,i
  character*10 :: time

  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Noah.3.3 model timestep: not defined')

     call LIS_parseTimeString(time,noah33_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Noah.3.3 restart output interval: not defined')

     call LIS_parseTimeString(time, noah33_struc(n)%rstInterval)
  enddo
  if (LIS_rc%startcode== "restart" ) then
     call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 restart file:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%rfile,rc=rc)
        call LIS_verify(rc,'Noah.3.3 restart file: not defined')
     enddo
  endif
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%vfile,rc=rc)
     call LIS_verify(rc,'Noah.3.3 vegetation parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 soil parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%sfile,rc=rc)
     call LIS_verify(rc,'Noah.3.3 soil parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 general parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%gfile,rc=rc)
     call LIS_verify(rc,'Noah.3.3 general parameter table: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 use PTF for mapping soil properties:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%useptf,rc=rc)
     call LIS_verify(rc,'Noah.3.3 use PTF for mapping soil properties: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 soils scheme:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%soilscheme,rc=rc)
     call LIS_verify(rc,'Noah.3.3 soils scheme: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 number of soil layers:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%nslay,rc=rc)
     call LIS_verify(rc,'Noah.3.3 number of soil layers: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 layer thicknesses: ",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah33_struc(n)%lyrthk(noah33_struc(n)%nslay))
     do i = 1,noah33_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%lyrthk(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.3 layer thicknesses: not defined')

  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 use distributed soil depth map: ",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%usedsoilmap,&
          default=0,rc=rc)
!     call LIS_verify(rc,'Noah.3.3 use distributed soil depth map: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 use distributed root depth map: ",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%usedrootmap,&
          default=0,rc=rc)
!     call LIS_verify(rc,'Noah.3.3 use distributed root depth map: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 initial skin temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%initskintemp,rc=rc)
     call LIS_verify(rc,'Noah.3.3 initial skin temperature: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 initial soil temperatures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah33_struc(n)%inittemp(noah33_struc(n)%nslay))
     do i = 1,noah33_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%inittemp(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.3 initial soil temperatures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 initial total soil moistures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah33_struc(n)%initsm(noah33_struc(n)%nslay))
     do i = 1,noah33_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%initsm(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.3 initial total soil moistures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 initial liquid soil moistures:",rc=rc)
  do n=1,LIS_rc%nnest
     allocate(noah33_struc(n)%initsmliq(noah33_struc(n)%nslay))
     do i = 1,noah33_struc(n)%nslay
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%initsmliq(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.3 initial liquid soil moistures: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 initial canopy water:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%initcanopywater,rc=rc)
     call LIS_verify(rc,'Noah.3.3 initial canopy water: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 initial snow depth:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%initsnowdepth,rc=rc)
     call LIS_verify(rc,'Noah.3.3 initial snow depth: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 initial snow equivalent:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%initsnowequiv,rc=rc)
     call LIS_verify(rc,'Noah.3.3 initial snow equivalent: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 fixed max snow albedo:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%fixedmxsnalb,rc=rc)
     call LIS_verify(rc,'Noah.3.3 fixed max snow albedo: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 fixed deep soil temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%fixedtbot,rc=rc)
     call LIS_verify(rc,'Noah.3.3 fixed deep soil temperature: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 fixed vegetation type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%fixedvegtype,rc=rc)
     call LIS_verify(rc,'Noah.3.3 fixed vegetation type: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 fixed soil type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%fixedsoiltype,rc=rc)
     call LIS_verify(rc,'Noah.3.3 fixed soil type: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 fixed slope type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%fixedslopetype,rc=rc)
     call LIS_verify(rc,'Noah.3.3 fixed slope type: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 sfcdif option:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%sfcdifoption,rc=rc)
     call LIS_verify(rc,'Noah.3.3 sfcdif option: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 z0 veg-type dependence option:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%iz0tlnd,rc=rc)
     call LIS_verify(rc,'Noah.3.3 z0 veg-type dependence option: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 greenness fraction:",rc=rc)
  do n=1,LIS_rc%nnest
     do i = 1,12
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%shdfac_monthly(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.3 greenness fraction: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 background albedo:",rc=rc)
  do n=1,LIS_rc%nnest
     do i = 1,12
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%albedo_monthly(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.3 background albedo: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 background roughness length:",rc=rc)
  do n=1,LIS_rc%nnest
     do i = 1,12
        call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%z0brd_monthly(i),rc=rc)
     enddo
     call LIS_verify(rc,'Noah.3.3 background roughness length: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 reference height for forcing T and q:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%zh,rc=rc)
     call LIS_verify(rc,'Noah.3.3 reference height for forcing T and q: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Noah.3.3 reference height for forcing u and v:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,noah33_struc(n)%zm,rc=rc)
     call LIS_verify(rc,'Noah.3.3 reference height for forcing u and v: not defined')
  enddo

  write(LIS_logunit,*) '[INFO] Running Noah3.3 LSM:'

end subroutine noah33_readcrd
