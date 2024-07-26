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
! !ROUTINE: readmoscrd
! \label{readmoscrd}
!
! !REVISION HISTORY:
! 15 Oct 2003; Sujay Kumar, Initial Code
! 25 Sep 2007: Sujay Kumar, Upgraded for LIS 5.0
! 12 Nov 2007: Chuck Alonge, Added entry # of soil textures
!
! !INTERFACE:    
subroutine readmoscrd()
! !USES:
  use ESMF 
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod, only : LIS_logunit, LIS_verify
  use LIS_coreMod, only : LIS_rc, LIS_config
  use mos_lsmMod, only : mos_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to Mosaic LSM from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc
  character*10 :: time

  write(LIS_logunit,*) 'Reading mosaic card file..'

  write(LIS_logunit,*)'Running MOS LSM:'
  
  call ESMF_ConfigFindLabel(LIS_config,"Mosaic model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'Mosaic model timestep: not defined')

     call LIS_parseTimeString(time,mos_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)


     call LIS_parseTimeString(time,mos_struc(n)%rstInterval)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic restart file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%mos_rfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%mos_vfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic monthly vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%mos_mvfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic soil parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%mos_sfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic number of soil classes:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%mos_nstxts,rc=rc)
  enddo

#if 0 
  call ESMF_ConfigFindLabel(LIS_config,"Mosaic Depth of Layer 1 (m):",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%dpthlyr1,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic Depth of Layer 2 (m):",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%dpthlyr2,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic Depth of Layer 3 (m):",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%dpthlyr3,rc=rc)
  enddo
#endif

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic initial soil moisture:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%mos_ism,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic initial soil temperature:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%mos_it,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "Mosaic use forcing data observation height:", rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%forcing_z,default=0,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "Mosaic use forcing data aerodynamic conductance:" , rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%forcing_ch,default=0,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Mosaic use distributed soil depth map: ",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,mos_struc(n)%usedsoilmap,&
          default=0,rc=rc)
!     call LIS_verify(rc,'MOS use distributed soil depth map: not defined')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'MOS Active Restart File: ', mos_struc(n)%MOS_RFILE
     mos_struc(n)%MOSopen = 0
     mos_struc(n)%mos_nvegp = 24
     mos_struc(n)%mos_nmvegp = 6
     mos_struc(n)%mos_nsoilp = 10
  enddo

end subroutine readmoscrd
