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
! !ROUTINE: snowmodel_readconfig
! \label{snowmodel_readconfig}
!
! !REVISION HISTORY:
!  14 Apr 2020; Kristi Arsenault, Initial Code
!
! !INTERFACE:    
subroutine snowmodel_readconfig()
! !USES:
  use ESMF

  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use snowmodel_lsmMod, only : snowmodel_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to SnowModel
!  from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n,i
  character*10 :: time

! LIS-related SnowModel input entries:

  call ESMF_ConfigFindLabel(LIS_config,"SnowModel model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'SnowModel model timestep: not defined')

     call LIS_parseTimeString(time,snowmodel_struc(n)%ts)
  enddo

  if (LIS_rc%startcode== "restart" ) then
     call ESMF_ConfigFindLabel(LIS_config,"SnowModel restart file:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%rfile,rc=rc)
        call LIS_verify(rc,'SnowModel restart file: not defined')
     enddo
  endif
  call ESMF_ConfigFindLabel(LIS_config,"SnowModel restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'SnowModel restart output interval: not defined')

     call LIS_parseTimeString(time, snowmodel_struc(n)%rstInterval)
  enddo

  ! set default restart format to netcdf
  do n=1,LIS_rc%nnest
     snowmodel_struc(n)%rformat = "netcdf"
  enddo
  ! restart run, read restart file
  if( trim(LIS_rc%startcode) == "restart" ) then
     Call ESMF_ConfigFindLabel(LIS_config, "SnowModel restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, snowmodel_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "SnowModel restart file: not defined")
        enddo

        Call ESMF_ConfigFindLabel(LIS_config, "SnowModel restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, snowmodel_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "SnowModel restart file format: not defined")
        enddo

  ! cold start run, read initial state variables
!  else
  endif

  ! SnowModel-specific parameters:

  call ESMF_ConfigFindLabel(LIS_config,"SnowModel parameter file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%parfile,rc=rc)
     call LIS_verify(rc,'SnowModel parameter file: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Write out SnowModel forcing file fields:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%write_sm_metfields,rc=rc)
     call LIS_verify(rc,'Write out SnowModel forcing file fields: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SnowModel parameters source option:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%sm_params_opt,rc=rc)
     call LIS_verify(rc,'SnowModel parameters source option: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SnowModel preprocess code option:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%call_sm_preproc,rc=rc)
     call LIS_verify(rc,'SnowModel preprocess code option: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SnowModel MicroMet input source:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%sm_micromet_opt,rc=rc)
     call LIS_verify(rc,'SnowModel MicroMet input source: not defined')
  enddo


  call ESMF_ConfigFindLabel(LIS_config,"SnowModel number of snow layers:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%nsnow,rc=rc)
     call LIS_verify(rc,'SnowModel number of snow layers: not defined')
     ! May move this allocate statement later ....
     allocate(snowmodel_struc(n)%lyrthk(snowmodel_struc(n)%nsnow))
  enddo

  ! Initial values:
  call ESMF_ConfigFindLabel(LIS_config,"SnowModel initial snow water equivalent:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%initsnowequiv,rc=rc)
     call LIS_verify(rc,'SnowModel initial snow water equivalent: not defined')
  enddo

! Future expansion to include SnowModel snow layers ...
!  Below is an example from soil layer thicknesses ...
!  call ESMF_ConfigFindLabel(LIS_config,"SnowModel layer thicknesses: ",rc=rc)
!  do n=1,LIS_rc%nnest
!     allocate(snowmodel_struc(n)%lyrthk(snowmodel_struc(n)%nsnow))
!     do i = 1,snowmodel_struc(n)%nsnow
!        call ESMF_ConfigGetAttribute(LIS_config,snowmodel_struc(n)%lyrthk(i),rc=rc)
!     enddo
!     call LIS_verify(rc,'SnowModel layer thicknesses: not defined')
!  enddo


  write(LIS_logunit,*) '[INFO] Running SnowModel (Liston):'

end subroutine snowmodel_readconfig
