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
! !ROUTINE: geowrsi2_readcrd
! \label{geowrsi2_readcrd}
!
! !REVISION HISTORY:
! 31 Jul 2011: Brad Wind; Initial Definition
! 14 Jan 2013: KR Arsenault; Updated config file input options
! 25 Oct 2013: KR Arsenault;  Added GeoWRSI2.0 model to LIS-7
!
! !INTERFACE:    
subroutine geowrsi2_readcrd()

! !USES:
  use ESMF
  use LIS_coreMod,     only : LIS_rc, LIS_config
  use LIS_logMod,      only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod,  only : LIS_parseTimeString
  use geowrsi2_lsmMod, only : geowrsi2_struc, num_growing_seasons, &
                              geowrsi2_lsmRunMode
!
! !DESCRIPTION:
!
!  This routine reads the options specific to WRSI LSM 
!  options from the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: rc
  integer :: n
  character*10 :: time
! _____________________

  write(LIS_logunit,*) "Reading in WRSI LSM config file options"

! 'SOS' | 'WRSI' model run mode:
  call ESMF_ConfigFindLabel(LIS_config,"WRSI CalcSOS model run mode:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_lsmRunMode,rc=rc)
  call LIS_verify(rc,'WRSI CalcSOS model run mode: not defined')

! WRSI user input settings file, which is the same as that used by GeoWRSI (VB):
  call ESMF_ConfigFindLabel(LIS_config,"WRSI user input settings file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_struc(n)%inputParmFile%str,rc=rc)
     call LIS_verify(rc,'WRSI user input settings file: not defined')
  enddo

! WRSI crop parameter directory path:
  call ESMF_ConfigFindLabel(LIS_config,"WRSI crop parameter directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_struc(n)%cropDir%str,rc=rc)
     call LIS_verify(rc,'WRSI crop parameter directory: not defined')
  enddo

! WRSI initial and final dekads of a regional growing season (in which a crop can grow):
  call ESMF_ConfigFindLabel(LIS_config,"WRSI initial dekad of season:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_struc(n)%initTstepSeason,rc=rc)
     call LIS_verify(rc,'WRSI initial dekad of season: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"WRSI final dekad of season:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_struc(n)%finalTstepSeason,rc=rc)
     call LIS_verify(rc,'WRSI final dekad of season: not defined')
  enddo

! InitialYear of the first growing season (instead of read-in from UserInputSettings file):
  call ESMF_ConfigFindLabel(LIS_config,"WRSI initial growing season year:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_struc(n)%InitialYear,rc=rc)
     call LIS_verify(rc,'WRSI initial growing season year: not defined')
  enddo
! LastCurrentYear is the end year of the run:
  call ESMF_ConfigFindLabel(LIS_config,"WRSI final growing season year:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_struc(n)%LastCurrentYear,rc=rc)
     call LIS_verify(rc,'WRSI final growing season year: not defined')
  enddo

!- To run multiple growing seasons (added by B.Wind):
  call ESMF_ConfigFindLabel(LIS_config,"WRSI number of growing seasons:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config,num_growing_seasons,rc=rc)
  call LIS_verify(rc,'WRSI number of growing seasons: not defined')
  if( num_growing_seasons < 1 ) num_growing_seasons = 1

! The WRSI model timestep (e.g., 1da, 1hr):
  call ESMF_ConfigFindLabel(LIS_config,"WRSI model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'WRSI model timestep: not defined')
     call LIS_parseTimeString(time,geowrsi2_struc(n)%ts)
  enddo

! Restart file entries:
  if( LIS_rc%startcode == "restart" ) then
     call ESMF_ConfigFindLabel(LIS_config,"WRSI restart output interval:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
        call LIS_verify(rc,'WRSI restart output interval: not defined')
        call LIS_parseTimeString(time, geowrsi2_struc(n)%rstInterval)
     enddo

     geowrsi2_struc%rfile = "none"
     call ESMF_ConfigFindLabel(LIS_config,"WRSI restart file:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,geowrsi2_struc(n)%rfile,rc=rc)
        call LIS_verify(rc,'WRSI restart file: not defined')
     enddo
  endif

end subroutine geowrsi2_readcrd
