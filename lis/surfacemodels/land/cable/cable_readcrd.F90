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
! !ROUTINE: cable_readcrd
! \label{cable_readcrd}
!
! !REVISION HISTORY:
!  21 Jul 2004: Sujay Kumar, Initial Code
!  23 Oct 2007: Kristi Arsenault, Updated for use in LISv5.0
!  25 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!
! !INTERFACE:
subroutine cable_readcrd
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_config
  use LIS_logMod,         only : LIS_logunit, LIS_verify
  use LIS_timeMgrMod,     only : LIS_parseTimeString
  use cable_lsmMod,       only : cable_struc
!
! !DESCRIPTION:
!  This subroutine reads the options specific to CABLE LSM
!  from the LIS configuration file.
!
!EOP
  implicit none
  
  integer :: rc
  integer :: n
  logical :: sliflag
  character*10 :: time
  
  sliflag = .false.
  
  call ESMF_ConfigFindLabel(LIS_config,"CABLE model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'CABLE model timestep: not defined')

     call LIS_parseTimeString(time,cable_struc(n)%ts)
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE restart output interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE restart output interval: not defined')
     call LIS_parseTimeString(time,cable_struc(n)%rstInterval)
  enddo
  
  if (LIS_rc%startcode == "restart") then
     call ESMF_ConfigFindLabel(LIS_config,                         &
          "CABLE restart file:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,                   &
             cable_struc(n)%rfile,rc=rc)
        call LIS_verify(rc,'CABLE restart file: not defined')
     enddo
  endif
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE vegetation parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%vfile,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE vegetation parameter table: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE canopy structure flag:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%canopyflag,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE canopy structure flag: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE photosynthesis structure flag:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%photosynflag,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE photosynthesis structure flag: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE soil structure flag:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%soilflag,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE soil structure flag: not defined')
     if (trim(cable_struc(n)%soilflag).eq.'sli') sliflag = .true.
  enddo
  
  if (sliflag) then
     
     call ESMF_ConfigFindLabel(LIS_config,                            &
          "CABLE sli soils litter structure flag:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,                      &
             cable_struc(n)%slilitterflag,rc=rc)
        call LIS_verify(rc,                                           &
             'CABLE sli soils litter structure flag: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LIS_config,                            &
          "CABLE sli soils isotope structure flag:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,                      &
             cable_struc(n)%sliisotopeflag,rc=rc)
        call LIS_verify(rc,                                           &
             'CABLE sli soils isotope structure flag: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LIS_config,                            &
          "CABLE sli soils coupled structure flag:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,                      &
             cable_struc(n)%slicoupledflag,rc=rc)
        call LIS_verify(rc,                                           &
             'CABLE sli soils coupled structure flag: not defined')
     enddo
     
  endif
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE soil parameter table:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%sfile,rc=rc)
     call LIS_verify(rc,'CABLE soil parameter table: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE fixed vegetation type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%fixedvegtype,rc=rc)
     call LIS_verify(rc,'CABLE fixed vegetation type: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE fixed soil type:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                       &
          cable_struc(n)%fixedsoiltype,rc=rc)
     call LIS_verify(rc,'CABLE fixed soil type: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE fixed snow-free soil albedo:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%fixedalbsoil,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE fixed snow-free soil albedo: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE fixed CO2 concentration:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%fixedco2,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE fixed CO2 concentration: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE reference height:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%refheight,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE reference height: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE maximum verbosity:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%verbose,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE maximum verbosity: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
       "CABLE tile to print:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          cable_struc(n)%tileprint,rc=rc)
     call LIS_verify(rc,                                           &
          'CABLE tile to print: not defined')
  enddo
  
  write(LIS_logunit,*) 'Running CABLE LSM:'
  
end subroutine cable_readcrd
