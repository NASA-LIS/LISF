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
! !ROUTINE: FLake1_readcrd
! \label{FLake1\_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   6/4/13 : Shugong Wang, initial implementation for LIS 7 and FLake1
!
! !INTERFACE:
subroutine FLake1_readcrd()
! !USES:
  use ESMF
  use LIS_coreMod, only    : LIS_rc , LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod, only     : LIS_logunit , LIS_verify
  use FLake1_Mod, only       : FLAKE1_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to FLake1 model from
!  the LIS configuration file.
!
!EOP
  implicit none
  
  integer      :: rc 
  integer      :: n, i
  character*10 :: time 
  
  write(LIS_logunit, *) "Start reading LIS configuration file for FLAKE1 model"
  
  call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 model timestep:", rc = rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
     call LIS_verify(rc, "FLAKE1 model timestep: not defined")
     call LIS_parseTimeString(time, FLAKE1_struc(n)%ts)

     FLAKE1_struc(n)%flake_dt = FLAKE1_struc(n)%ts
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 restart interval:", rc = rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
     call LIS_verify(rc,"FLAKE1 restart interval: not defined")
     call LIS_parseTimeString(time, FLAKE1_struc(n)%rstInterval)
  enddo

  ! height where wind speed is measured
  call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 height where wind speed is measured:", rc = rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%height_wind, rc=rc)
     call LIS_verify(rc, "FLAKE1 height where wind speed is measured: not defiled")
  enddo

  ! height where temperature and humidity are measured
  call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 height where temperature and humidity are measured:", rc = rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%height_tq, rc=rc)
     call LIS_verify(rc, "FLAKE1 height where temperature and humidity are measured: not defiled")
  enddo

  ! restart run, read restart file
  if (trim(LIS_rc%startcode) == "restart") then 
     Call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 restart file:", rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%rfile, rc=rc)
        call LIS_verify(rc, "FLAKE1 restart file: not defined")
     enddo

     !        Call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 restart file format:", rc=rc)
     !        do n=1,LIS_rc%nnest
     !            call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%rformat, rc=rc)
     !            call LIS_verify(rc, "FLAKE1 restart file format: not defined")
     !        enddo
     ! cold start run, read initial state variables
  else 
     ! temperature at the air-snow interface
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial temperature at the air-snow interface:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_T_snow, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial temperature at the air-snow interface: not defined")
     enddo

     ! temperature at the snow-ice interface
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial temperature at the snow-ice interface:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_T_ice, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial temperature at the snow-ice interface: not defined")
     enddo

     ! mean temperature of the water column
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial mean temperature of the water column:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_T_mnw, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial mean temperature of the water column: not defined")
     enddo

     ! temperature of mixed layer
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial temperature of mixed layer:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_T_wML, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial temperature of mixed layer: not defined")
     enddo

     ! temperature at the water-bottom sediment interface
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial temperature at the water-bottom sediment interface:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_T_bot, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial temperature at the water-bottom sediment interface: not defined")
     enddo

     ! temperature at the bottom of the upper layer of the sediments
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial temperature at the bottom of the upper layer of the sediments:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_T_b1, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial temperature at the bottom of the upper layer of the sediments: not defined")
     enddo

     ! thermocline shape factor
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial thermocline shape factor:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_C_T, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial thermocline shape factor: not defined")
     enddo

     ! snow thickness
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial snow thickness:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_H_snow, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial snow thickness: not defined")
     enddo

     ! ice thickness
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial ice thickness:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_H_ice, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial ice thickness: not defined")
     enddo

     ! thickness of mixed layer
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial thickness of mixed layer:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_H_ML, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial thickness of mixed layer: not defined")
     enddo

     ! thickness of the upper layer of bottom sediments
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial thickness of the upper layer of bottom sediments:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_H_B1, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial thickness of the upper layer of bottom sediments: not defined")
     enddo

     ! surface temperature
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial surface temperature:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_T_sfc, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial surface temperature: not defined")
     enddo

     ! water surface albedo with respect to solar radiation
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial water surface albedo with respect to solar radiation:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_albedo_water, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial water surface albedo with respect to solar radiation: not defined")
     enddo

     ! ice surface albedo with respect to the solar radiation
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial ice surface albedo with respect to the solar radiation:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_albedo_ice, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial ice surface albedo with respect to the solar radiation: not defined")
     enddo

     ! snow surface albedo with respect to the solar radiation
     call ESMF_ConfigFindLabel(LIS_config, "FLAKE1 initial snow surface albedo with respect to the solar radiation:", rc = rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, FLAKE1_struc(n)%init_albedo_snow, rc=rc)
        call LIS_verify(rc, "FLAKE1 initial snow surface albedo with respect to the solar radiation: not defined")
     enddo

  end if

  do n=1,LIS_rc%nnest
     FLAKE1_struc(n)%rformat = "netcdf"
  enddo
  write(LIS_logunit, *) "Finished reading LIS configuration file for FLAKE1 model"
     
end subroutine FLAKE1_readcrd
