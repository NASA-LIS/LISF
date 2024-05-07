!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine readConfigFile(configfile)

! !USES:
  use ESMF 
  use preprocMod
!
! !DESCRIPTION:
!
!EOP
  implicit none
  
  character(len=*), intent(in) :: configfile

  integer       :: rc
  integer       :: n, i
  character*100 :: temp1
  character*10  :: lis_map_proj
  character*1   :: fproc(4)
!____________________________________________________________

  LVT_config = ESMF_ConfigCreate(rc=rc)
  call LVT_verify(rc,'problem in creating LVT_config object')

  call ESMF_ConfigLoadFile(LVT_config,trim(configfile),rc=rc)
  call LVT_verify(rc,'problem in loading preproc.config')

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%lis_map_proj,&
       label="Map projection of the LIS domain:",rc=rc)
  call LVT_verify(rc,'Map projection of the LIS domain: option not specified in the config file')

  LVT_rc%nnest = 1
  allocate(LVT_rc%gridDesc(LVT_rc%nnest,20))
  allocate(LVT_rc%mfile(LVT_rc%nnest))
  allocate(LVT_rc%vfile(LVT_rc%nnest))
  allocate(LVT_rc%vfile_form(LVT_rc%nnest))
  LVT_rc%udef = -9999.

!- Assign values for surface type model options:
  LVT_rc%max_model_types = 3
  LVT_rc%lsm_index = 1
  LVT_rc%lake_index = 2
  LVT_rc%glacier_index = 3
  
  allocate(LVT_rc%sf_model_type_name(LVT_rc%max_model_types))
  allocate(LVT_rc%sf_model_type(LVT_rc%max_model_types))
  
  LVT_rc%sf_model_type_name(1) = "LSM"
  LVT_rc%sf_model_type(1)      = LVT_rc%lsm_index
  LVT_rc%sf_model_type_name(2) = "Lake"
  LVT_rc%sf_model_type(2)      = LVT_rc%lake_index
  LVT_rc%sf_model_type_name(3) = "Glacier"
  LVT_rc%sf_model_type(3)      = LVT_rc%glacier_index


  allocate(LVT_rc%gnc(LVT_rc%nnest))
  allocate(LVT_rc%gnr(LVT_rc%nnest))
  allocate(LVT_rc%lnc(LVT_rc%nnest))
  allocate(LVT_rc%lnr(LVT_rc%nnest))
  allocate(LVT_rc%lnc_red(LVT_rc%nnest))
  allocate(LVT_rc%lnr_red(LVT_rc%nnest))
  allocate(LVT_domain(LVT_rc%nnest))

  if(LVT_rc%lis_map_proj.eq."latlon") then 
     call readinput_latlon()
  elseif(LVT_rc%lis_map_proj.eq."lambert") then 
     call readinput_lambert()
  elseif(LVT_rc%lis_map_proj.eq."polar") then 
     call readinput_polar()
  endif
  
end subroutine readConfigFile

