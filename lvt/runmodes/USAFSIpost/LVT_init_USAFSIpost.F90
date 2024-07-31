!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine LVT_init_USAFSIpost()

   ! Imports
   use LVT_logmod, only: LVT_logunit, LVT_endrun

   ! Defaults
   implicit none

   ! Read the USAFSIpost specific config settings
   call read_usafsipost_settings()

   ! Other initialize steps
   flush(LVT_logunit)

contains

   ! Internal subroutine for reading USAFSIpost-specific config settings
   subroutine read_usafsipost_settings()

      ! Imports
      use ESMF
      use LVT_coreMod, only: LVT_config, LVT_rc
      use LVT_logMod, only: LVT_endrun, LVT_verify, LVT_logunit

      ! Defaults
      implicit none
      
      ! Local variables
      character(len=255) :: cfgline
      integer :: rc
      
      write(LVT_rc%yyyymmddhh, '(I4.4,I2.2,I2.2,I2.2)') &
           LVT_rc%syr, LVT_rc%smo, LVT_rc%sda, LVT_rc%shr

      cfgline = "USAFSI output GRIB2 native grid:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_native, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      cfgline = "USAFSI output GRIB1 global 0.25deg lat/lon:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_global_ll0p25, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")
      
      cfgline = &
           "USAFSI output GRIB1 NH 16th mesh polar stereographic:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_nh_ps16, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")
      
      cfgline = &
           "USAFSI output GRIB1 SH 16th mesh polar stereographic:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_sh_ps16, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      cfgline = &
           "USAFSI output GRIB1 SNODEP NH 16th mesh polar stereographic:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_nh_ps16_snodep, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")
      
      cfgline = &
           "USAFSI output GRIB1 SNODEP SH 16th mesh polar stereographic:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_sh_ps16_snodep, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      
      if ( .not. LVT_rc%output_native .and. &
           .not. LVT_rc%output_global_ll0p25 .and. &
           .not. LVT_rc%output_nh_ps16 .and. &
           .not. LVT_rc%output_sh_ps16 .and. &
           .not. LVT_rc%output_nh_ps16_snodep .and. &
           .not. LVT_rc%output_sh_ps16_snodep) then
         write(LVT_logunit,*) "[ERR] No output selected for USAFpost mode!"
         write(LVT_logunit,*) "[ERR] Check the lvt.config file settings!"
         write(LVT_logunit,*) "[ERR] LVT will exit gracefully."
         call LVT_endrun()
      end if
      
      cfgline = "USAFSI input netcdf directory:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%input_dir, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      cfgline = "USAFSI input netcdf file prefix:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%input_prefix, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      cfgline = "USAFSI output grib directory:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_dir, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      ! Sanity check earlier lvt.config setting
!      if (trim(LVT_rc%lvt_out_format) .ne. "grib2") then
!         write(LVT_logunit,*) &
!              '[ERR] LVT output format must be set to "grib2" for ' // &
!              '"USAFSIpost" runmode'
!         write(LVT_logunit,*) &
!              '[ERR] Instead of "grib2", found "' // &
!              trim(LVT_rc%lvt_out_format) // '"'
!         write(LVT_logunit,*) '[ERR] Update lvt.config and try again!'
!         write(LVT_logunit,*) '[ERR] LVT will stop'
!         call LVT_endrun()
!      end if

      ! Hard code these settings for 557WW
      LVT_rc%security_class = 'U'
      LVT_rc%distribution_class = 'C'
      LVT_rc%data_category = 'C'
      LVT_rc%area_of_data = 'GLOBAL'

   end subroutine read_usafsipost_settings
   
end subroutine LVT_init_USAFSIpost
