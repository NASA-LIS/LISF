!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------

subroutine LVT_init_LDTSIpost()

   ! Imports
   use LVT_logmod, only: LVT_logunit, LVT_endrun, LVT_flush

   ! Defaults
   implicit none

   ! Read the LDTSIpost specific config settings
   call read_ldtsipost_settings()

   ! Other initialize steps
   call LVT_flush(LVT_logunit)

contains

   ! Internal subroutine for reading LDTSIpost-specific config settings
   subroutine read_ldtsipost_settings()

      ! Imports
      use ESMF
      use LVT_coreMod, only: LVT_config, LVT_rc
      use LVT_logMod, only: LVT_endrun, LVT_verify, LVT_logunit

      ! Defaults
      implicit none
      
      ! Local variables
      character(len=255) :: cfgline
      integer :: rc
      
      !cfgline = "LDTSI valid date (YYYYMMDDHH):"
      !call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%yyyymmddhh, &
      !     label=trim(cfgline), rc=rc)
      !call LVT_verify(rc, trim(cfgline)//" not defined")
      write(LVT_rc%yyyymmddhh, '(I4.4,I2.2,I2.2,I2.2)') &
           LVT_rc%syr, LVT_rc%smo, LVT_rc%sda, LVT_rc%shr

      cfgline = "LDTSI output global 0.25deg lat/lon:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_global_ll0p25, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")
      
      cfgline = &
           "LDTSI output NH 16th mesh polar stereographic:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_nh_ps16, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")
      
      cfgline = &
           "LDTSI output SH 16th mesh polar stereographic:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_sh_ps16, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")
      
      if ( .not. LVT_rc%output_global_ll0p25 .and. &
           .not. LVT_rc%output_nh_ps16 .and. &
           .not. LVT_rc%output_sh_ps16) then
         write(LVT_logunit,*) "[ERR] No output selected for LDTSIpost mode!"
         write(LVT_logunit,*) "[ERR] Check the lvt.config file settings!"
         write(LVT_logunit,*) "[ERR] LVT will exit gracefully."
         call LVT_endrun()
      end if
      
      cfgline = "LDTSI input netcdf directory:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%input_dir, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      cfgline = "LDTSI output grib2 directory:"
      call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%output_dir, &
           label=trim(cfgline), rc=rc)
      call LVT_verify(rc, trim(cfgline)//" not defined")

      ! Sanity check earlier lvt.config setting
      if (trim(LVT_rc%lvt_out_format) .ne. "grib2") then
         write(LVT_logunit,*) &
              '[ERR] LVT output format must be set to "grib2" for ' // &
              '"LDTSI post" runmode'
         write(LVT_logunit,*) &
              '[ERR] Instead of "grib2", found "' // &
              trim(LVT_rc%lvt_out_format) // '"'
         write(LVT_logunit,*) '[ERR] Update lvt.config and try again!'
         write(LVT_logunit,*) '[ERR] LVT will stop'
         call LVT_endrun()
      end if
   end subroutine read_ldtsipost_settings
   
end subroutine LVT_init_LDTSIpost
