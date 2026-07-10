!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readconfig_chirps
! \label{readconfig_chirps}
!
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 06 Jul 2015; KR Arsenault, Added support for latest CHIRPS dataset
! 29 June 2025: J Nattala, Updated to support both CHIRPS v2.0 and v3.0
!                          added ESMF config query for CHIRPS version label;
!                          generic forcing directory and resolution labels
! 13 Mar 2026: J Nattala,  Added boundary exclusion lats (CHIRPS3 only)
!
! !INTERFACE:    
subroutine readconfig_chirps()

! !USES:
  use ESMF
  use chirps_forcingMod, only : chirps_struc
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod,  only : LDT_logunit, LDT_endrun, LDT_verify
!
! !DESCRIPTION:
!
!  This routine reads the options specific to CHIRPS forcing 
!  from the LDT configuration file and detects which version to use.
!  
!EOP
  implicit none
  integer :: n, rc, findex
  character(len=20) :: chirps_version_str

   ! Find which forcing index corresponds to CHIRPS
   findex = 0
   do n = 1, LDT_rc%nmetforc
      if (trim(LDT_rc%metforc(n)) == "CHIRPS" .or. &
          trim(LDT_rc%metforc(n)) == "CHIRPS2" .or. &
          trim(LDT_rc%metforc(n)) == "CHIRPS3") then
         findex = n
         exit
      endif
   enddo
   
   if (findex == 0) then
      write(LDT_logunit,*) "[ERR] CHIRPS not found in Met forcing sources"
      call LDT_endrun
   endif
   
   ! Query ESMF config for CHIRPS version label (introduced in unified reader)
   call ESMF_ConfigFindLabel(LDT_config,"CHIRPS version:",rc=rc)
   if (rc == ESMF_SUCCESS) then
      call ESMF_ConfigGetAttribute(LDT_config,chirps_version_str,rc=rc)
      call LDT_verify(rc,"CHIRPS version: attribute not defined")
      
      ! Set version based on the string
      if (trim(chirps_version_str) == "CHIRPS3") then
         do n=1,LDT_rc%nnest
            chirps_struc(n)%version = 3
         enddo
         write(LDT_logunit,*) "[INFO] CHIRPS version: CHIRPS3 (version 3.0)"
      elseif (trim(chirps_version_str) == "CHIRPS2") then
         do n=1,LDT_rc%nnest
            chirps_struc(n)%version = 2
         enddo
         write(LDT_logunit,*) "[INFO] CHIRPS version: CHIRPS2 (version 2.0)"
      else
         write(LDT_logunit,*) "[ERR] Invalid CHIRPS version: ", trim(chirps_version_str)
         write(LDT_logunit,*) "[ERR] Must be 'CHIRPS2' or 'CHIRPS3'"
         call LDT_endrun
      endif
   else
      write(LDT_logunit,*) "[ERR] CHIRPS version: not found in config file"
      write(LDT_logunit,*) "[ERR] Please specify 'CHIRPS version:' as either &
        &'CHIRPS2' or 'CHIRPS3'"
      call LDT_endrun
   endif
   
   ! Query ESMF config for generic forcing directory label (version removed
   ! from label)
   call ESMF_ConfigFindLabel(LDT_config,"CHIRPS forcing directory:",rc=rc)
   if (rc /= ESMF_SUCCESS) then
      write(LDT_logunit,*) "[ERR] CHIRPS forcing directory: not found in config file"
      call LDT_endrun
   endif
   
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,chirps_struc(n)%directory,rc=rc)
      call LDT_verify(rc,"CHIRPS forcing directory: attribute not defined")
      write(LDT_logunit,*) "[INFO] CHIRPS forcing directory: ", &
        trim(chirps_struc(n)%directory)
   enddo
   
   ! Query ESMF config for generic forcing resolution
   call ESMF_ConfigFindLabel(LDT_config,"CHIRPS forcing resolution:",rc=rc)
   if (rc /= ESMF_SUCCESS) then
      write(LDT_logunit,*) "[ERR] CHIRPS forcing resolution: not found in config file"
      call LDT_endrun
   endif
   
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,chirps_struc(n)%xres,rc=rc)
      call LDT_verify(rc,"CHIRPS forcing resolution: attribute not defined")
      chirps_struc(n)%yres = chirps_struc(n)%xres
      write(LDT_logunit,*) "[INFO] CHIRPS forcing resolution: ", &
        chirps_struc(n)%xres
   enddo

   ! Read boundary exclusion latitudes (CHIRPS3 only)
   if( chirps_struc(1)%version == 3 ) then

      call ESMF_ConfigFindLabel(LDT_config,"CHIRPS northern exclusion lat:",rc=rc)
      if (rc /= ESMF_SUCCESS) then
         write(LDT_logunit,*) "[ERR] CHIRPS northern exclusion lat: " // &
           "not found in config file"
         call LDT_endrun
      endif

      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,chirps_struc(n)%north_excl_lat,rc=rc)
         call LDT_verify(rc,"CHIRPS northern exclusion lat: attribute not defined")
         write(LDT_logunit,*) "[INFO] CHIRPS northern exclusion lat: ", &
           chirps_struc(n)%north_excl_lat
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"CHIRPS southern exclusion lat:",rc=rc)
      if (rc /= ESMF_SUCCESS) then
         write(LDT_logunit,*) "[ERR] CHIRPS southern exclusion lat: " // &
           "not found in config file"
         call LDT_endrun
      endif

      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,chirps_struc(n)%south_excl_lat,rc=rc)
         call LDT_verify(rc,"CHIRPS southern exclusion lat: attribute not defined")
         write(LDT_logunit,*) "[INFO] CHIRPS southern exclusion lat: ", &
           chirps_struc(n)%south_excl_lat
      enddo

   endif

   chirps_struc(:)%chirpstime1 = 0.
   chirps_struc(:)%chirpstime2 = 0.

end subroutine readconfig_chirps
