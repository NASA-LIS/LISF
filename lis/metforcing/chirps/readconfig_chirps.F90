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
! 28 Jun 2026: J Nattala, Updated to support both CHIRPS v2.0 and v3.0;
!                         added ESMF config query for CHIRPS version label;
!                         generic forcing directory and resolution labels;
!                         added boundary exclusion lats (CHIRPS3 only)
!
! !INTERFACE:
subroutine readconfig_chirps()

! !USES:
  use ESMF
  use chirps_forcingMod, only : chirps_struc
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_endrun, LIS_verify
!
! !DESCRIPTION:
!
!  This routine reads options specific to CHIRPS forcing from the
!  LIS config and detects the applicable version
!
!EOP
  implicit none
  integer :: n, rc, findex
  character(len=20) :: chirps_version_str

   ! Find which forcing index corresponds to CHIRPS
   findex = 0
   do n = 1, LIS_rc%nmetforc
      if (trim(LIS_rc%metforc(n)) == "CHIRPS" .or. &
          trim(LIS_rc%metforc(n)) == "CHIRPS2" .or. &
          trim(LIS_rc%metforc(n)) == "CHIRPS3") then
         findex = n
         exit
      endif
   enddo

   if (findex == 0) then
      write(LIS_logunit,*) "[ERR] CHIRPS not found in Met forcing sources"
      call LIS_endrun
   endif

   ! Query ESMF config for CHIRPS version label (introduced in unified reader)
   call ESMF_ConfigFindLabel(LIS_config,"CHIRPS version:",rc=rc)
   if (rc == ESMF_SUCCESS) then
      call ESMF_ConfigGetAttribute(LIS_config,chirps_version_str,rc=rc)
      call LIS_verify(rc,"CHIRPS version: attribute not defined")

      ! Set version based on the string
      if (trim(chirps_version_str) == "CHIRPS3") then
         do n=1,LIS_rc%nnest
            chirps_struc(n)%version = 3
         enddo
         write(LIS_logunit,*) "[INFO] CHIRPS version: CHIRPS3 (version 3.0)"
      elseif (trim(chirps_version_str) == "CHIRPS2") then
         do n=1,LIS_rc%nnest
            chirps_struc(n)%version = 2
         enddo
         write(LIS_logunit,*) "[INFO] CHIRPS version: CHIRPS2 (version 2.0)"
      else
         write(LIS_logunit,*) "[ERR] Invalid CHIRPS version: ", trim(chirps_version_str)
         write(LIS_logunit,*) "[ERR] Must be 'CHIRPS2' or 'CHIRPS3'"
         call LIS_endrun
      endif
   else
      write(LIS_logunit,*) "[ERR] CHIRPS version: not found in config file"
      write(LIS_logunit,*) "[ERR] Please specify 'CHIRPS version:' as either &
        &'CHIRPS2' or 'CHIRPS3'"
      call LIS_endrun
   endif

   ! Query ESMF config for generic forcing directory label (version removed
   ! from label)
   call ESMF_ConfigFindLabel(LIS_config,"CHIRPS forcing directory:",rc=rc)
   if (rc /= ESMF_SUCCESS) then
      write(LIS_logunit,*) "[ERR] CHIRPS forcing directory: not found in config file"
      call LIS_endrun
   endif

   do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,chirps_struc(n)%directory,rc=rc)
      call LIS_verify(rc,"CHIRPS forcing directory: attribute not defined")
      write(LIS_logunit,*) "[INFO] CHIRPS forcing directory: ", &
        trim(chirps_struc(n)%directory)
   enddo

   ! Read generic resolution
   call ESMF_ConfigFindLabel(LIS_config,"CHIRPS forcing resolution:",rc=rc)
   if (rc /= ESMF_SUCCESS) then
      write(LIS_logunit,*) "[ERR] CHIRPS forcing resolution: not found in config file"
      call LIS_endrun
   endif

   do n=1,LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,chirps_struc(n)%xres,rc=rc)
      call LIS_verify(rc,"CHIRPS forcing resolution: attribute not defined")
      chirps_struc(n)%yres = chirps_struc(n)%xres
      write(LIS_logunit,*) "[INFO] CHIRPS forcing resolution: ", &
        chirps_struc(n)%xres
   enddo

   ! Read boundary exclusion latitudes (CHIRPS3 only)
   if( chirps_struc(1)%version == 3 ) then

      call ESMF_ConfigFindLabel(LIS_config,"CHIRPS northern exclusion lat:",rc=rc)
      if (rc /= ESMF_SUCCESS) then
         write(LIS_logunit,*) "[ERR] CHIRPS northern exclusion lat: " // &
           "not found in config file"
         call LIS_endrun
      endif

      do n=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,chirps_struc(n)%north_excl_lat,rc=rc)
         call LIS_verify(rc,"CHIRPS northern exclusion lat: attribute not defined")
         write(LIS_logunit,*) "[INFO] CHIRPS northern exclusion lat: ", &
           chirps_struc(n)%north_excl_lat
      enddo

      call ESMF_ConfigFindLabel(LIS_config,"CHIRPS southern exclusion lat:",rc=rc)
      if (rc /= ESMF_SUCCESS) then
         write(LIS_logunit,*) "[ERR] CHIRPS southern exclusion lat: " // &
           "not found in config file"
         call LIS_endrun
      endif

      do n=1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,chirps_struc(n)%south_excl_lat,rc=rc)
         call LIS_verify(rc,"CHIRPS southern exclusion lat: attribute not defined")
         write(LIS_logunit,*) "[INFO] CHIRPS southern exclusion lat: ", &
           chirps_struc(n)%south_excl_lat
      enddo

   endif

   chirps_struc(:)%chirpstime1 = 0.
   chirps_struc(:)%chirpstime2 = 0.

end subroutine readconfig_chirps
