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
! !ROUTINE: readinput_gaussian
!  \label{readinput_gaussian}
!
! !REVISION HISTORY:
!
!  11 Apr 2000: Brian Cosgrove; Added Elevation correction and Forcing
!               Mask read statements
!  14 Nov 2003: Sujay Kumar; Modified card file that includes regional
!               modeling options
! !INTERFACE:
subroutine readinput_gaussian(nest)

! !USES:
  use ESMF
  use LDT_coreMod,   only : LDT_rc, LDT_domain, LDT_config
  use LDT_domainMod, only : LDT_quilt_domain
  use LDT_logMod,    only : LDT_logunit, LDT_verify
  use map_utils

  implicit none

  integer, intent(in) :: nest
! !DESCRIPTION:
!
!  This routine reads the options specified in the LIS configuration
!  file to determine the region of interest in the LIS simulation, in
!  a quasi-regular Gaussian lat/lon projection.
!  The routine reads the extents of the
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest.
!  Based on the number of processors and the processor layout specified,
!  this routine also performs the domain decomposition.
!
!  The routines invoked are:
!  \begin{description}
!  \item[ESMF\_ConfigFindLabel] (\ref{ESMF_ConfigFindLabel}) \newline
!    routine to process LDT\_config file
!  \item[ESMF\_ConfigGetAttribute] (\ref{ESMF_ConfigGetAttribute}) \newline
!    routine to process LDT\_config file
!  \item[LDT\_verify] (\ref{LDT_verify}) \newline
!    routine to check return codes
!  \item[gaussian\_comp\_lats] (\ref{gaussian_comp_lats}) \newline
!    routine to compute Gaussian latitudes
!  \item[gaussian\_find\_row] (\ref{gaussian_find_row}) \newline
!    funtion to find index of a given latitude within in a given
!    array of Gaussian latitudes
!  \item[diff\_lon] (\ref{diff_lon}) \newline
!    funtion to compute the difference between two given longitude values
!  \item[LDT\_quilt\_domain] (\ref{LDT_quilt_domain}) \newline
!    routine to quilt the given domain
!  \item[map\_set] (\ref{map_set}) \newline
!    routine to set the map projection data structure for the given domain
!  \end{description}
!EOP

  integer              :: i, j, k, n
  real, allocatable    :: run_dd(:,:)
  real, allocatable    :: gaussian_lat_array(:)
  integer              :: lnc,lnr
  integer              :: nc, nr, s, e
  logical              :: sflag
  integer              :: ierr, rc
  real(ESMF_KIND_R8)   :: dx
  real                 :: diff_lon
  integer              :: ncircles, nlats
  integer              :: gaussian_find_row
! __________________________________________________________________________

  LDT_rc%lis_map_resfactor = 1.

  allocate(run_dd(LDT_rc%nnest,20))

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lat:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,1),rc=rc)
     call LDT_verify(rc, 'Run domain lower left lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lon:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,2),rc=rc)
     call LDT_verify(rc, 'Run domain lower left lon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lat:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,3),rc=rc)
     call LDT_verify(rc, 'Run domain upper right lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lon:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,4),rc=rc)
     call LDT_verify(rc, 'Run domain upper right lon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution dlon:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,5),rc=rc)
     call LDT_verify(rc, 'Run domain resolution dlon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain number of lat circles:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,6),rc=rc)
     call LDT_verify(rc, 'Run domain number of lat circles: not defined')
  enddo

  do n=1,nest
   ! "NUMBER OF LAT CIRCLES" represents the number of latitudes from
   ! a pole to the equator.
   ! The gausslat routine wants the total number of latitudes.
   ncircles = int(run_dd(n,6))
   nlats = ncircles*2

   allocate(gaussian_lat_array(nlats))
   call gaussian_comp_lats(nlats, gaussian_lat_array)

   s = gaussian_find_row(nlats, gaussian_lat_array, run_dd(n,1))
   e = gaussian_find_row(nlats, gaussian_lat_array, run_dd(n,3))
   nr = e - s + 1

   run_dd(n,1) = gaussian_lat_array(s)
   run_dd(n,3) = gaussian_lat_array(e)

   deallocate(gaussian_lat_array)

   ! To avoid problems regarding wrapping around the date line
   ! adjust longitudes to (0, 360).
   !nc = nint(run_dd(n,4)-run_dd(n,2))/run_dd(n,5)
   nc = nint(diff_lon(run_dd(n,4),run_dd(n,2))/run_dd(n,5)) + 1
   dx = run_dd(n,5)

   LDT_rc%gridDesc(n,1)  = 4
   LDT_rc%gridDesc(n,2)  = nc
   LDT_rc%gridDesc(n,3)  = nr
   LDT_rc%gridDesc(n,4)  = run_dd(n,1)
   LDT_rc%gridDesc(n,5)  = run_dd(n,2)
   LDT_rc%gridDesc(n,6)  = 128
   LDT_rc%gridDesc(n,7)  = run_dd(n,3)
   LDT_rc%gridDesc(n,8)  = run_dd(n,4)
   LDT_rc%gridDesc(n,9)  = run_dd(n,5)
   ! The Gaussian projection support in LIS expects the number of circles
   ! in a hemisphere.  Use ncircles. (ncircles=run_dd(n,6))
   LDT_rc%gridDesc(n,10) = ncircles
   LDT_rc%gridDesc(n,20) = 0

   call LDT_quilt_domain(n,nc,nr)

   write(unit=LDT_logunit,fmt=*) 'local domain ', &
      LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,7),  &
      LDT_rc%gridDesc(n,5),LDT_rc%gridDesc(n,8)

   LDT_rc%gnc(n) = nc
   LDT_rc%gnr(n) = nr

   write(LDT_logunit,*)'-------------------- LDT/LIS Domain ----------------------'
   do k=1,13
      write(unit=LDT_logunit,fmt=*) '(',k,',',LDT_rc%gridDesc(n,k),')'
   enddo

   ! map_set wants the total number of latitudes.  Use nlats.
   ! (nlats=2*run_dd(n,6)).
   call map_set(PROJ_GAUSS,LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,5), &
                LDT_rc%gridDesc(n,9),real(nlats),                     &
                LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,7),            &
                LDT_rc%lnc(n),LDT_rc%lnr(n),                          &
                LDT_domain(n)%ldtproj)

   call map_set(PROJ_GAUSS,run_dd(n,1),run_dd(n,2), &
                run_dd(n,5),real(nlats),            &
                run_dd(n,1),run_dd(n,3),            &
                LDT_rc%gnc(n),LDT_rc%gnr(n),        &
                LDT_domain(n)%ldtglbproj)

   write(LDT_logunit,*)'----------------------------------------------------------'

 enddo  ! end nest loop

 deallocate(run_dd)

end subroutine readinput_gaussian
