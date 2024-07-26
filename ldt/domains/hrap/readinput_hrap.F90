!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readinput_hrap
!  \label{readinput_hrap}
! 
! !REVISION HISTORY:
!  15 Jul 2005: Sujay Kumar; Initial Specification
!  15 Dec 2013: KR Arsenault; Updated with latest RDHM code 
!  28 Jan 2014: Shugong Wang; Replace run domain specifications from lat-lon
!
! !INTERFACE:
subroutine readinput_hrap(nest)

! !USES:
  use ESMF
  use LDT_coreMod,   only : LDT_rc, LDT_domain, LDT_config, LDT_localPet, &
       LDT_masterproc, LDT_npes, &
       LDT_ews_halo_ind, LDT_nss_halo_ind, &
       LDT_ewe_halo_ind, LDT_nse_halo_ind
  use LDT_domainMod, only : LDT_quilt_domain, LDT_setParameterDomainSpecs
  use LDT_logMod,    only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
  use LDT_xmrg_reader, only : hrap_to_latlon
  use map_utils

  implicit none

  integer, intent(in) :: nest

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LDT configuration 
!  file to determine the region of interest in the LDT simulation, in 
!  a hrap projection. The routine reads the extents of the 
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest. The domain 
!  specifications of the parameter maps are also read in by this routine. 
!  Based on the number of processors and the processor layout specified, 
!  this routine also performs the domain decomposition. 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[ldttask\_for\_point](\ref{LDT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
!EOP
  integer           :: i, k, n, rc
  integer           :: nc,nr,count
  real              :: lat_str,lon_str
  real              :: stlat,stlon,truelat,stdlon,orient,xmesh
!  real(ESMF_KIND_R8)   :: stlat, stlon, dx, dy
  real              :: ll_lon, ll_lat
  real              :: ur_lon, ur_hrapx
  real              :: ur_lat, ur_hrapy
  real              :: hrap_res
  real, allocatable :: lat(:,:)
  real, allocatable :: lon(:,:)
  real, allocatable :: run_dd(:,:)   ! LIS Target grid/domain
  real, allocatable :: ll_hrapx(:)
  real, allocatable :: ll_hrapy(:)
  type(proj_info)   :: proj_temp
! ________________________________________________________________


  ! Perform LSM check - do not run for other options:
  if( LDT_rc%lsm .ne. "RDHM.3.5.6"    .and. &
      LDT_rc%lsm .ne. "SACHTET.3.5.6" .and. &
      LDT_rc%lsm .ne. "SNOW17" ) then 

     write(LDT_logunit,*) &
        " [ERR] HRAP projection DISABLED for given options at this time!"
     write(LDT_logunit,*) &
        "  -- Please contact the LIS team for further assistance to "
     write(LDT_logunit,*) &
        "      implement required features involving the HRAP grid. "
     call LDT_endrun()

  ! EMK...HRAP support is DISABLED for all other options.
  !  - To later, best handle multi-processor support and subsetting.
  endif

  allocate(run_dd(LDT_rc%nnest,8))
  allocate(ll_hrapx(LDT_rc%nnest))
  allocate(ll_hrapy(LDT_rc%nnest))

  LDT_rc%lis_map_resfactor(nest) = 100.

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left hrap y:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config, ll_hrapy(n), rc=rc)
     call LDT_verify(rc, 'Run domain lower left hrap y: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left hrap x:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config, ll_hrapx(n), rc=rc)
     call LDT_verify(rc, 'Run domain lower left hrap x: not defined')
  enddo

  do n=1, nest
     call hrap_to_latlon(ll_hrapx(n), ll_hrapy(n), ll_lon, ll_lat)
     run_dd(n, 1) = ll_lat
     run_dd(n, 2) = ll_lon
  enddo

  !HRAP true lat = 60 N
  !Run domain true lat
  do n=1,nest
     run_dd(n, 3) = 60.0
  enddo

  !HRAP standard lon = 105.0 West
  !Run domain standard lon
  do n=1,nest
     run_dd(n, 4) = -105.0
  enddo

  ! HRAP orientiation: 0.0
  ! Run domain orientation
  do n=1,nest
     run_dd(n, 5) = 0.0
  enddo

  ! HRAP resolution at true lat: 4.7625 KM when resolution is 1.0 
  call ESMF_ConfigFindLabel(LDT_config,"Run domain hrap resolution:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config, hrap_res,rc=rc)
     run_dd(n, 6) = hrap_res * 4.7625
     call LDT_verify(rc, 'Run domain hrap resolution: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain x-dimension size:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,7),rc=rc)
     call LDT_verify(rc, 'Run domain x-dimension size: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain y-dimension size:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,8),rc=rc)
     call LDT_verify(rc, 'Run domain y-dimension size: not defined')
  enddo

  write(LDT_logunit,*)"DOMAIN details: Running in hrap projection"
  if(LDT_npes.eq.1) then
     LDT_rc%npesx = 1
     LDT_rc%npesy = 1
  endif

  if(LDT_rc%npesx*LDT_rc%npesy.ne.LDT_npes) then
     write(LDT_logunit,*) 'Layout does not match the number of processors...'
     write(LDT_logunit,*) 'npex, npey, ',LDT_rc%npesx,'x',LDT_rc%npesy,'!=',LDT_npes
     write(LDT_logunit,*) 'Stopping program ..'
     call LDT_endrun()
  endif
  
  do n=1,nest
     stlat   = run_dd(n,1)
     stlon   = run_dd(n,2)
     truelat = run_dd(n,3)
     stdlon  = run_dd(n,4)
     orient  = run_dd(n,5)
     xmesh   = run_dd(n,6)
     nc = nint(run_dd(n,7))
     nr = nint(run_dd(n,8))
     
     ur_hrapx = ll_hrapx(n) + (nc-1)*(xmesh/4.7625)
     ur_hrapy = ll_hrapy(n) + (nr-1)*(xmesh/4.7625)

     call hrap_to_latlon(ur_hrapx, ur_hrapy, ur_lon, ur_lat)

     call LDT_quilt_domain(n,nc,nr)

! Call map_set on the overall running domain

     call map_set( PROJ_HRAP, stlat, stlon,&
          xmesh*1000.0, stdlon, truelat,&
          0.0, nc, nr, proj_temp)

     call ij_to_latlon(proj_temp,float(LDT_ews_halo_ind(n,LDT_localPet+1)),&
          float(LDT_nss_halo_ind(n,LDT_localPet+1)), lat_str,lon_str)

     LDT_rc%gnc(n) = nc
     LDT_rc%gnr(n) = nr

     LDT_rc%gridDesc(n,1) = 8
     LDT_rc%gridDesc(n,2) = LDT_rc%lnc(n)
     LDT_rc%gridDesc(n,3) = LDT_rc%lnr(n)
     LDT_rc%gridDesc(n,4) = lat_str
     LDT_rc%gridDesc(n,5) = lon_str
     LDT_rc%gridDesc(n,6) = 8
     LDT_rc%gridDesc(n,7) = orient
     LDT_rc%gridDesc(n,8) = xmesh
     LDT_rc%gridDesc(n,9) = xmesh
     LDT_rc%gridDesc(n,10) = truelat
     LDT_rc%gridDesc(n,11) = stdlon
     
     write(LDT_logunit,*)'-------------------- LDT/LIS Domain ----------------------'
     do k=1,13
        write(unit=LDT_logunit,fmt=*) '(',k,',',LDT_rc%gridDesc(n,k),')'
     enddo

   ! Local grid domain (for parallel components):
     call map_set(PROJ_HRAP,LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,5),&
          LDT_rc%gridDesc(n,8)*1000.0,LDT_rc%gridDesc(n,11),LDT_rc%gridDesc(n,10),&
          0.0,LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_domain(n)%ldtproj)

   ! Added by Shugong Wang to setup the global projection
     call map_set(PROJ_HRAP, stlat, stlon, &
          LDT_rc%gridDesc(n,8)*1000.0,LDT_rc%gridDesc(n,11),LDT_rc%gridDesc(n,10),&
          0.0,LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_domain(n)%ldtglbproj)

     write(LDT_logunit,*)'----------------------------------------------------------'

  enddo
  deallocate(run_dd)
  deallocate(ll_hrapx)
  deallocate(ll_hrapy)

end subroutine readinput_hrap
