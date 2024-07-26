!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: readinput_hrap
!  \label{readinput_hrap}
! 
! !REVISION HISTORY:
!  15 Jul 2005: Sujay Kumar; Initial Specification
!  28 Jan 2014: Shugong Wang; Replace run domain specifications from lat-lon
!               into HRAP X and Y 
! 
! !INTERFACE:
subroutine readinput_hrap

! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_domain, LIS_config, LIS_localPet, &
       LIS_masterproc, LIS_npes, &
       LIS_ews_halo_ind, LIS_nss_halo_ind, LIS_ewe_halo_ind, LIS_nse_halo_ind
  use LIS_domainMod, only : LIS_quilt_domain, LIS_setParameterDomainSpecs
  use LIS_logMod,    only : LIS_logunit, LIS_endrun, LIS_verify
  use LIS_fileIOMod, only : LIS_readDomainConfigSpecs
  use LIS_XMRG_Reader, only : hrap_to_latlon
  use map_utils
  use netcdf 

  implicit none

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LIS configuration 
!  file to determine the region of interest in the LIS simulation, in 
!  a hrap projection. The routine reads the extents of the 
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest. The domain 
!  specifications of the parameter maps are also read in by this routine. 
!  Based on the number of processors and the processor layout specified, 
!  this routine also performs the domain decomposition. 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[listask\_for\_point](\ref{LIS_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
!EOP
  integer           :: i, k, n, rc
  integer           :: nc,nr,count
  real              :: lat_str,lon_str
  real              :: stlat,stlon,truelat,stdlon,orient,xmesh
  real, allocatable :: run_dd(:,:)   ! LIS Target grid/domain
  type(proj_info)   :: proj_temp
  real, allocatable :: ll_hrapx(:)
  real, allocatable :: ll_hrapy(:)
  real              :: ll_lon
  real              :: ll_lat
  real              :: ur_lon, ur_hrapx
  real              :: ur_lat, ur_hrapy
  real              :: hrap_res
  real, allocatable     :: lat(:,:)
  real, allocatable     :: lon(:,:)
  integer           :: ios
  integer           :: ftn
  logical           :: file_exists
  integer           :: latId, lonId
  integer           :: ncId, nrId
! ________________________________________________________________

  allocate(run_dd(LIS_rc%nnest,8))
  allocate(ll_hrapx(LIS_rc%nnest))
  allocate(ll_hrapy(LIS_rc%nnest))

  call ESMF_ConfigFindLabel(LIS_config,"Run domain lower left hrap y:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, ll_hrapy(n), rc=rc)
     call LIS_verify(rc, 'Run domain lower left hrap y: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain lower left hrap x:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, ll_hrapx(n), rc=rc)
     call LIS_verify(rc, 'Run domain lower left hrap x: not defined')
  enddo

  do n=1, LIS_rc%nnest
     call hrap_to_latlon(ll_hrapx(n), ll_hrapy(n), ll_lon, ll_lat)
     run_dd(n, 1) = ll_lat
     run_dd(n, 2) = ll_lon
  enddo

  !HRAP true lat = 60 N
  !Run domain true lat
  do n=1,LIS_rc%nnest
     run_dd(n, 3) = 60.0 
  enddo

  !HRAP standard lon = 105.0 West
  !Run domain standard lon
  do n=1,LIS_rc%nnest
     run_dd(n, 4) = -105.0
  enddo

  ! HRAP orientiation: 0.0
  ! Run domain orientation
  do n=1,LIS_rc%nnest
     run_dd(n, 5) = 0.0
  enddo

  ! HRAP resolution at true lat: 4.7625 KM when resolution is 1.0 
  call ESMF_ConfigFindLabel(LIS_config,"Run domain hrap resolution:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, hrap_res,rc=rc)
     run_dd(n, 6) = hrap_res * 4.7625 
     call LIS_verify(rc, 'Run domain hrap resolution: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain x-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,7),rc=rc)
     call LIS_verify(rc, 'Run domain x-dimension size: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain y-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,8),rc=rc)
     call LIS_verify(rc, 'Run domain y-dimension size: not defined')
  enddo

  write(LIS_logunit,*)'DOMAIN details: Running in hrap projection'
  if(LIS_npes.eq.1) then 
     LIS_rc%npesx = 1
     LIS_rc%npesy = 1
  endif
  
  if(LIS_rc%npesx*LIS_rc%npesy.ne.LIS_npes) then
     write(LIS_logunit,*) 'Layout does not match the number of processors...'
     write(LIS_logunit,*) 'npex, npey, ',LIS_rc%npesx,'x',LIS_rc%npesy,'!=',LIS_npes
     write(LIS_logunit,*) 'Stopping program ..'
     call LIS_endrun()
  endif
  
  do n=1,LIS_rc%nnest
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
     
     call LIS_quilt_domain(n,nc,nr)

! Call map_set on the overall running domain

     call map_set(PROJ_HRAP,stlat,stlon,&
          xmesh*1000.0,stdlon,truelat,&
          0.0,nc,nr,proj_temp)

     call ij_to_latlon(proj_temp,float(LIS_ews_halo_ind(n,LIS_localPet+1)),&
          float(LIS_nss_halo_ind(n,LIS_localPet+1)), lat_str,lon_str)

     LIS_rc%gnc(n) = nc
     LIS_rc%gnr(n) = nr
    
     LIS_rc%pnc(n) = nc
     LIS_rc%pnr(n) = nr
     ! added by Shugong Wang for reading lat-lon
     allocate(lat(LIS_rc%gnc(n), LIS_rc%gnr(n)))
     allocate(lon(LIS_rc%gnc(n), LIS_rc%gnr(n)))
     
     inquire(file=LIS_rc%paramfile(n), exist=file_exists)
     if(file_exists) then 
         ios = nf90_open(path=LIS_rc%paramfile(n),&
              mode=NF90_NOWRITE,ncid=ftn)
         call LIS_verify(ios,'Error in nf90_open in readinput_hrap')
         
         ios = nf90_inq_dimid(ftn,"east_west",ncId)
         call LIS_verify(ios,'Error in nf90_inq_dimid in readinput_hrap')
         
         ios = nf90_inq_dimid(ftn,"north_south",nrId)
         call LIS_verify(ios,'Error in nf90_inq_dimid in readinput_hrap')
         
         ios = nf90_inquire_dimension(ftn,ncId, len=LIS_rc%gnc(n))
         call LIS_verify(ios,'Error in nf90_inquire_dimension in readinput_hrap')
         
         ios = nf90_inquire_dimension(ftn,nrId, len=LIS_rc%gnr(n))
         call LIS_verify(ios,'Error in nf90_inquire_dimension in readinput_hrap')
        
         ios = nf90_inq_varid(ftn,'lat',latid)
         call LIS_verify(ios,'lat field not found in the LIS param file')
         
         ios = nf90_inq_varid(ftn,'lon',lonid)
         call LIS_verify(ios,'lon field not found in the LIS param file')
         
         ios = nf90_get_var(ftn,latid,lat)
         call LIS_verify(ios,'Error in nf90_get_var for lat in readinput_hrap')
         
         ios = nf90_get_var(ftn,lonid,lon)
         call LIS_verify(ios,'Error in nf90_get_var for lon in readinput_hrap')
         
         if(abs(run_dd(n,1)-lat(1,1)).gt.0.0001.or.&
            abs(run_dd(n,2)-lon(1,1)).gt.0.0001.or.&
            abs(ur_lat-lat(LIS_rc%gnc(n),LIS_rc%gnr(n))).gt.0.0001.or.&
            abs(ur_lon-lon(LIS_rc%gnc(n),LIS_rc%gnr(n))).gt.0.0001) then 
            write(LIS_logunit,*) run_dd(n, 1), lat(1,1), run_dd(n, 2), lon(1,1)
            write(LIS_logunit,*) 'The selected running domain is outside the '
            write(LIS_logunit,*) 'domain provided through the parameter attributes file '
            write(LIS_logunit,*) 'program stopping ...'
            call LIS_endrun
         endif
        
        LIS_rc%gridDesc(n,34) = lat(1,1)
        LIS_rc%gridDesc(n,35) = lon(1,1)
        LIS_rc%gridDesc(n,37) = lat(LIS_rc%gnc(n), LIS_rc%gnr(n))
        LIS_rc%gridDesc(n,38) = lon(LIS_rc%gnc(n), LIS_rc%gnr(n))
        
        deallocate(lat)
        deallocate(lon)
     else
        write(LIS_logunit,*) LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit,*) 'program stopping ...'
        call LIS_endrun
     endif
     
     LIS_rc%gridDesc(n,1) = 8
     LIS_rc%gridDesc(n,2) = LIS_rc%lnc(n)
     LIS_rc%gridDesc(n,3) = LIS_rc%lnr(n)
     LIS_rc%gridDesc(n,4) = lat_str
     LIS_rc%gridDesc(n,5) = lon_str
     LIS_rc%gridDesc(n,6) = 8
     LIS_rc%gridDesc(n,7) = orient
     LIS_rc%gridDesc(n,8) = xmesh
     LIS_rc%gridDesc(n,9) = xmesh
     LIS_rc%gridDesc(n,10) = truelat
     LIS_rc%gridDesc(n,11) = stdlon
     
     write(LIS_logunit,*)'-------------------- LIS/LIS Domain ----------------------'
     do k=1,13
        write(unit=LIS_logunit,fmt=*) '(',k,',',LIS_rc%gridDesc(n,k),')'
     enddo

     call map_set(PROJ_HRAP,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
          LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),&
          0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)

     ! added by Shugong Wang to setup the global projection
     call map_set(PROJ_HRAP, stlat, stlon, &
          LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),&
          0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisglbproj)
     ! added by Shugong Wang to setup the parameter projection
     call map_set(PROJ_HRAP, LIS_rc%gridDesc(n,34),LIS_rc%gridDesc(n,35), &
          LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),&
          0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisparamproj)
  enddo
  deallocate(run_dd)
  deallocate(ll_hrapx)
  deallocate(ll_hrapy)
end subroutine readinput_hrap
