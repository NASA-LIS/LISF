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
! !ROUTINE: readinput_polar
!  \label{readinput_polar}
! 
! !REVISION HISTORY:
!  15 Jul 2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine readinput_polar
! !USES:
  use ESMF 
  use map_utils
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_config, &
       LIS_localPet, LIS_masterproc, LIS_npes, &
       LIS_ews_halo_ind, LIS_nss_halo_ind, LIS_ewe_halo_ind, LIS_nse_halo_ind
  use LIS_domainMod, only : LIS_quilt_domain, LIS_setParameterDomainSpecs
  use LIS_fileIOMod, only : LIS_readDomainConfigSpecs
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !DESCRIPTION: 
!
!  This routine reads the options specified in the LIS configuration 
!  file to determine the region of interest in the LIS simulation, in 
!  a polar stereographic projection. The routine reads the extents of the 
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

  integer           :: i
  integer           :: k
  real, allocatable :: run_dd(:,:)
  real              :: lat_str,lon_str
  real              :: stlat,stlon,truelat,stdlon,orient,xmesh
  integer           :: nc,nr
  type(proj_info)   :: proj_temp
  integer           :: n , rc
  integer              :: ios
  integer              :: ftn
  logical              :: file_exists
  integer              :: latId, lonId
  integer              :: ncId, nrId
  integer, allocatable     :: gnc(:), gnr(:)
  real, allocatable        :: lat(:,:)
  real, allocatable        :: lon(:,:)

  allocate(run_dd(LIS_rc%nnest,8))
  allocate(gnc(LIS_rc%nnest))
  allocate(gnr(LIS_rc%nnest))

  call ESMF_ConfigFindLabel(LIS_config,"Run domain lower left lat:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,1),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain lower left lon:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,2),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Run domain true lat:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,3),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain standard lon:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,4),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain orientation:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,5),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain resolution:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,6),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain x-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,7),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain y-dimension size:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(n,8),rc=rc)
  enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  do n=1,LIS_rc%nnest
     inquire(file=LIS_rc%paramfile(n), exist=file_exists)
     if(file_exists) then 
        
        ios = nf90_open(path=LIS_rc%paramfile(n),&
             mode=NF90_NOWRITE,ncid=ftn)
        call LIS_verify(ios,'Error in nf90_open in readinput_latlon')
        
        ios = nf90_inq_dimid(ftn,"east_west",ncId)
        call LIS_verify(ios,'Error in nf90_inq_dimid in readinput_latlon')
        
        ios = nf90_inq_dimid(ftn,"north_south",nrId)
        call LIS_verify(ios,'Error in nf90_inq_dimid in readinput_latlon')
        
        ios = nf90_inquire_dimension(ftn,ncId, len=gnc(n))
        call LIS_verify(ios,'Error in nf90_inquire_dimension in readinput_latlon')
        
        ios = nf90_inquire_dimension(ftn,nrId, len=gnr(n))
        call LIS_verify(ios,'Error in nf90_inquire_dimension in readinput_latlon')
        
        LIS_rc%pnc(n) = gnc(n)
        LIS_rc%pnr(n) = gnr(n)

        allocate(lat(gnc(n),gnr(n)))
        allocate(lon(gnc(n),gnr(n)))
        
        ios = nf90_inq_varid(ftn,'lat',latid)
        call LIS_verify(ios,'lat field not found in the LIS param file')
        
        ios = nf90_inq_varid(ftn,'lon',lonid)
        call LIS_verify(ios,'lon field not found in the LIS param file')
        
        ios = nf90_get_var(ftn,latid,lat)
        call LIS_verify(ios,'Error in nf90_get_var for lat in readinput_latlon')
        
        ios = nf90_get_var(ftn,lonid,lon)
        call LIS_verify(ios,'Error in nf90_get_var for lon in readinput_latlon')
        
        if(abs(run_dd(n,1)-lat(1,1)).gt.0.0001.or.&
             abs(run_dd(n,3)-lat(gnc(n),gnr(n))).gt.0.0001.or.&
             abs(run_dd(n,2)-lon(1,1)).gt.0.0001.or.&
             abs(run_dd(n,4)-lon(gnc(n),gnr(n))).gt.0.0001) then 
           write(LIS_logunit,*) 'The selected running domain is outside the '
           write(LIS_logunit,*) 'domain provided through the parameter attributes file '
           write(LIS_logunit,*) 'program stopping ...'
           call LIS_endrun
        endif
        
        LIS_rc%gridDesc(n,34) = lat(1,1)
        LIS_rc%gridDesc(n,35) = lon(1,1)
        LIS_rc%gridDesc(n,37) = lat(gnc(n),gnr(n))
        LIS_rc%gridDesc(n,38) = lon(gnc(n),gnr(n))
        
        deallocate(lat)
        deallocate(lon)

     else
        write(LIS_logunit,*) LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit,*) 'program stopping ...'
        call LIS_endrun
     endif
  enddo
#endif

  write(LIS_logunit,*)'DOMAIN details: Running in polar stereographic projection'
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
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     truelat = run_dd(n,3)
     stdlon = run_dd(n,4)
     orient = run_dd(n,5)
     xmesh = run_dd(n,6)
     nc = nint(run_dd(n,7))
     nr = nint(run_dd(n,8))

     call LIS_quilt_domain(n,nc,nr)

!call map_set on the overall running domain

     call map_set(PROJ_PS,stlat,stlon,&
          xmesh*1000.0,stdlon,truelat,&
          0.0,nc,nr,proj_temp)
     call ij_to_latlon(proj_temp,float(LIS_ews_halo_ind(n,LIS_localPet+1)),&
          float(LIS_nss_halo_ind(n,LIS_localPet+1)),&
          lat_str,lon_str)


     LIS_rc%gridDesc(n,1) = 5
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

!------------------------------------------------------------------------
! Read namelist of parameters depending on the domain
!------------------------------------------------------------------------
! for now set the global domain the same as local domain
! need to be modified to work on subsetted domains

     do k=1,13
        write(LIS_logunit,*)'(',k,',',LIS_rc%gridDesc(n,k),')'
     enddo

     LIS_rc%gnc(n) = nc
     LIS_rc%gnr(n) = nr
     LIS_rc%pnc(n) = LIS_rc%gridDesc(n,2)
     LIS_rc%pnr(n) = LIS_rc%gridDesc(n,3)

     call map_set(PROJ_PS,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
          LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),&
          0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)

!global domain
     call map_set(PROJ_LC, stlat, stlon, &
          LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),&
          LIS_rc%gridDesc(n,10),&
          0.0,LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_domain(n)%lisglbproj)

!parameter data 
     call map_set(PROJ_LC, LIS_rc%gridDesc(n,34),LIS_rc%gridDesc(n,35),&
          LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),&
          LIS_rc%gridDesc(n,10),&
          0.0,LIS_rc%pnc(n),LIS_rc%pnr(n),LIS_domain(n)%lisparamproj)

  enddo
  deallocate(run_dd)
  deallocate(gnc)
  deallocate(gnr)

#if 0 
  call ESMF_ConfigFindLabel(LIS_config,"tile coord file:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%tile_coord_file(i),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"tile veg file:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%tile_veg_file(i),rc=rc)
  enddo
#endif

end subroutine readinput_polar
      










