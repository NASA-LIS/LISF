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
! !ROUTINE: readinput_gaussian
!  \label{readinput_gaussian}
!  
! !REVISION HISTORY:
! 
!  13 Nov 2006: James Geiger; Initial specification
!
! !INTERFACE:
subroutine readinput_gaussian
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config, LIS_deltas, LIS_offsets, &
       LIS_localPet, LIS_masterproc, LIS_npes, LIS_ews_halo_ind, &
       LIS_ewe_halo_ind,LIS_nss_halo_ind, LIS_nse_halo_ind
  use LIS_logMod, only :LIS_logunit, LIS_endrun
  use LIS_domainMod, only : LIS_quilt_domain, LIS_setParameterDomainSpecs
  use LIS_fileIOMod, only : LIS_readDomainConfigSpecs
  use gaussian_mod

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LIS configuration 
!  file to determine the region of interest in the LIS simulation, in 
!  a quasi-regular Gaussian lat/lon projection. 
!  The routine reads the extents of the 
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
  implicit none

  real, allocatable :: run_dd(:,:)
  real, allocatable :: param_dd(:,:)
  integer :: nc,nr,cindex,rindex
  integer :: i,j,k,n,rc,imax,jmax
  real              :: stlat, stlon, dx
  real,allocatable,dimension(:) :: lat_str,lat_end,lon_str,lon_end

  if(LIS_rc%halox.gt.0.and.LIS_rc%haloy.gt.0) then 
     write(LIS_logunit,*) 'Halo support for gaussian domain needs testing'
     write(LIS_logunit,*) 'Program stopping '
     call LIS_endrun
  endif

  allocate(run_dd(LIS_rc%nnest,6))
  allocate(param_dd(LIS_rc%nnest,6))
  call ESMF_ConfigFindLabel(LIS_config,"Run domain first grid point lat:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,1),rc=rc)
  enddo  
  call ESMF_ConfigFindLabel(LIS_config,"Run domain first grid point lon:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,2),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain last grid point lat:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,3),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain last grid point lon:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,4),rc=rc)
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"Run domain resolution dlon:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,5),rc=rc)
  enddo  

  call ESMF_ConfigFindLabel(LIS_config,"Run domain number of lat circles:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,6),rc=rc)
  enddo  

  call LIS_readDomainConfigSpecs("param domain", param_dd)

  write(LIS_logunit,fmt=*)'DOMAIN details:'
  if(LIS_npes.eq.1) then 
     LIS_rc%npesx = 1
     LIS_rc%npesy = 1
  endif

  if(LIS_rc%npesx*LIS_rc%npesy.ne.LIS_npes) then 
     write(unit=LIS_logunit,fmt=*) 'Layout does not match the number of processors...'
     write(unit=LIS_logunit,fmt=*) 'npex, npey, ',LIS_rc%npesx,'x',LIS_rc%npesy,'!=',LIS_npes
     write(unit=LIS_logunit,fmt=*) 'Stopping program ..'
     call LIS_endrun()
  endif

  allocate(lat_str(LIS_rc%npesy))
  allocate(lon_str(LIS_rc%npesx))
  allocate(lat_end(LIS_rc%npesy))
  allocate(lon_end(LIS_rc%npesx))

  do n=1,LIS_rc%nnest

     dx    = run_dd(n,5)
     jmax  = 2*run_dd(n,6) ! total number of latitude circles
     imax = nint(360/dx)

     call gaussian_comp_lats(jmax)
     call gaussian_comp_lons(param_dd(n,2), param_dd(n,5))
     !call gaussian_read_grid(imax,jmax)

     nc = gaussian_find_col(run_dd(n,4)) - &
          gaussian_find_col(run_dd(n,2)) + 1
     nr = gaussian_find_row(run_dd(n,3)) - &
          gaussian_find_row(run_dd(n,1)) + 1

     stlat = run_dd(n,1)
     stlon = run_dd(n,2)

     LIS_rc%gridDesc(n,1)  = 4
     !LIS_rc%gridDesc(n,4)  = run_dd(n,1)
     !LIS_rc%gridDesc(n,5)  = run_dd(n,2)
     !LIS_rc%gridDesc(n,7)  = run_dd(n,3)
     !LIS_rc%gridDesc(n,8)  = run_dd(n,4)
     LIS_rc%gridDesc(n,9)  = run_dd(n,5)
     LIS_rc%gridDesc(n,10) = run_dd(n,6)

     !LIS_rc%gridDesc(n,2) = nint(360/LIS_rc%gridDesc(n,9))
     !LIS_rc%gridDesc(n,3) = 2*LIS_rc%gridDesc(n,10)

     LIS_rc%gridDesc(n,6)  = 128
     LIS_rc%gridDesc(n,11) = 64
     LIS_rc%gridDesc(n,20) = 64

     call LIS_setParameterDomainSpecs(param_dd(n,:),LIS_rc%gridDesc(n,:))

     LIS_rc%gnc(n) = nc
     LIS_rc%gnr(n) = nr

     call LIS_quilt_domain(n,nc,nr)

     LIS_rc%gridDesc(n,2) = nc
     LIS_rc%gridDesc(n,3) = nr
     
     rindex = gaussian_find_row(stlat) + LIS_nss_halo_ind(n,LIS_localPet+1) - 1
     cindex = gaussian_find_col(stlon) + LIS_ews_halo_ind(n,LIS_localPet+1) - 1
     !LIS_rc%gridDesc(n,4) = gaussian_latitudes(cindex,rindex)
     !LIS_rc%gridDesc(n,5) = gaussian_longitudes(cindex,rindex)
     LIS_rc%gridDesc(n,4) = gaussian_lat_array(rindex)
     LIS_rc%gridDesc(n,5) = gaussian_lon_array(cindex)

     rindex = gaussian_find_row(stlat) + LIS_nse_halo_ind(n,LIS_localPet+1) - 1
     cindex = gaussian_find_col(stlon) + LIS_ewe_halo_ind(n,LIS_localPet+1) - 1
     !LIS_rc%gridDesc(n,7) = gaussian_latitudes(cindex,rindex)
     !LIS_rc%gridDesc(n,8) = gaussian_longitudes(cindex,rindex)
     LIS_rc%gridDesc(n,7) = gaussian_lat_array(rindex)
     LIS_rc%gridDesc(n,8) = gaussian_lon_array(cindex)
     
     write(unit=LIS_logunit,fmt=*) 'local domain ',LIS_localPet, ',', &
                               LIS_rc%gridDesc(n,4),        &
                               LIS_rc%gridDesc(n,7),        &
                               LIS_rc%gridDesc(n,5),        &
                               LIS_rc%gridDesc(n,8),        &
                               'd: ',n

  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'--------------------Domain ',n,'----------------------'
     do k=1,20
        write(unit=LIS_logunit,fmt=*) '(',k,',',LIS_rc%gridDesc(n,k),')'
     enddo
     do k=41,50
        write(unit=LIS_logunit,fmt=*) '(',k,',',LIS_rc%gridDesc(n,k),')'
     enddo
     write(LIS_logunit,*)'--------------------------------------------------------------'

     LIS_rc%pnc(n) = LIS_rc%gridDesc(n,42)
     LIS_rc%pnr(n) = LIS_rc%gridDesc(n,43)
  enddo

  deallocate(lat_str)
  deallocate(lon_str)
  deallocate(lat_end)
  deallocate(lon_end)
  
  deallocate(run_dd)
  deallocate(param_dd) 

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

  return

end subroutine readinput_gaussian

