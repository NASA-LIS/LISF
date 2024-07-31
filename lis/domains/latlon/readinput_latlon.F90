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
! !ROUTINE: readinput_latlon
!  \label{readinput_latlon}
!  
! !REVISION HISTORY:
! 
!  11 Apr 2000: Brian Cosgrove; Added Elevation correction and Forcing
!               Mask read statements
!  14 Nov 2003: Sujay Kumar; Modified card file that includes regional 
!               modeling options
! !INTERFACE:
subroutine readinput_latlon
! !USES:
  use ESMF 
  use LIS_coreMod
  use LIS_domainMod
  use LIS_logMod
  use LIS_fileIOMod
  use map_utils
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LIS configuration 
!  file to determine the region of interest in the LIS simulation, in 
!  a lat/lon projection. The routine reads the extents of the 
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

  integer              :: i
  integer              :: k,n
  real, allocatable    :: run_dd(:,:)
  integer              :: nc, nr
  integer              :: rc
  real                 :: stlat, stlon, dx, dy
  integer              :: ftn
  integer              :: ios
  logical              :: file_exists
  integer              :: latId, lonId
  integer              :: ncId, nrId
  integer, allocatable :: gnc(:), gnr(:)
  real, allocatable    :: lat(:,:)
  real, allocatable    :: lon(:,:)

  allocate(run_dd(LIS_rc%nnest,6))
  allocate(gnc(LIS_rc%nnest))
  allocate(gnr(LIS_rc%nnest))

  call ESMF_ConfigFindLabel(LIS_config,"Run domain lower left lat:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,1),rc=rc)
     call LIS_verify(rc,'Run domain lower left lat: not defined')
  enddo  
  call ESMF_ConfigFindLabel(LIS_config,"Run domain lower left lon:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,2),rc=rc)
     call LIS_verify(rc,'Run domain lower left lon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain upper right lat:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,3),rc=rc)
     call LIS_verify(rc,'Run domain upper right lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Run domain upper right lon:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,4),rc=rc)
     call LIS_verify(rc,'Run domain upper right lon: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"Run domain resolution (dx):",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,5),rc=rc)
     call LIS_verify(rc,'Run domain resolution (dx): not defined')
  enddo  

  call ESMF_ConfigFindLabel(LIS_config,"Run domain resolution (dy):",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,6),rc=rc)
     call LIS_verify(rc,'Run domain resolution (dy): not defined')
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
        
        if(abs(run_dd(n,1) - lat(1,1)).gt.0.0001.or.&
             abs(run_dd(n,3)-lat(gnc(n),gnr(n))).gt.0.0001.or.&
             abs(run_dd(n,2)-lon(1,1)).gt.0.0001.or.&
             abs(run_dd(n,4)-lon(gnc(n),gnr(n))).gt.0.0001) then 
           write(LIS_logunit,*) '[ERR] The selected running domain is outside the '
           write(LIS_logunit,*) '[ERR] domain provided through the parameter attributes file '
           write(LIS_logunit,*) '[ERR] ',abs(run_dd(n,1) - lat(1,1))
           write(LIS_logunit,*) '[ERR] ',abs(run_dd(n,3) - lat(gnc,gnr))
           write(LIS_logunit,*) '[ERR] ',abs(run_dd(n,2) - lon(1,1))
           write(LIS_logunit,*) '[ERR] ',abs(run_dd(n,4) - lon(gnc,gnr))
           
           write(LIS_logunit,*) '[ERR] program stopping ...'
           call LIS_endrun
        endif
        
        LIS_rc%gridDesc(n,34) = lat(1,1)
        LIS_rc%gridDesc(n,35) = lon(1,1)
        LIS_rc%gridDesc(n,37) = lat(gnc(n),gnr(n))
        LIS_rc%gridDesc(n,38) = lon(gnc(n),gnr(n))
        
        deallocate(lat)
        deallocate(lon)

     else
        write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit,*) '[ERR] program stopping ...'
        call LIS_endrun
     endif
  enddo
#endif

  write(LIS_logunit,*)'[INFO] DOMAIN details:'

  if(LIS_rc%npesx*LIS_rc%npesy.ne.LIS_npes) then 
     write(LIS_logunit,*) '[ERR] Layout does not match the number of processors...'
     write(LIS_logunit,*) '[ERR] npex, npey, ',LIS_rc%npesx,'x',LIS_rc%npesy,'!=',LIS_npes
     write(LIS_logunit,*) '[ERR] Stopping program ..'
     call LIS_endrun()
  endif

  do n=1,LIS_rc%nnest
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     dx = run_dd(n,5)
     dy = run_dd(n,6)

     nc = (nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))) + 1
     nr = (nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))) + 1

     call LIS_quilt_domain(n, nc,nr)
     
     LIS_rc%gridDesc(n,4) = stlat + (LIS_nss_halo_ind(n,LIS_localPet+1)-1)*dy
     LIS_rc%gridDesc(n,5) = stlon + (LIS_ews_halo_ind(n,LIS_localPet+1)-1)*dx
     LIS_rc%gridDesc(n,7) = stlat + (LIS_nse_halo_ind(n,LIS_localPet+1)-1)*dy
     LIS_rc%gridDesc(n,8) = stlon + (LIS_ewe_halo_ind(n,LIS_localPet+1)-1)*dx
     
     write(LIS_logunit,*) '[INFO] local domain ',&
          LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,7),&
          LIS_rc%gridDesc(n,5),LIS_rc%gridDesc(n,8), 'd: ',n

     LIS_rc%gnc(n) = nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))+1
     LIS_rc%gnr(n) = nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))+1

     LIS_rc%gridDesc(n,1) = 0
     LIS_rc%gridDesc(n,9) = run_dd(n,5)

     if(LIS_rc%gridDesc(n,1).eq.0) then 
        LIS_rc%gridDesc(n,10) = run_dd(n,6)
        LIS_rc%gridDesc(n,6) = 128
        LIS_rc%gridDesc(n,11) = 64
        LIS_rc%gridDesc(n,20) = 64
     endif
     if(LIS_rc%gridDesc(n,7).lt.LIS_rc%gridDesc(n,4)) then
        write(LIS_logunit,*) '[ERR] lat2 must be greater than lat1'
        write(LIS_logunit,*) '[ERR] ',LIS_rc%gridDesc(n,7),LIS_rc%gridDesc(n,4)
        write(LIS_logunit,*) '[ERR] Stopping run...'
        call LIS_endrun
     endif
     if(LIS_rc%gridDesc(n,8).lt.LIS_rc%gridDesc(n,5)) then
        write(LIS_logunit,*) '[ERR] lon2 must be greater than lon1'
        write(LIS_logunit,*) '[ERR] ',LIS_rc%gridDesc(n,8),LIS_rc%gridDesc(n,5)
        write(LIS_logunit,*) '[ERR] Stopping run...'
        call LIS_endrun
     endif
     
     LIS_rc%gridDesc(n,2) = nint((LIS_rc%gridDesc(n,8)-LIS_rc%gridDesc(n,5))&
          /LIS_rc%gridDesc(n,9))+ 1
     LIS_rc%gridDesc(n,3) = nint((LIS_rc%gridDesc(n,7)-LIS_rc%gridDesc(n,4))&
          /LIS_rc%gridDesc(n,10)) + 1  
     write(LIS_logunit,*)'[INFO] --------------------Domain ',n,'----------------------'
     do k=1,13
        write(LIS_logunit,*) '[INFO] (',k,',',LIS_rc%gridDesc(n,k),')'
     enddo

!local domain of each processor
     call map_set(PROJ_LATLON, LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
          0.0, LIS_rc%gridDesc(n,9),LIS_rc%gridDesc(n,10), 0.0,&
          LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)

!global domain
     call map_set(PROJ_LATLON, stlat, stlon, &
          0.0, LIS_rc%gridDesc(n,9),LIS_rc%gridDesc(n,10), 0.0,&
          LIS_rc%gnc(n),LIS_rc%gnr(n),LIS_domain(n)%lisglbproj)

     LIS_rc%pnc(n) = gnc(n)
     LIS_rc%pnr(n) = gnr(n)

!parameter data 
     call map_set(PROJ_LATLON, LIS_rc%gridDesc(n,34),LIS_rc%gridDesc(n,35),&
          0.0, LIS_rc%gridDesc(n,9),LIS_rc%gridDesc(n,10), 0.0,&
          LIS_rc%pnc(n),LIS_rc%pnr(n),LIS_domain(n)%lisparamproj)

     write(LIS_logunit,*)'[INFO] --------------------------------------------------------------'

     LIS_rc%minLat(n) = LIS_rc%gridDesc(n,4) 
     LIS_rc%minLon(n) = LIS_rc%gridDesc(n,5) 
     LIS_rc%maxLat(n) = LIS_rc%gridDesc(n,7) 
     LIS_rc%maxLon(n) = LIS_rc%gridDesc(n,8) 
  enddo
  deallocate(run_dd)
  deallocate(gnc)
  deallocate(gnr)

  return

end subroutine readinput_latlon
      










