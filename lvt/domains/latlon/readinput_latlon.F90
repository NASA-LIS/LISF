!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readinput_latlon
!  \label{readinput_latlon}
!
! !INTERFACE:
subroutine readinput_latlon
! 
! !USES:
  use ESMF 
  use LVT_coreMod
  use LVT_domainMod
  use LVT_logMod
  use map_utils
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
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
!  \item[listask\_for\_point](\ref{LVT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY:
! 
!  14 Nov 2003: Sujay Kumar; Modified card file that includes regional 
!               modeling options
! 
!EOP

  implicit none

  integer              :: i, j,count
  integer              :: k,n
  real, allocatable    :: run_dd(:)
  real, allocatable    :: input_run_dd(:)
!  real, allocatable    :: param_dd(:)
  integer              :: lnc,lnr
  integer              :: nc, nr
  integer              :: ierr
  integer              :: rc
  integer              :: ips, ipe, jps, jpe
  integer              :: Px, Py, P
  integer              :: mytask_x, mytask_y
  real(ESMF_KIND_R8)   :: stlat, stlon, dx, dy
  integer              :: status
  integer              :: ftn
  integer              :: ios
  logical              :: file_exists
  integer              :: latId, lonId
  integer              :: ncId, nrId
  integer              :: gnc, gnr
  character*50         :: map_proj
  real, allocatable        :: lat(:,:)
  real, allocatable        :: lon(:,:)


  LVT_rc%gridDesc = 0 
  LVT_rc%input_gridDesc = 0 

  allocate(run_dd(6))
  allocate(input_run_dd(6))
!  allocate(param_dd(6))

  call ESMF_ConfigFindLabel(LVT_config,"Run domain lower left lat:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(1),rc=rc)
  call LVT_verify(rc, 'Run domain lower left lat: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain lower left lon:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(2),rc=rc)
  call LVT_verify(rc, 'Run domain lower left lon: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain upper right lat:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(3),rc=rc)
  call LVT_verify(rc, 'Run domain upper right lat: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain upper right lon:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(4),rc=rc)
  call LVT_verify(rc, 'Run domain upper right lon: not defined')
  
  call ESMF_ConfigFindLabel(LVT_config,"Run domain resolution (dx):",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(5),rc=rc)
  call LVT_verify(rc, 'Run domain resolution (dx): not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain resolution (dy):",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(6),rc=rc)
  call LVT_verify(rc, 'Run domain resolution (dy): not defined')

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=LVT_rc%paramfile, exist=file_exists)
  if(file_exists) then 
     
     ios = nf90_open(path=LVT_rc%paramfile,&
          mode=NF90_NOWRITE,ncid=ftn)
     call LVT_verify(ios,'Error in nf90_open in readinput_latlon')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',map_proj)
     call LVT_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')
     
     if(map_proj.ne."EQUIDISTANT CYLINDRICAL") then 
        write(LVT_logunit,*) '[ERR] The map projection in the domain and parameter'
        write(LVT_logunit,*) '[ERR] file is different from the map projection '
        write(LVT_logunit,*) '[ERR] intended for the LVT run (latlon)'
        call LVT_endrun()
     endif

     ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
          input_run_dd(1))
     call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
          input_run_dd(2))
     call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',input_run_dd(5))
     call LVT_verify(ios, 'Error in nf90_get_att: DX')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',input_run_dd(6))
     call LVT_verify(ios, 'Error in nf90_get_att: DY')
     
     ios = nf90_inq_dimid(ftn,"east_west",ncId)
     call LVT_verify(ios,'Error in nf90_inq_dimid ')
     
     ios = nf90_inq_dimid(ftn,"north_south",nrId)
     call LVT_verify(ios,'Error in nf90_inq_dimid ')
     
     ios = nf90_inquire_dimension(ftn,ncId, len=gnc)
     call LVT_verify(ios,'Error in nf90_inquire_dimension')
     
     ios = nf90_inquire_dimension(ftn,nrId, len=gnr)
     call LVT_verify(ios,'Error in nf90_inquire_dimension ')
     
     input_run_dd(3) = input_run_dd(1) + (gnr-1)*input_run_dd(6)
     input_run_dd(4) = input_run_dd(2) + (gnc-1)*input_run_dd(5)
     

     ios = nf90_close(ftn)
     call LVT_verify(ios,'Error in nf90_close')

     LVT_rc%pnc = gnc
     LVT_rc%pnr = gnr 

     
  else
     write(LVT_logunit,*) '[ERR] ',LVT_rc%paramfile, ' does not exist'
     call LVT_endrun
  endif
#endif

  stlat = run_dd(1)
  stlon = run_dd(2)
  dx = run_dd(5)
  dy = run_dd(6)
  
  nc = (nint((run_dd(4)-run_dd(2))/run_dd(5))) + 1
  nr = (nint((run_dd(3)-run_dd(1))/run_dd(6))) + 1
    
  call LVT_quilt_domain(nc,nr)

  LVT_rc%gridDesc(4) = stlat + (LVT_nss_halo_ind(LVT_localPet+1)-1)*dy
  LVT_rc%gridDesc(5) = stlon + (LVT_ews_halo_ind(LVT_localPet+1)-1)*dx
  LVT_rc%gridDesc(7) = stlat + (LVT_nse_halo_ind(LVT_localPet+1)-1)*dy
  LVT_rc%gridDesc(8) = stlon + (LVT_ewe_halo_ind(LVT_localPet+1)-1)*dx
  
  write(unit=LVT_logunit,fmt=*) '[INFO] Local domain ',&
       LVT_rc%gridDesc(4),LVT_rc%gridDesc(7),&
       LVT_rc%gridDesc(5),LVT_rc%gridDesc(8)
  
  LVT_rc%gnc = nint((run_dd(4)-run_dd(2))/run_dd(5))+1
  LVT_rc%gnr = nint((run_dd(3)-run_dd(1))/run_dd(6))+1

  LVT_rc%gridDesc(1) = 0
  LVT_rc%gridDesc(9) = run_dd(5)
  
  if(LVT_rc%gridDesc(1).eq.0) then 
     LVT_rc%gridDesc(10) = run_dd(6)
     LVT_rc%gridDesc(6) = 128
     LVT_rc%gridDesc(11) = 64
     LVT_rc%gridDesc(20) = 64
  endif
  if(LVT_rc%gridDesc(7).lt.LVT_rc%gridDesc(4)) then
     write(unit=LVT_logunit,fmt=*) '[ERR] lat2 must be greater than lat1'
     write(LVT_logunit,*) '[ERR] ', LVT_rc%gridDesc(7),LVT_rc%gridDesc(4)
     call LVT_endrun
  endif
  if(LVT_rc%gridDesc(8).lt.LVT_rc%gridDesc(5)) then
     write(unit=LVT_logunit,fmt=*) '[ERR] lon2 must be greater than lon1'
     write(LVT_logunit,*) '[ERR] ',LVT_rc%gridDesc(8),LVT_rc%gridDesc(5)
     call LVT_endrun
  endif
  
  LVT_rc%gridDesc(2) = nint((LVT_rc%gridDesc(8)-LVT_rc%gridDesc(5))&
       /LVT_rc%gridDesc(9))+ 1
  LVT_rc%gridDesc(3) = nint((LVT_rc%gridDesc(7)-LVT_rc%gridDesc(4))&
       /LVT_rc%gridDesc(10)) + 1  
  write(LVT_logunit,*)'[INFO] --------------------Domain ','----------------------'
  do k=1,13
     write(unit=LVT_logunit,fmt=*) '[INFO] (',k,',',LVT_rc%gridDesc(k),')'
  enddo

  call map_set(PROJ_LATLON, LVT_rc%gridDesc(4),LVT_rc%gridDesc(5),&
       0.0, LVT_rc%gridDesc(9),LVT_rc%gridDesc(10), 0.0,&
       LVT_rc%lnc,LVT_rc%lnr,LVT_domain%lvtproj)
  
  stlat = input_run_dd(1)
  stlon = input_run_dd(2)
  dx = input_run_dd(5)
  dy = input_run_dd(6)
  
  nc = (nint((input_run_dd(4)-input_run_dd(2))/input_run_dd(5))) + 1
  nr = (nint((input_run_dd(3)-input_run_dd(1))/input_run_dd(6))) + 1

  LVT_rc%input_gridDesc(4) = stlat 
  LVT_rc%input_gridDesc(5) = stlon 
  LVT_rc%input_gridDesc(7) = stlat + (nr-1)*dy
  LVT_rc%input_gridDesc(8) = stlon + (nc-1)*dx
  
  write(unit=LVT_logunit,fmt=*) '[INFO] Input domain ',&
       LVT_rc%input_gridDesc(4),LVT_rc%input_gridDesc(7),&
       LVT_rc%input_gridDesc(5),LVT_rc%input_gridDesc(8)
  
  LVT_rc%input_lnc = nint((input_run_dd(4)-input_run_dd(2))/input_run_dd(5))+1
  LVT_rc%input_lnr = nint((input_run_dd(3)-input_run_dd(1))/input_run_dd(6))+1
  LVT_rc%input_gridDesc(1) = 0
  LVT_rc%input_gridDesc(9) = input_run_dd(5)
  
  if(LVT_rc%input_gridDesc(1).eq.0) then 
     LVT_rc%input_gridDesc(10) = input_run_dd(6)
     LVT_rc%input_gridDesc(6) = 128
     LVT_rc%input_gridDesc(11) = 64
     LVT_rc%input_gridDesc(20) = 64
  endif
  if(LVT_rc%input_gridDesc(7).lt.LVT_rc%input_gridDesc(4)) then
     write(unit=LVT_logunit,fmt=*) '[ERR] lat2 must be greater than lat1'
     write(LVT_logunit,*) '[ERR] ',LVT_rc%input_gridDesc(7),LVT_rc%input_gridDesc(4)
     call LVT_endrun
  endif
  if(LVT_rc%input_gridDesc(8).lt.LVT_rc%input_gridDesc(5)) then
     write(unit=LVT_logunit,fmt=*) '[ERR] lon2 must be greater than lon1'
     write(LVT_logunit,*) '[ERR] ',LVT_rc%input_gridDesc(8),LVT_rc%input_gridDesc(5)
     call LVT_endrun
  endif
  
  call map_set(PROJ_LATLON, LVT_rc%input_gridDesc(34), LVT_rc%input_gridDesc(35),&
       0.0, LVT_rc%input_gridDesc(9),LVT_rc%input_gridDesc(10), 0.0,&
       LVT_rc%pnc,LVT_rc%pnr,LVT_domain%lvtparamproj)

  LVT_rc%input_gridDesc(2) = nint((LVT_rc%input_gridDesc(8)-LVT_rc%input_gridDesc(5))&
       /LVT_rc%input_gridDesc(9))+ 1
  LVT_rc%input_gridDesc(3) = nint((LVT_rc%input_gridDesc(7)-LVT_rc%input_gridDesc(4))&
       /LVT_rc%input_gridDesc(10)) + 1  

  LVT_rc%input_lnc = nint(LVT_rc%input_gridDesc(2))
  LVT_rc%input_lnr = nint(LVT_rc%input_gridDesc(3))

  write(LVT_logunit,*)'[INFO] --------------------Input Domain ','----------------------'
  do k=1,13
     write(unit=LVT_logunit,fmt=*) '[INFO] (',k,',',LVT_rc%input_gridDesc(k),')'
  enddo
!  call map_set(PROJ_LATLON, LVT_rc%input_gridDesc(4),LVT_rc%input_gridDesc(5),&
!       0.0, LVT_rc%input_gridDesc(9), LVT_rc%input_gridDesc(10),0.0,&
!       LVT_rc%input_lnc, LVT_rc%input_lnr,LVT_domain%lisproj)

!  do k=30,40
!     write(unit=LVT_logunit,fmt=*) '(',k,',',LVT_rc%input_gridDesc(k),')'
!  enddo
  write(LVT_logunit,*)'[INFO] --------------------------------------------------------------'

  deallocate(input_run_dd)
  deallocate(run_dd)
!  deallocate(param_dd)
end subroutine readinput_latlon
      










