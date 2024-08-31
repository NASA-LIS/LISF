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
! !ROUTINE: readinput_lambert
!  \label{readinput_lambert}
!
! !INTERFACE:
subroutine readinput_lambert
! 
! !USES:
  use ESMF 
  use LVT_logMod
  use map_utils
  use LVT_coreMod
  use LVT_domainMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!  This routine reads the options specified in the LVT configuration 
!  file to determine the region of interest in the LVT simulation, in 
!  a lambert conformal projection. The routine reads the extents of the 
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest. The domain 
!  specifications of the parameter maps are also read in by this routine. 
!  Based on the number of processors and the processor layout specified, 
!  this routine also performs the domain decomposition. 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[lvttask\_for\_point](\ref{LVT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  15Jul2005: Sujay Kumar; Initial Specification
! 
!EOP

  integer           :: i
  integer           :: k
  real, allocatable :: run_dd(:)
  real, allocatable :: input_run_dd(:)
  real              :: lat_str,lon_str
  real              :: stlat,stlon,truelat1,truelat2,stdlon,xmesh
  integer           :: nc,nr
  type(proj_info)   :: proj_temp
  integer           :: n, rc, ios
  integer           :: latId, lonId
  integer           :: ncId, nrId
  integer           :: gnc, gnr
  integer           :: ftn
  character*50      :: map_proj
  logical           :: file_exists

  LVT_rc%gridDesc = 0 
  LVT_rc%input_gridDesc = 0 

  allocate(run_dd(8))
  allocate(input_run_dd(8))

  call ESMF_ConfigFindLabel(LVT_config,"Run domain lower left lat:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(1),rc=rc)
  call LVT_verify(rc,'Run domain lower left lat: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain lower left lon:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(2),rc=rc)
  call LVT_verify(rc,'Run domain lower left lon: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain true lat1:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(3),rc=rc)
  call LVT_verify(rc,'Run domain true lat1: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain true lat2:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(4),rc=rc)
  call LVT_verify(rc,'Run domain true lat2: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain standard lon:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(5),rc=rc)
  call LVT_verify(rc,'Run domain standard lon: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain resolution:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(6),rc=rc)
  call LVT_verify(rc,'Run domain resolution: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain x-dimension size:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(7),rc=rc)
  call LVT_verify(rc,'Run domain x-dimension size: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain y-dimension size:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(8),rc=rc)
  call LVT_verify(rc,'Run domain y-dimension size: not defined')


#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=LVT_rc%paramfile, exist=file_exists)
  if(file_exists) then 
     
     ios = nf90_open(path=LVT_rc%paramfile,&
          mode=NF90_NOWRITE,ncid=ftn)
     call LVT_verify(ios,'Error in nf90_open in readinput_latlon')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',map_proj)
     call LVT_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')
     
     if(map_proj.ne."LAMBERT CONFORMAL") then 
        write(LVT_logunit,*) 'The map projection in the domain and parameter'
        write(LVT_logunit,*) 'file is different from the map projection '
        write(LVT_logunit,*) 'intended for the LVT run (lambert)'
        write(LVT_logunit,*) 'Program stopping ....'
        call LVT_endrun()
     endif

     ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
          input_run_dd(1))
     call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
          input_run_dd(2))
     call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')

     ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',&
          input_run_dd(3))
     call LVT_verify(ios, 'Error in nf90_get_att: TRUELAT1')

     ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT2',&
          input_run_dd(4))
     call LVT_verify(ios, 'Error in nf90_get_att: TRUELAT2')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',&
          input_run_dd(5))
     call LVT_verify(ios, 'Error in nf90_get_att: STANDARD_LON')

     ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',input_run_dd(6))
     call LVT_verify(ios, 'Error in nf90_get_att: DX')
     
     ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',input_run_dd(6))
     call LVT_verify(ios, 'Error in nf90_get_att: DY')
     
     ios = nf90_inq_dimid(ftn,"east_west",ncId)
     call LVT_verify(ios,'Error in nf90_inq_dimid ')
     
     ios = nf90_inq_dimid(ftn,"north_south",nrId)
     call LVT_verify(ios,'Error in nf90_inq_dimid ')
     
     ios = nf90_inquire_dimension(ftn,ncId, len=gnc)
     call LVT_verify(ios,'Error in nf90_inquire_dimension')
     input_run_dd(7) = gnc

     ios = nf90_inquire_dimension(ftn,nrId, len=gnr)
     call LVT_verify(ios,'Error in nf90_inquire_dimension ')
     input_run_dd(8) = gnr

     ios = nf90_close(ftn)
     call LVT_verify(ios,'Error in nf90_close')
  else
     write(LVT_logunit,*) '[ERR] ',LVT_rc%paramfile, ' does not exist'
     call LVT_endrun
  endif
#endif

  write(LVT_logunit,*)'DOMAIN details: Running in lambert projection'
  if(LVT_npes.eq.1) then 
     LVT_rc%npesx = 1
     LVT_rc%npesy = 1
  endif

  if(LVT_rc%npesx*LVT_rc%npesy.ne.LVT_npes) then 
     write(LVT_logunit,*) 'Layout does not match the number of processors...'
     write(LVT_logunit,*) 'npex, npey, ',LVT_rc%npesx,'x',LVT_rc%npesy,'!=',LVT_npes
     write(LVT_logunit,*) 'Stopping program ..'
     call LVT_endrun()
  endif

  stlat = run_dd(1)
  stlon = run_dd(2)
  truelat1 = run_dd(3)
  truelat2 = run_dd(4)
  stdlon = run_dd(5)
  xmesh = run_dd(6)
  nc = nint(run_dd(7))
  nr = nint(run_dd(8))
  
  call LVT_quilt_domain(nc,nr)
  
!call map_set on the overall running domain

  call map_set(PROJ_LC,stlat,stlon,&
       xmesh*1000.0,stdlon,truelat1,&
       truelat2,nc,nr,proj_temp)
  call ij_to_latlon(proj_temp,float(LVT_ews_halo_ind(LVT_localPet+1)),&
       float(LVT_nss_halo_ind(LVT_localPet+1)),lat_str,lon_str)     
  
  LVT_rc%gridDesc(1) = 3
  LVT_rc%gridDesc(2) = LVT_rc%lnc
  LVT_rc%gridDesc(3) = LVT_rc%lnr
  LVT_rc%gridDesc(4) = lat_str
  LVT_rc%gridDesc(5) = lon_str
  LVT_rc%gridDesc(6) = 8
  LVT_rc%gridDesc(7) = truelat2
  LVT_rc%gridDesc(8) = xmesh
  LVT_rc%gridDesc(9) = xmesh
  LVT_rc%gridDesc(10) = truelat1
  LVT_rc%gridDesc(11) = stdlon
  
!------------------------------------------------------------------------
! Read namelvtt of parameters depending on the domain
!------------------------------------------------------------------------
! for now set the global domain the same as local domain
! need to be modified to work on subsetted domains

  do k=1,13
     write(LVT_logunit,*)'(',k,',',LVT_rc%gridDesc(k),')'
  enddo

  LVT_rc%gnc = nc
  LVT_rc%gnr = nr
  
  call map_set(PROJ_LC,LVT_rc%gridDesc(4),LVT_rc%gridDesc(5),&
       LVT_rc%gridDesc(8)*1000.0,LVT_rc%gridDesc(11),&
       LVT_rc%gridDesc(10),&
       truelat2,LVT_rc%lnc,LVT_rc%lnr,LVT_domain%lvtproj)

  call map_set(PROJ_LC,LVT_rc%gridDesc(4),LVT_rc%gridDesc(5),&
       LVT_rc%gridDesc(8)*1000.0,LVT_rc%gridDesc(11),&
       LVT_rc%gridDesc(10),&
       truelat2,LVT_rc%lnc,LVT_rc%lnr,LVT_domain%lvtparamproj)

  stlat = input_run_dd(1)
  stlon = input_run_dd(2)
  truelat1 = input_run_dd(3)
  truelat2 = input_run_dd(4)
  stdlon = input_run_dd(5)
  xmesh = input_run_dd(6)
  nc = nint(input_run_dd(7))
  nr = nint(input_run_dd(8))
  
!call map_set on the overall running domain

  call map_set(PROJ_LC,stlat,stlon,&
       xmesh*1000.0,stdlon,truelat1,&
       truelat2,nc,nr,proj_temp)
  call ij_to_latlon(proj_temp,1.0,1.0,lat_str,lon_str)
  
  LVT_rc%input_gridDesc(1) = 3
  LVT_rc%input_gridDesc(2) = nc
  LVT_rc%input_gridDesc(3) = nr
  LVT_rc%input_gridDesc(4) = lat_str
  LVT_rc%input_gridDesc(5) = lon_str
  LVT_rc%input_gridDesc(6) = 8
  LVT_rc%input_gridDesc(7) = truelat2
  LVT_rc%input_gridDesc(8) = xmesh
  LVT_rc%input_gridDesc(9) = xmesh
  LVT_rc%input_gridDesc(10) = truelat1
  LVT_rc%input_gridDesc(11) = stdlon
  
!------------------------------------------------------------------------
! Read namelvtt of parameters depending on the domain
!------------------------------------------------------------------------
! for now set the global domain the same as local domain
! need to be modified to work on subsetted domains

  write(LVT_logunit,*)'--------------------LIS Domain ','----------------------'     
  do k=1,13
     write(LVT_logunit,*)'(',k,',',LVT_rc%input_gridDesc(k),')'
  enddo
  
!  do k=30,40
!     write(LVT_logunit,*)'(',k,',',LVT_rc%input_gridDesc(k),')'
!  enddo
  
  LVT_rc%input_lnc = nc
  LVT_rc%input_lnr = nr
  LVT_rc%pnc = LVT_rc%input_gridDesc(2)
  LVT_rc%pnr = LVT_rc%input_gridDesc(3)
  
!  call map_set(PROJ_LC,LVT_rc%input_gridDesc(4),LVT_rc%input_gridDesc(5),&
!       LVT_rc%input_gridDesc(8)*1000.0,LVT_rc%input_gridDesc(11),&
!       LVT_rc%input_gridDesc(10),&
!       truelat2,LVT_rc%input_lnc,LVT_rc%input_lnr,LVT_domain%lisproj)
  
  deallocate(input_run_dd)
  deallocate(run_dd)
!  deallocate(param_dd)

end subroutine readinput_lambert
      










