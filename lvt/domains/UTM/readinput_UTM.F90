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
! !ROUTINE: readinput_UTM
!  \label{readinput_UTM}
!
! !INTERFACE:
subroutine readinput_UTM
! 
! !USES:
  use ESMF 
  use LVT_coreMod
  use LVT_domainMod
  use LVT_logMod
  use map_utils
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
!  28 Jun 2009: Sujay Kumar; Initial Specification
! 
!EOP

  implicit none

  integer              :: i
  integer              :: k,n
  real, allocatable    :: run_dd(:)
  real, allocatable    :: input_run_dd(:)
!  real, allocatable    :: param_dd(:)
  integer              :: nc, nr
  integer              :: rc
  real(ESMF_KIND_R8)   :: stnorthing, steasting, dx

  allocate(run_dd(6))
  allocate(input_run_dd(6))
!  allocate(param_dd(6))

  call ESMF_ConfigFindLabel(LVT_config,"Run domain UTM zone:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(1),rc=rc)
  call LVT_verify(rc,'Run domain UTM zone: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain northing of SW corner:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(2),rc=rc)
  call LVT_verify(rc,'Run domain northing of SW corner: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain easting of SW corner:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(3),rc=rc)
  call LVT_verify(rc,'Run domain easting of SW corner: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain x-dimension size:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(4),rc=rc)
  call LVT_verify(rc,'Run domain x-dimension size: not defined')

  
  call ESMF_ConfigFindLabel(LVT_config,"Run domain y-dimension size:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(5),rc=rc)
  call LVT_verify(rc,'Run domain y-dimension size: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"Run domain resolution:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,run_dd(6),rc=rc)
  call LVT_verify(rc,'Run domain resolution: not defined')

  print*, 'TBD: readinput_UTM routine needs updating'
  stop
#if 0 
  call ESMF_ConfigFindLabel(LVT_config,"LIS run domain UTM zone:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,input_run_dd(1),rc=rc)
  call LVT_verify(rc,'LIS run domain UTM zone: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"LIS run domain northing of SW corner:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,input_run_dd(2),rc=rc)
  call LVT_verify(rc,'LIS run domain northing of SW corner: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"LIS run domain easting of SW corner:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,input_run_dd(3),rc=rc)
  call LVT_verify(rc,'LIS run domain easting of SW corner: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"LIS run domain x-dimension size:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,input_run_dd(4),rc=rc)
  call LVT_verify(rc,'LIS run domain x-dimension size: not defined')

  
  call ESMF_ConfigFindLabel(LVT_config,"LIS run domain y-dimension size:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,input_run_dd(5),rc=rc)
  call LVT_verify(rc,'LIS run domain y-dimension size: not defined')

  call ESMF_ConfigFindLabel(LVT_config,"LIS run domain resolution:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,input_run_dd(6),rc=rc)
  call LVT_verify(rc,'LIS run domain resolution: not defined')


  write(LVT_logunit,fmt=*)'DOMAIN details:'

  if(LVT_rc%npesx*LVT_rc%npesy.ne.LVT_npes) then 
     write(unit=LVT_logunit,fmt=*) 'Layout does not match the number of processors...'
     write(unit=LVT_logunit,fmt=*) 'npex, npey, ',LVT_rc%npesx,'x',LVT_rc%npesy,'!=',LVT_npes
     write(unit=LVT_logunit,fmt=*) 'Stopping program ..'
     call LVT_endrun()
  endif

  stnorthing = run_dd(2)
  steasting  = run_dd(3)
  
  dx = run_dd(6)
  nc = nint(run_dd(4))
  nr = nint(run_dd(5))
  
  call LVT_quilt_domain( nc,nr)
  
  LVT_rc%gridDesc(1) = 7
  
  LVT_rc%gridDesc(4) = stnorthing+ (LVT_nss_halo_ind(LVT_localPet+1)-1)*dx
  LVT_rc%gridDesc(5) = steasting + (LVT_ews_halo_ind(LVT_localPet+1)-1)*dx
  LVT_rc%gridDesc(7) = stnorthing+ (LVT_nse_halo_ind(LVT_localPet+1)-1)*dx
  LVT_rc%gridDesc(8) = steasting + (LVT_ewe_halo_ind(LVT_localPet+1)-1)*dx
  
  write(unit=LVT_logunit,fmt=*) 'local domain ',&
       LVT_rc%gridDesc(4),LVT_rc%gridDesc(7),&
       LVT_rc%gridDesc(5),LVT_rc%gridDesc(8), 'd: ',n
  
  LVT_rc%gnc = nint(run_dd(4))
  LVT_rc%gnr = nint(run_dd(5))
  
  LVT_rc%gridDesc(9) = run_dd(6)  !resolution
  
  LVT_rc%gridDesc(10) = run_dd(1) !zone
  LVT_rc%gridDesc(6) = 128
  LVT_rc%gridDesc(11) = 64
  LVT_rc%gridDesc(20) = 255
  
  if(LVT_rc%gridDesc(7).lt.LVT_rc%gridDesc(4)) then
     write(unit=LVT_logunit,fmt=*) 'northing2 must be greater than northing1'
     write(LVT_logunit,*) LVT_rc%gridDesc(7),LVT_rc%gridDesc(4)
     write(unit=LVT_logunit,fmt=*) 'Stopping run...'
     call LVT_endrun
  endif
  if(LVT_rc%gridDesc(8).lt.LVT_rc%gridDesc(5)) then
     write(unit=LVT_logunit,fmt=*) 'easting2 must be greater than easting1'
     write(LVT_logunit,*) LVT_rc%gridDesc(8),LVT_rc%gridDesc(5)
     write(unit=LVT_logunit,fmt=*) 'Stopping run...'
     call LVT_endrun
  endif
  
  LVT_rc%gridDesc(2) = nint((LVT_rc%gridDesc(8)-LVT_rc%gridDesc(5))&
       /dx)+ 1
  LVT_rc%gridDesc(3) = nint((LVT_rc%gridDesc(7)-LVT_rc%gridDesc(4))&
       /dx) + 1  
  write(LVT_logunit,*)'--------------------Domain ','----------------------'
  do k=1,13
     write(unit=LVT_logunit,fmt=*) '(',k,',',LVT_rc%gridDesc(k),')'
  enddo
  
  !     call map_set(PROJ_UTM, LVT_rc%gridDesc(4),LVT_rc%gridDesc(5),&
  !          0.0, LVT_rc%gridDesc(9),LVT_rc%gridDesc(10), 0.0,&
  !          LVT_rc%lnc,LVT_rc%lnr,LVT_domain%lisproj)
  
  ! parameter projection information
  
  stnorthing = input_run_dd(2)
  steasting  = input_run_dd(3)
  
  dx = input_run_dd(6)
  nc = nint(input_run_dd(4))
  nr = nint(input_run_dd(5))
  
  call LVT_quilt_lis_domain( nc,nr)
  
  LVT_rc%lis_gridDesc(1) = 7
  
  LVT_rc%lis_gridDesc(4) = stnorthing+ (LVT_lis_nss_halo_ind(LVT_localPet+1)-1)*dx
  LVT_rc%lis_gridDesc(5) = steasting + (LVT_lis_ews_halo_ind(LVT_localPet+1)-1)*dx
  LVT_rc%lis_gridDesc(7) = stnorthing+ (LVT_lis_nse_halo_ind(LVT_localPet+1)-1)*dx
  LVT_rc%lis_gridDesc(8) = steasting + (LVT_lis_ewe_halo_ind(LVT_localPet+1)-1)*dx

  write(unit=LVT_logunit,fmt=*) 'LIS domain ',&
       LVT_rc%lis_gridDesc(4),LVT_rc%lis_gridDesc(7),&
       LVT_rc%lis_gridDesc(5),LVT_rc%lis_gridDesc(8), 'd: ',n  

  LVT_rc%lis_gnc = nint(input_run_dd(4))
  LVT_rc%lis_gnr = nint(input_run_dd(5))
  
  LVT_rc%lis_gridDesc(9) = input_run_dd(6)  !resolution
  
  LVT_rc%lis_gridDesc(10) = input_run_dd(1) !zone
  LVT_rc%lis_gridDesc(6) = 128
  LVT_rc%lis_gridDesc(11) = 64
  LVT_rc%lis_gridDesc(20) = 255
  
  if(LVT_rc%lis_gridDesc(7).lt.LVT_rc%lis_gridDesc(4)) then
     write(unit=LVT_logunit,fmt=*) 'northing2 must be greater than northing1'
     write(LVT_logunit,*) LVT_rc%lis_gridDesc(7),LVT_rc%lis_gridDesc(4)
     write(unit=LVT_logunit,fmt=*) 'Stopping run...'
     call LVT_endrun
  endif
  if(LVT_rc%lis_gridDesc(8).lt.LVT_rc%lis_gridDesc(5)) then
     write(unit=LVT_logunit,fmt=*) 'easting2 must be greater than easting1'
     write(LVT_logunit,*) LVT_rc%lis_gridDesc(8),LVT_rc%lis_gridDesc(5)
     write(unit=LVT_logunit,fmt=*) 'Stopping run...'
     call LVT_endrun
  endif
  
  LVT_rc%lis_gridDesc(2) = nint((LVT_rc%lis_gridDesc(8)-LVT_rc%lis_gridDesc(5))&
       /dx)+ 1
  LVT_rc%lis_gridDesc(3) = nint((LVT_rc%lis_gridDesc(7)-LVT_rc%lis_gridDesc(4))&
       /dx) + 1  
  write(LVT_logunit,*)'--------------------Domain ','----------------------'
  do k=1,13
     write(unit=LVT_logunit,fmt=*) '(',k,',',LVT_rc%lis_gridDesc(k),')'
  enddo

!  do k=30,40
!     write(unit=LVT_logunit,fmt=*) '(',k,',',LVT_rc%lis_gridDesc(k),')'
!  enddo
  write(LVT_logunit,*)'--------------------------------------------------------------'
  LVT_rc%pnc = LVT_rc%lis_gridDesc(2)
  LVT_rc%pnr = LVT_rc%lis_gridDesc(3)
  

  deallocate(run_dd)
  deallocate(input_run_dd)
!  deallocate(param_dd)

  return
#endif
end subroutine readinput_UTM
      










