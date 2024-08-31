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
subroutine readinput_latlon(nest)

! !USES:
  use ESMF 
  use LDT_coreMod
  use LDT_domainMod
  use LDT_logMod
  use LDT_fileIOMod
  use map_utils

  implicit none

  integer, intent(in) :: nest

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LIS configuration 
!  file to determine the region of interest in the LIS simulation, in 
!  a lat/lon projection. The routine reads the extents of the 
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest. 
!  Based on the number of processors and the processor layout specified, 
!  this routine also performs the domain decomposition. 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[diff\_lon] (\ref{diff_lon}) \newline
!    funtion to compute the difference between two given longitude values
!  \item[ldttask\_for\_point](\ref{LDT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
!EOP

  integer              :: i, j, k, n,count
  real, allocatable    :: run_dd(:,:)  ! LIS Target grid/domain
  integer              :: lnc,lnr
  integer              :: nc, nr
  real                 :: diff_lon
  integer              :: ierr, rc
  integer              :: ips, ipe, jps, jpe
  integer              :: Px, Py, P
  integer              :: mytask_x, mytask_y
  real(ESMF_KIND_R8)   :: stlat, stlon, dx, dy
  integer              :: status
  integer              :: nDupl
  logical              :: checkDuplicate
  logical              :: newFormatCheck 
! __________________________________________________________________________

  allocate(run_dd(LDT_rc%nnest,20))

  LDT_rc%lis_map_resfactor(nest) = 1.

  newFormatCheck = .false. 
! The following is not a full proof test of duplicates
  checkDuplicate = .false. 
  nDupl = 1
  do n=2,LDT_rc%nnest
     if(LDT_rc%lis_map_proj(n).eq.LDT_rc%lis_map_proj(1)) then 
        checkDuplicate = .true. 
        nDupl = nDupl + 1
     endif
  enddo
  nDupl = min(nest,nDupl)

  call ESMF_ConfigFindLabel(LDT_config,"LL run domain lower left lat:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,1),rc=rc)
     if(rc.eq.0) newFormatCheck = .true.
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"LL run domain lower left lon:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,2),rc=rc)
     if(rc.eq.0) newFormatCheck = .true.
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"LL run domain upper right lat:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,3),rc=rc)
     if(rc.eq.0) newFormatCheck = .true.
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"LL run domain upper right lon:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,4),rc=rc)
     if(rc.eq.0) newFormatCheck = .true.
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"LL run domain resolution (dx):",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,5),rc=rc)
     if(rc.eq.0) newFormatCheck = .true.
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"LL run domain resolution (dy):",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,6),rc=rc)
     if(rc.eq.0) newFormatCheck = .true.
  enddo

  if(.not.newFormatCheck) then 
     call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lat:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,1),rc=rc)
        call LDT_verify(rc, 'Run domain lower left lat: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lon:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,2),rc=rc)
        call LDT_verify(rc, 'Run domain lower left lon: not defined')
     enddo
     call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lat:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,3),rc=rc)
        call LDT_verify(rc, 'Run domain upper right lat: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lon:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,4),rc=rc)
        call LDT_verify(rc, 'Run domain upper right lon: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution (dx):",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,5),rc=rc)
        call LDT_verify(rc, 'Run domain resolution (dx): not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution (dy):",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,6),rc=rc)
        call LDT_verify(rc, 'Run domain resolution (dy): not defined')
     enddo
  endif

  do n=1,nDupl
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     dx    = run_dd(n,5)
     dy    = run_dd(n,6)

     ! Check if lon values are reverse (e.g., crosses International Dateline)
     if( run_dd(n,4)-run_dd(n,2) < 0. ) then   ! Check long. extents
        nc = nint(diff_lon(run_dd(n,4),run_dd(n,2))/run_dd(n,5)) + 1
     else
        nc = (nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))) + 1
     endif
     nr = (nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))) + 1
     LDT_rc%gnc(nest) = nc
     LDT_rc%gnr(nest) = nr
     
     ! Quilt domain - decompose global domain:
     call LDT_quilt_domain(nest,nc,nr)
  
     ! Assign LIS domain gridDesc elements, based on decomposed subdomains:
     LDT_rc%gridDesc(nest,4) = stlat + (LDT_nss_halo_ind(nest,LDT_localPet+1)-1)*dy
     LDT_rc%gridDesc(nest,5) = stlon + (LDT_ews_halo_ind(nest,LDT_localPet+1)-1)*dx
     LDT_rc%gridDesc(nest,7) = stlat + (LDT_nse_halo_ind(nest,LDT_localPet+1)-1)*dy

     LDT_rc%gridDesc(nest,8) = stlon + (LDT_ewe_halo_ind(nest,LDT_localPet+1)-1)*dx
     if( LDT_rc%gridDesc(nest,8) > 180. ) then
       LDT_rc%gridDesc(nest,8) = LDT_rc%gridDesc(nest,8) - 360.0 
     endif 

     write(unit=LDT_logunit,fmt=*) 'local domain: ',&
          LDT_rc%gridDesc(nest,4),LDT_rc%gridDesc(nest,7),&
          LDT_rc%gridDesc(nest,5),LDT_rc%gridDesc(nest,8)
     
     LDT_rc%gridDesc(nest,1) = 0                ! latlon grid
     LDT_rc%gridDesc(nest,9) = run_dd(n,5)      ! dx
     if(LDT_rc%gridDesc(nest,1).eq.0) then 
        LDT_rc%gridDesc(nest,10) = run_dd(n,6)  ! dy
        LDT_rc%gridDesc(nest,6)  = 128
        LDT_rc%gridDesc(nest,11) = 64
        LDT_rc%gridDesc(nest,20) = 64
     endif

     ! Original checks for when LAT2<LAT1; LON2<LON1:

     if(LDT_rc%gridDesc(nest,7).lt.LDT_rc%gridDesc(nest,4)) then
        write(LDT_logunit,*) '[ERR] lat2 must be greater than lat1 ...'
        write(LDT_logunit,*) LDT_rc%gridDesc(nest,7),LDT_rc%gridDesc(nest,4)
        call LDT_endrun
     endif
     if(LDT_rc%gridDesc(nest,8).lt.LDT_rc%gridDesc(nest,5)) then
        write(LDT_logunit,*) '[INFO] lon2 < lon1 ... ', &
                             LDT_rc%gridDesc(nest,8),LDT_rc%gridDesc(nest,5)
     endif
     
     ! Difference in number of longitudes (nx):
     LDT_rc%gridDesc(nest,2) = nint(diff_lon(LDT_rc%gridDesc(nest,8),LDT_rc%gridDesc(nest,5))&
                          / LDT_rc%gridDesc(nest,9)) + 1   

     ! Difference in number of latitudes (ny):
     LDT_rc%gridDesc(nest,3) = nint((LDT_rc%gridDesc(nest,7)-LDT_rc%gridDesc(nest,4))&
                          / LDT_rc%gridDesc(nest,10)) + 1    

     write(LDT_logunit,*)'-------------------- LDT/LIS Domain ----------------------'
     do k=1,13
        write(unit=LDT_logunit,fmt=*) '(',k,',',LDT_rc%gridDesc(nest,k),')'
     enddo
     
   ! local domain of each processor
     call map_set( PROJ_LATLON, LDT_rc%gridDesc(nest,4),LDT_rc%gridDesc(nest,5), &
                   0.0, LDT_rc%gridDesc(nest,9),LDT_rc%gridDesc(nest,10), 0.0,   &
                   LDT_rc%lnc(nest),LDT_rc%lnr(nest),LDT_domain(nest)%ldtproj )

   ! global domain
     call map_set( PROJ_LATLON,real(stlat),real(stlon), &
                   0.0, real(dx),real(dy), 0.0,   &
                   LDT_rc%gnc(nest),LDT_rc%gnr(nest),LDT_domain(nest)%ldtglbproj )
     
     write(LDT_logunit,*)'----------------------------------------------------------'

  enddo  ! end nest loop
  
  deallocate(run_dd)

end subroutine readinput_latlon
      
