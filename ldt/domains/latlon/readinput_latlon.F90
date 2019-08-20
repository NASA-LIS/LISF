!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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
subroutine readinput_latlon

! !USES:
  use ESMF 
  use LDT_coreMod
  use LDT_domainMod
  use LDT_logMod
  use LDT_fileIOMod
  use map_utils

  implicit none

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
!  \item[listask\_for\_point](\ref{LDT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
!EOP

  integer              :: i, j, k, n, count
  real, allocatable    :: run_dd(:,:)  ! LIS Target grid/domain
  integer              :: lnc,lnr
  integer              :: nc, nr
  integer              :: ierr, rc
  integer              :: ips, ipe, jps, jpe
  integer              :: Px, Py, P
  integer              :: mytask_x, mytask_y
  real(ESMF_KIND_R8)   :: stlat, stlon, dx, dy
  integer              :: status
! __________________________________________________________________________

  allocate(run_dd(LDT_rc%nnest,20))

  LDT_rc%lis_map_resfactor = 1.

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lat:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,1),rc=rc)
     call LDT_verify(rc, 'Run domain lower left lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lon:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,2),rc=rc)
     call LDT_verify(rc, 'Run domain lower left lon: not defined')
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lat:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,3),rc=rc)
     call LDT_verify(rc, 'Run domain upper right lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain upper right lon:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,4),rc=rc)
     call LDT_verify(rc, 'Run domain upper right lon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution (dx):",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,5),rc=rc)
     call LDT_verify(rc, 'Run domain resolution (dx): not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution (dy):",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,6),rc=rc)
     call LDT_verify(rc, 'Run domain resolution (dy): not defined')
  enddo

  do n=1,LDT_rc%nnest
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     dx    = run_dd(n,5)
     dy    = run_dd(n,6)
     nc    = (nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))) + 1
     nr    = (nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))) + 1
     
     call LDT_quilt_domain(n,nc,nr)

     LDT_rc%gridDesc(n,4) = stlat + (LDT_nss_halo_ind(n,LDT_localPet+1)-1)*dy
     LDT_rc%gridDesc(n,5) = stlon + (LDT_ews_halo_ind(n,LDT_localPet+1)-1)*dx
     LDT_rc%gridDesc(n,7) = stlat + (LDT_nse_halo_ind(n,LDT_localPet+1)-1)*dy
     LDT_rc%gridDesc(n,8) = stlon + (LDT_ewe_halo_ind(n,LDT_localPet+1)-1)*dx

     write(unit=LDT_logunit,fmt=*) 'local domain ',&
          LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,7),&
          LDT_rc%gridDesc(n,5),LDT_rc%gridDesc(n,8)
     
     LDT_rc%gnc(n) = nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))+1
     LDT_rc%gnr(n) = nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))+1
     LDT_rc%gridDesc(n,1) = 0
     LDT_rc%gridDesc(n,9) = run_dd(n,5)
     
     if(LDT_rc%gridDesc(n,1).eq.0) then 
        LDT_rc%gridDesc(n,10) = run_dd(n,6)
        LDT_rc%gridDesc(n,6)  = 128
        LDT_rc%gridDesc(n,11) = 64
        LDT_rc%gridDesc(n,20) = 64
     endif
     if(LDT_rc%gridDesc(n,7).lt.LDT_rc%gridDesc(n,4)) then
        write(unit=LDT_logunit,fmt=*) 'lat2 must be greater than lat1'
        write(LDT_logunit,*) LDT_rc%gridDesc(n,7),LDT_rc%gridDesc(n,4)
        write(unit=LDT_logunit,fmt=*) 'Stopping run...'
        call LDT_endrun
     endif
     if(LDT_rc%gridDesc(n,8).lt.LDT_rc%gridDesc(n,5)) then
        write(unit=LDT_logunit,fmt=*) 'lon2 must be greater than lon1'
        write(LDT_logunit,*) LDT_rc%gridDesc(n,8),LDT_rc%gridDesc(n,5)
        write(unit=LDT_logunit,fmt=*) 'Stopping run...'
        call LDT_endrun
     endif
     
     LDT_rc%gridDesc(n,2) = nint((LDT_rc%gridDesc(n,8)-LDT_rc%gridDesc(n,5))&
          /LDT_rc%gridDesc(n,9))+ 1      ! dx
     LDT_rc%gridDesc(n,3) = nint((LDT_rc%gridDesc(n,7)-LDT_rc%gridDesc(n,4))&
          /LDT_rc%gridDesc(n,10)) + 1    ! dy

     write(LDT_logunit,*)'-------------------- LDT/LIS Domain ----------------------'
     do k=1,13
        write(unit=LDT_logunit,fmt=*) '(',k,',',LDT_rc%gridDesc(n,k),')'
     enddo
     
   ! local domain of each processor
     call map_set( PROJ_LATLON, LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,5), &
                   0.0, LDT_rc%gridDesc(n,9),LDT_rc%gridDesc(n,10), 0.0,   &
                   LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_domain(n)%ldtproj )

   ! global domain
     call map_set( PROJ_LATLON,real(stlat),real(stlon), &
                   0.0, real(dx),real(dy), 0.0,   &
                   LDT_rc%gnc(n),LDT_rc%gnr(n),LDT_domain(n)%ldtglbproj )
     
     write(LDT_logunit,*)'----------------------------------------------------------'

  enddo  ! end nest loop
  
  deallocate(run_dd)

end subroutine readinput_latlon
      
