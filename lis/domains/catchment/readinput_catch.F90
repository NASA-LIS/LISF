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
! !ROUTINE: readinput_catch
!  \label{readinput_catch}
!  
! !REVISION HISTORY:
! 
!  14 Nov 2003: Sujay Kumar; Initial Code
!
! !INTERFACE:
subroutine readinput_catch
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config, LIS_localPet, LIS_masterproc, &
       LIS_npes,LIS_ews_halo_ind, LIS_nss_halo_ind, LIS_ewe_halo_ind, &
       LIS_nse_halo_ind
  use LIS_logMod, only :LIS_logunit, LIS_endrun
  use LIS_domainMod, only : LIS_quilt_domain, LIS_setParameterDomainSpecs
  use LIS_fileIOMod, only : LIS_readDomainConfigSpecs

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

  integer :: i
  integer :: k,n
  real, allocatable :: run_dd(:,:)
  real, allocatable :: param_dd(:,:)
  integer :: nc, nr
  integer :: rc
  real(ESMF_KIND_R8) :: stlat, stlon, dx, dy

  allocate(run_dd(LIS_rc%nnest,6))
  allocate(param_dd(LIS_rc%nnest,6))
  call ESMF_ConfigFindLabel(LIS_config,"run domain lower left lat:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,1),rc=rc)
  enddo  
  call ESMF_ConfigFindLabel(LIS_config,"run domain lower left lon:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,2),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"run domain upper right lat:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,3),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"run domain upper right lon:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,4),rc=rc)
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"run domain resolution (dx):",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,5),rc=rc)
  enddo  

  call ESMF_ConfigFindLabel(LIS_config,"run domain resolution (dy):",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,run_dd(i,6),rc=rc)
  enddo  

  call LIS_readDomainConfigSpecs("param domain", param_dd)

  write(LIS_logunit,fmt=*)'DOMAIN details:'

  if(LIS_rc%npesx*LIS_rc%npesy.ne.LIS_npes) then 
     write(unit=LIS_logunit,fmt=*) 'Layout does not match the number of processors...'
     write(unit=LIS_logunit,fmt=*) 'npex, npey, ',LIS_rc%npesx,'x',LIS_rc%npesy,'!=',LIS_npes
     write(unit=LIS_logunit,fmt=*) 'Stopping program ..'
     call LIS_endrun()
  endif

  do n=1,LIS_rc%nnest
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     dx = run_dd(n,5)
     dy = run_dd(n,6)

     nc = (nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))) + 1
     nr = (nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))) + 1

     call LIS_quilt_domain(n,nc,nr)
     
     LIS_rc%gridDesc(n,4) = stlat + (LIS_nss_halo_ind(n,LIS_localPet+1)-1)*dy
     LIS_rc%gridDesc(n,5) = stlon + (LIS_ews_halo_ind(n,LIS_localPet+1)-1)*dx
     LIS_rc%gridDesc(n,7) = stlat + (LIS_nse_halo_ind(n,LIS_localPet+1)-1)*dy
     LIS_rc%gridDesc(n,8) = stlon + (LIS_ewe_halo_ind(n,LIS_localPet+1)-1)*dx
     
     write(unit=LIS_logunit,fmt=*) 'local domain ',&
          LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,7),&
          LIS_rc%gridDesc(n,5),LIS_rc%gridDesc(n,8), 'd: ',n

     LIS_rc%gnc(n) = nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,5))+1
     LIS_rc%gnr(n) = nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,6))+1
     LIS_rc%gridDesc(n,1) = 0
     LIS_rc%gridDesc(n,9) = run_dd(n,5)
     if(LIS_rc%gridDesc(n,1).eq.0) then 
        LIS_rc%gridDesc(n,10) = run_dd(n,6)
     endif
  
     call LIS_setParameterDomainSpecs(param_dd(n,:),LIS_rc%gridDesc(n,:))

     if(LIS_rc%gridDesc(n,1).eq.0) then 
        LIS_rc%gridDesc(n,6) = 128
        LIS_rc%gridDesc(n,11) = 64
        LIS_rc%gridDesc(n,20) = 64
     endif
     if(LIS_rc%gridDesc(n,7).le.LIS_rc%gridDesc(n,4)) then
        write(unit=LIS_logunit,fmt=*) 'lat2 must be greater than lat1'
        write(unit=LIS_logunit,fmt=*) 'Stopping run...'
        call LIS_endrun
     endif
     if(LIS_rc%gridDesc(n,8).le.LIS_rc%gridDesc(n,5)) then
        write(unit=LIS_logunit,fmt=*) 'lon2 must be greater than lon1'
        write(unit=LIS_logunit,fmt=*) 'Stopping run...'
        call LIS_endrun
     endif
     
     LIS_rc%gridDesc(n,2) = nint((LIS_rc%gridDesc(n,8)-LIS_rc%gridDesc(n,5))&
          /LIS_rc%gridDesc(n,9))+ 1
     LIS_rc%gridDesc(n,3) = nint((LIS_rc%gridDesc(n,7)-LIS_rc%gridDesc(n,4))&
          /LIS_rc%gridDesc(n,10)) + 1

     write(LIS_logunit,*)'--------------------Domain ',n,'----------------------'
     do k=1,13
        write(unit=LIS_logunit,fmt=*) '(',k,',',LIS_rc%gridDesc(n,k),')'
     enddo
     do k=30,40
        write(unit=LIS_logunit,fmt=*) '(',k,',',LIS_rc%gridDesc(n,k),')'
     enddo
     write(LIS_logunit,*)'--------------------------------------------------------------'

     LIS_rc%pnc(n) = LIS_rc%gridDesc(n,32)
     LIS_rc%pnr(n) = LIS_rc%gridDesc(n,33)

  enddo

  deallocate(run_dd)
  deallocate(param_dd)

  call ESMF_ConfigFindLabel(LIS_config,"bottom temperature map:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%tbotfile(i),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"greenness maximum map:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%shdmaxfile(i),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"greenness minimum map:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%shdminfile(i),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"slope type map:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%slopetypefile(i),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"tile coord file:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%tile_coord_file(i),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"tile veg file:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%tile_veg_file(i),rc=rc)
  enddo

  return

end subroutine readinput_catch
      










