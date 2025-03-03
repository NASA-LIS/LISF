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
! !ROUTINE: readinput_UTM
!  \label{readinput_UTM}
!  
! !REVISION HISTORY:
! 
!  28 Jun 2009: Sujay Kumar; Initial Specification
! !INTERFACE:
subroutine readinput_UTM(nest)

! !USES:
  use ESMF
  use LDT_coreMod,   only : LDT_rc, LDT_domain, LDT_config, &
       LDT_localPet, LDT_npes, LDT_ews_halo_ind, LDT_ewe_halo_ind, &
       LDT_nss_halo_ind, LDT_nse_halo_ind
  use LDT_domainMod, only : LDT_quilt_domain, LDT_setParameterDomainSpecs
  use LDT_logMod,    only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
  use map_utils

  integer, intent(in) :: nest

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LDT configuration 
!  file to determine the region of interest in the LDT simulation, in 
!  a lat/lon projection. The routine reads the extents of the 
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest. The domain 
!  specifications of the parameter maps are also read in by this routine. 
!  Based on the number of processors and the processor layout specified, 
!  this routine also performs the domain decomposition. 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[ldttask\_for\_point](\ref{LDT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
!EOP
  implicit none

  integer              :: i,k,n
  real, allocatable    :: run_dd(:,:)
  integer              :: nc, nr
  integer              :: rc
  real(ESMF_KIND_R8)   :: stnorthing, steasting, dx
! _____________________________________________________________________

  allocate(run_dd(LDT_rc%nnest,20))

  LDT_rc%lis_map_resfactor(nest) = 100000.   ! if in meters

  call ESMF_ConfigFindLabel(LDT_config,"Run domain UTM zone:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,1),rc=rc)
     call LDT_verify(rc,'Run domain UTM zone: not defined')
  enddo  

  call ESMF_ConfigFindLabel(LDT_config,"Run domain northing of SW corner:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,2),rc=rc)
     call LDT_verify(rc,'Run domain northing of SW corner: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain easting of SW corner:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,3),rc=rc)
     call LDT_verify(rc,'Run domain easting of SW corner: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain x-dimension size:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,4),rc=rc)
     call LDT_verify(rc,'Run domain x-dimension size: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"Run domain y-dimension size:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,5),rc=rc)
     call LDT_verify(rc,'Run domain y-dimension size: not defined')
  enddo  

! Resolution:  in meters?
  call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution:",rc=rc)
  do i=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(i,6),rc=rc)
     call LDT_verify(rc,'Run domain resolution: not defined')
  enddo  

  write(LDT_logunit,fmt=*) 'DOMAIN details:'

  if(LDT_rc%npesx*LDT_rc%npesy.ne.LDT_npes) then 
     write(unit=LDT_logunit,fmt=*) 'Layout does not match the number of processors...'
     write(unit=LDT_logunit,fmt=*) 'npex, npey, ',LDT_rc%npesx,'x',LDT_rc%npesy,'!=',LDT_npes
     write(unit=LDT_logunit,fmt=*) 'Stopping program ..'
     call LDT_endrun()
  endif

  do n=1,nest

     stnorthing = run_dd(n,2)
     steasting  = run_dd(n,3)

     dx = run_dd(n,6)
     nc = nint(run_dd(n,4))
     nr = nint(run_dd(n,5))

     call LDT_quilt_domain(n, nc,nr)
     
     LDT_rc%gridDesc(n,1) = 7

     LDT_rc%gridDesc(n,4) = stnorthing+ (LDT_nss_halo_ind(n,LDT_localPet+1)-1)*dx
     LDT_rc%gridDesc(n,5) = steasting + (LDT_ews_halo_ind(n,LDT_localPet+1)-1)*dx
     LDT_rc%gridDesc(n,7) = stnorthing+ (LDT_nse_halo_ind(n,LDT_localPet+1)-1)*dx
     LDT_rc%gridDesc(n,8) = steasting + (LDT_ewe_halo_ind(n,LDT_localPet+1)-1)*dx
     
     write(unit=LDT_logunit,fmt=*) 'local domain ',&
          LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,7),&
          LDT_rc%gridDesc(n,5),LDT_rc%gridDesc(n,8), 'd: ',n

     LDT_rc%gnc(n) = nint(run_dd(n,4))
     LDT_rc%gnr(n) = nint(run_dd(n,5))

     LDT_rc%gridDesc(n,9) = run_dd(n,6)  !resolution

     LDT_rc%gridDesc(n,10) = run_dd(n,1) !zone
     LDT_rc%gridDesc(n,6) = 128
     LDT_rc%gridDesc(n,11) = 64
     LDT_rc%gridDesc(n,20) = 64
     
     if(LDT_rc%gridDesc(n,7).lt.LDT_rc%gridDesc(n,4)) then
        write(unit=LDT_logunit,fmt=*) 'northing2 must be greater than northing1'
        write(LDT_logunit,*) LDT_rc%gridDesc(n,7),LDT_rc%gridDesc(n,4)
        write(unit=LDT_logunit,fmt=*) 'Stopping run...'
        call LDT_endrun
     endif
     if(LDT_rc%gridDesc(n,8).lt.LDT_rc%gridDesc(n,5)) then
        write(unit=LDT_logunit,fmt=*) 'easting2 must be greater than easting1'
        write(LDT_logunit,*) LDT_rc%gridDesc(n,8),LDT_rc%gridDesc(n,5)
        write(unit=LDT_logunit,fmt=*) 'Stopping run...'
        call LDT_endrun
     endif
     
     LDT_rc%gridDesc(n,2) = nint((LDT_rc%gridDesc(n,8)-LDT_rc%gridDesc(n,5))&
          /dx)+ 1
     LDT_rc%gridDesc(n,3) = nint((LDT_rc%gridDesc(n,7)-LDT_rc%gridDesc(n,4))&
          /dx) + 1  
     write(LDT_logunit,*)'-------------------- LDT/LDT Domain ----------------------'
     do k=1,13
        write(unit=LDT_logunit,fmt=*) '(',k,',',LDT_rc%gridDesc(n,k),')'
     enddo

     LDT_rc%pnc(n) = LDT_rc%gridDesc(n,32)
     LDT_rc%pnr(n) = LDT_rc%gridDesc(n,33)

     write(LDT_logunit,*)'----------------------------------------------------------'

  enddo
  deallocate(run_dd)

end subroutine readinput_UTM
      
