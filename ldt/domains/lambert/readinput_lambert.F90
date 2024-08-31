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
! !ROUTINE: readinput_lambert
!  \label{readinput_lambert}
! 
! !REVISION HISTORY:
!  15Jul2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine readinput_lambert(nest)
! !USES:
  use ESMF 
  use LDT_logMod, only : LDT_logunit, LDT_endrun, LDT_verify
  use map_utils
  use LDT_coreMod, only : LDT_rc, LDT_domain, LDT_config, &
       LDT_localPet, LDT_masterproc, LDT_npes, &
       LDT_ews_halo_ind, LDT_nss_halo_ind, LDT_ewe_halo_ind, &
       LDT_nse_halo_ind
  use LDT_domainMod, only : LDT_quilt_domain, LDT_setParameterDomainSpecs
  use LDT_fileIOMod, only : LDT_readDomainConfigSpecs

  implicit none

  integer, intent(in) :: nest

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LDT configuration 
!  file to determine the region of interest in the LDT simulation, in 
!  a lambert conformal projection. The routine reads the extents of the 
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

  integer           :: i, k
  real, allocatable :: run_dd(:,:)
  real              :: lat_str,lon_str
  real              :: stlat,stlon,truelat1,truelat2,stdlon,xmesh
  integer           :: nc,nr
  type(proj_info)   :: proj_temp
  integer           :: n, rc
  integer           :: nDupl
  logical           :: checkDuplicate
  logical           :: newFormatCheck 
! __________________________________________

  allocate(run_dd(LDT_rc%nnest,20))

  newFormatCheck = .false. 
  checkDuplicate = .false. 
  nDupl = 1
  do n=2,LDT_rc%nnest
     if(LDT_rc%lis_map_proj(n).eq.LDT_rc%lis_map_proj(1)) then 
        checkDuplicate = .true. 
        nDupl = nDupl + 1
     endif
  enddo

  LDT_rc%lis_map_resfactor(nest) = 100.

  call ESMF_ConfigFindLabel(LDT_config,"LM run domain lower left lat:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,1),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"LM run domain lower left lon:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,2),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"LM run domain true lat1:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,3),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"LM run domain true lat2:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,4),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"LM run domain standard lon:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,5),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"LM run domain resolution:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,6),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"LM run domain x-dimension size:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,7),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  call ESMF_ConfigFindLabel(LDT_config,"LM run domain y-dimension size:",rc=rc)
  do n=1,nDupl
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,8),rc=rc)
     if(rc.eq.0) newFormatCheck = .true. 
  enddo
  
  if(.not.newFormatCheck) then 
     call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lat:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,1),rc=rc)
        call LDT_verify(rc,'Run domain lower left lat: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lon:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,2),rc=rc)
        call LDT_verify(rc,'Run domain lower left lon: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain true lat1:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,3),rc=rc)
        call LDT_verify(rc,'Run domain true lat1: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain true lat2:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,4),rc=rc)
        call LDT_verify(rc,'Run domain true lat2: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain standard lon:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,5),rc=rc)
        call LDT_verify(rc,'Run domain standard lon: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,6),rc=rc)
        call LDT_verify(rc,'Run domain resolution: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain x-dimension size:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,7),rc=rc)
        call LDT_verify(rc,'Run domain x-dimension size: not defined')
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Run domain y-dimension size:",rc=rc)
     do n=1,nDupl
        call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,8),rc=rc)
        call LDT_verify(rc,'Run domain y-dimension size: not defined')
     enddo
  endif

  write(LDT_logunit,*)'DOMAIN details: Running in lambert projection'
  if(LDT_npes.eq.1) then 
     LDT_rc%npesx = 1
     LDT_rc%npesy = 1
  endif

  if(LDT_rc%npesx*LDT_rc%npesy.ne.LDT_npes) then 
     write(LDT_logunit,*) 'Layout does not match the number of processors...'
     write(LDT_logunit,*) 'npex, npey, ',LDT_rc%npesx,'x',LDT_rc%npesy,'!=',LDT_npes
     write(LDT_logunit,*) 'Stopping program ..'
     call LDT_endrun()
  endif

  do n = 1,nDupl
     stlat    = run_dd(n,1)
     stlon    = run_dd(n,2)
     truelat1 = run_dd(n,3)
     truelat2 = run_dd(n,4)
     stdlon   = run_dd(n,5)
     xmesh    = run_dd(n,6)
     nc       = nint(run_dd(n,7))
     nr       = nint(run_dd(n,8))
  
     call LDT_quilt_domain(nest,nc,nr)
  
!call map_set on the overall running domain

     call map_set( PROJ_LC, stlat, stlon,&
              xmesh*1000.0, stdlon, truelat1,&
              truelat2, nc, nr, proj_temp)

     call ij_to_latlon(proj_temp,float(LDT_ews_halo_ind(nest,LDT_localPet+1)),&
          float(LDT_nss_halo_ind(nest,LDT_localPet+1)),lat_str,lon_str)     
     
     LDT_rc%gridDesc(nest,1) = 3
     LDT_rc%gridDesc(nest,2) = LDT_rc%lnc(nest)
     LDT_rc%gridDesc(nest,3) = LDT_rc%lnr(nest)
     LDT_rc%gridDesc(nest,4) = lat_str
     LDT_rc%gridDesc(nest,5) = lon_str
     LDT_rc%gridDesc(nest,6) = 8
     LDT_rc%gridDesc(nest,7) = truelat2
     LDT_rc%gridDesc(nest,8) = xmesh
     LDT_rc%gridDesc(nest,9) = xmesh
     LDT_rc%gridDesc(nest,10) = truelat1
     LDT_rc%gridDesc(nest,11) = stdlon
     
!------------------------------------------------------------------------
! Read namelist of parameters depending on the domain
!------------------------------------------------------------------------
! for now set the global domain the same as local domain
! need to be modified to work on subsetted domains

     write(LDT_logunit,*)'-------------------- LDT/LIS Domain ----------------------'
     do k=1,13
        write(LDT_logunit,*)'(',k,',',LDT_rc%gridDesc(nest,k),')'
     enddo
     
     LDT_rc%gnc(nest) = nc
     LDT_rc%gnr(nest) = nr
     
     call map_set(PROJ_LC, LDT_rc%gridDesc(nest,4), LDT_rc%gridDesc(nest,5),&
              LDT_rc%gridDesc(nest,8)*1000.0, LDT_rc%gridDesc(nest,11),&
              LDT_rc%gridDesc(nest,10), truelat2, &
              LDT_rc%lnc(nest), LDT_rc%lnr(nest), LDT_domain(nest)%ldtproj)
     write(LDT_logunit,*)'----------------------------------------------------------'

     ! EMK: Copy to global structure.
     LDT_domain(nest)%ldtglbproj = LDT_domain(nest)%ldtproj

  enddo
  deallocate(run_dd)

 end subroutine readinput_lambert

