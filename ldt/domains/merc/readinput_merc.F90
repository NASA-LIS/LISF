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
! !ROUTINE: readinput_merc
! \label{redinput_merc}
! 
! !REVISION HISTORY:
!  15 Jul 2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine readinput_merc(nest)
! !USES:
  use ESMF 
  use map_utils
  use LDT_coreMod, only :LDT_rc, LDT_domain, LDT_config, LDT_localPet, &
          LDT_masterproc, LDT_npes, &
          LDT_ews_halo_ind, LDT_nss_halo_ind, LDT_ewe_halo_ind, LDT_nse_halo_ind
  use LDT_domainMod, only : LDT_quilt_domain
  use LDT_fileIOMod, only : LDT_readDomainConfigSpecs       
  use LDT_logMod,    only : LDT_logunit, LDT_endrun, LDT_verify


  implicit none 

  integer, intent(in) :: nest

! !DESCRIPTION: 
!
!  This routine reads the options specified in the LDT configuration 
!  file to determine the region of interest in the LDT simulation, in 
!  a mercator projection. The routine reads the extents of the 
!  running domain. The land surface parameters used for the simulation
!  are read from a file that spans the area of interest. The domain 
!  specifications of the parameter maps are also read in by this routine. 
!  Based on the number of processors and the processor layout specified, 
!  this routine also performs the domain decomposition. 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[listask\_for\_point](\ref{LDT_mpDecomp}) \newline
!   routine to perform domain decomposition
!  \end{description}
!EOP

  integer           :: i
  integer           :: k
  real, allocatable :: run_dd(:,:)
  real              :: lat_str,lon_str
  real              :: stlat,stlon,truelat1,stdlon,xmesh
  integer           :: nc,nr
  type(proj_info)   :: proj_temp
  integer           :: n, rc
  integer           :: ios
  integer           :: ftn
  logical           :: file_exists
  integer           :: latId, lonId
  integer           :: ncId, nrId
  integer, allocatable  :: gnc(:), gnr(:)
  real, allocatable     :: lat(:,:)
  real, allocatable     :: lon(:,:)
! _________________________________________________

  allocate(run_dd(LDT_rc%nnest,7))
  allocate(gnc(LDT_rc%nnest))
  allocate(gnr(LDT_rc%nnest))

  LDT_rc%lis_map_resfactor(nest) = 100.

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lat:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,1),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain lower left lon:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,2),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain true lat1:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,3),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain standard lon:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,4),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain resolution:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,5),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain x-dimension size:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,6),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"Run domain y-dimension size:",rc=rc)
  do n=1,nest
     call ESMF_ConfigGetAttribute(LDT_config,run_dd(n,7),rc=rc)
  enddo

  do n=1,nest

     write(LDT_logunit,*) "DOMAIN details: Running in Mercator projection."
     if(LDT_npes.eq.1) then 
        LDT_rc%npesx = 1
        LDT_rc%npesy = 1
     endif
     
     if(LDT_rc%npesx*LDT_rc%npesy.ne.LDT_npes) then 
        write(LDT_logunit,*) "Layout does not match the number of processors..."
        write(LDT_logunit,*) "npex, npey, ",LDT_rc%npesx,"x",LDT_rc%npesy,"!=",LDT_npes
        write(LDT_logunit,*) "Stopping program ..."
        call LDT_endrun()
     endif
     
     stlat    = run_dd(n,1)
     stlon    = run_dd(n,2)
     truelat1 = run_dd(n,3)
     stdlon   = run_dd(n,4)
     xmesh    = run_dd(n,5)
     nc       = nint(run_dd(n,6))
     nr       = nint(run_dd(n,7))

     call LDT_quilt_domain(n,nc,nr)

   ! Call map_set on the overall running domain:
     call map_set(PROJ_MERC,stlat,stlon,&
          xmesh*1000.0,stdlon,truelat1,&
          0.0,nc,nr,proj_temp)

     call ij_to_latlon(proj_temp,float(LDT_ews_halo_ind(n,LDT_localPet+1)),&
          float(LDT_nss_halo_ind(n,LDT_localPet+1)),lat_str,lon_str)

     LDT_rc%gridDesc(n,1) = 1
     LDT_rc%gridDesc(n,2) = LDT_rc%lnc(n)
     LDT_rc%gridDesc(n,3) = LDT_rc%lnr(n)
     LDT_rc%gridDesc(n,4) = lat_str
     LDT_rc%gridDesc(n,5) = lon_str
     LDT_rc%gridDesc(n,6) = 8
     LDT_rc%gridDesc(n,7) = 0.0
     LDT_rc%gridDesc(n,8) = xmesh
     LDT_rc%gridDesc(n,9) = xmesh
     LDT_rc%gridDesc(n,10) = truelat1
     LDT_rc%gridDesc(n,11) = stdlon
     
!------------------------------------------------------------------------
! Read namelist of parameters depending on the domain
!------------------------------------------------------------------------
! for now set the global domain the same as local domain
! need to be modified to work on subsetted domains

     write(LDT_logunit,*)'-------------------- LDT/LIS Domain ----------------------'
     do k=1,13
        write(LDT_logunit,*)'(',k,',',LDT_rc%gridDesc(n,k),')'
     enddo
     
     LDT_rc%gnc(n) = nc
     LDT_rc%gnr(n) = nr
     
     call map_set(PROJ_MERC,LDT_rc%gridDesc(n,4),LDT_rc%gridDesc(n,5),&
          LDT_rc%gridDesc(n,8)*1000.0,LDT_rc%gridDesc(n,11),&
          LDT_rc%gridDesc(n,10),&
          0.0,LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_domain(n)%ldtproj)
     write(LDT_logunit,*)'----------------------------------------------------------'

     ! EMK...Copy to global structure
     LDT_domain(n)%ldtglbproj = LDT_domain(n)%ldtproj

  enddo
  deallocate(run_dd)

end subroutine readinput_merc
      
