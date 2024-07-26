!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readinput_lambert
!  \label{readinput_lambert}
! 
! !REVISION HISTORY:
!  15Jul2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine readinput_lambert
! !USES:
  use ESMF 
  use map_utils
  use preprocMod

  implicit none

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
!EOP

  integer           :: i
  integer           :: k
  real, allocatable :: run_dd(:,:)
!  real, allocatable :: param_dd(:,:)
  real              :: lat_str,lon_str
  real              :: stlat,stlon,truelat1,truelat2,stdlon,xmesh
  integer           :: nc,nr
  type(proj_info)   :: proj_temp
  integer           :: n, rc


  allocate(run_dd(LVT_rc%nnest,20))
!  allocate(param_dd(LVT_rc%nnest,6))

  call ESMF_ConfigFindLabel(LVT_config,"run domain lower left lat:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,1),rc=rc)
     call LVT_verify(rc,'run domain lower left lat: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain lower left lon:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,2),rc=rc)
     call LVT_verify(rc,'run domain lower left lon: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LVT_config,"run domain true lat1:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,3),rc=rc)
     call LVT_verify(rc,'run domain true lat1: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain true lat2:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,4),rc=rc)
     call LVT_verify(rc,'run domain true lat2: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain standard lon:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,5),rc=rc)
     call LVT_verify(rc,'run domain standard lon: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain resolution:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,6),rc=rc)
     call LVT_verify(rc,'run domain resolution: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain x-dimension size:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,7),rc=rc)
     call LVT_verify(rc,'run domain x-dimension size: not defined')
  enddo

  call ESMF_ConfigFindLabel(LVT_config,"run domain y-dimension size:",rc=rc)
  do n=1,LVT_rc%nnest
     call ESMF_ConfigGetAttribute(LVT_config,run_dd(n,8),rc=rc)
     call LVT_verify(rc,'run domain y-dimension size: not defined')
  enddo

!  call LVT_readDomainConfigSpecs("param domain",param_dd)

  write(*,*)'DOMAIN details: Running in lambert projection'

  do n=1,LVT_rc%nnest
     stlat = run_dd(n,1)
     stlon = run_dd(n,2)
     truelat1 = run_dd(n,3)
     truelat2 = run_dd(n,4)
     stdlon = run_dd(n,5)
     xmesh = run_dd(n,6)
     nc = nint(run_dd(n,7))
     nr = nint(run_dd(n,8))

!call map_set on the overall running domain

     call map_set(PROJ_LC,stlat,stlon,&
          xmesh*1000.0,stdlon,truelat1,&
          truelat2,nc,nr,proj_temp)
     call ij_to_latlon(proj_temp,1.0,1.0,&
          lat_str,lon_str)     

     LVT_rc%lnc(n) = nc
     LVT_rc%lnr(n) = nr
     
     LVT_rc%gridDesc(n,1) = 3
     LVT_rc%gridDesc(n,2) = LVT_rc%lnc(n)
     LVT_rc%gridDesc(n,3) = LVT_rc%lnr(n)
     LVT_rc%gridDesc(n,4) = lat_str
     LVT_rc%gridDesc(n,5) = lon_str
     LVT_rc%gridDesc(n,6) = 8
     LVT_rc%gridDesc(n,7) = truelat2
     LVT_rc%gridDesc(n,8) = xmesh
     LVT_rc%gridDesc(n,9) = xmesh
     LVT_rc%gridDesc(n,10) = truelat1
     LVT_rc%gridDesc(n,11) = stdlon
     
!------------------------------------------------------------------------
! Read namelvtt of parameters depending on the domain
!------------------------------------------------------------------------
! for now set the global domain the same as local domain
! need to be modified to work on subsetted domains

!     call LVT_setParameterDomainSpecs(param_dd(n,:),LVT_rc%gridDesc(n,:))
     
     do k=1,13
        write(*,*)'(',k,',',LVT_rc%gridDesc(n,k),')'
     enddo
     
!     do k=30,40
!        write(*,*)'(',k,',',LVT_rc%gridDesc(n,k),')'
!     enddo
     
     LVT_rc%gnc(n) = nc
     LVT_rc%gnr(n) = nr
!     LVT_rc%pnc(n) = LVT_rc%gridDesc(n,32)
!     LVT_rc%pnr(n) = LVT_rc%gridDesc(n,33)
     
     call map_set(PROJ_LC,LVT_rc%gridDesc(n,4),LVT_rc%gridDesc(n,5),&
          LVT_rc%gridDesc(n,8)*1000.0,LVT_rc%gridDesc(n,11),&
          LVT_rc%gridDesc(n,10),&
          truelat2,LVT_rc%lnc(n),LVT_rc%lnr(n),LVT_domain(n)%lvtproj)
  enddo
  deallocate(run_dd)
!  deallocate(param_dd)
     
 end subroutine readinput_lambert

