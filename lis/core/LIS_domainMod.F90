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
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_domainMod
!BOP
!
! !MODULE: LIS_domainMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage the creation of 
!   the domain (grid, tile, patch and ensemble spaces) used in the LIS run
!
!   \subsubsection{Overview}
!   The `domain' is defined in LIS as a cartesian grid in two dimensions
!   with an ordering from the lower left corner to the upper right corner. 
!   This rule is followed for the running domain of LIS that may use different
!   map projections. 
!  
!   When the running domain is created, LIS generates a 1-d vector 
!   representation of the 2-d input grid based on the input landmask (where
!   open water points are typically excluded). The basic computational unit 
!   in LIS is called a ``tile''. If options for subgrid variability are
!   employed, LIS will create multiple number of tiles per each grid cell. 
!   The subgrid tiles can be defined based on the distribution of 
!   vegetation types, soil characteristics or topography. Each tile can 
!   further define a specified number of ensembles.
!
!   LIS also supports the use of multiple surface types. For example, 
!   the domain could consist of `land' surfaces and `lake' surfaces and
!   LSMs will be executed over land points and Lake models will be executed
!   over lake points. Within LIS, these surface models are defined to 
!   run over ``patches''. The patches across different surface models 
!   sum to the tile dimension. 
!    
! 
! !REVISION HISTORY: 
!  17 Feb 2004    Sujay Kumar  Initial Specification
!  24 Aug 2008    Sujay Kumar  Implemented halo support 
!   3 Aug 2012    Sujay Kumar  Added support for flexible tiling
!   3 Mar 2022    Kristi Arsenault  Added support for curvature tiles
! 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_LMLCMod
  use LIS_topoMod
  use LIS_soilsMod
  use LIS_histDataMod
  use LIS_mpiMod
  use map_utils
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_domain_init          !initialize specified domains
  public :: LIS_domain_setup        !setup domain related structures
  public :: LIS_quilt_domain         !generate quilted domains
  public :: LIS_domain_finalize      !cleanup allocated structures
  public :: decompose_nx_ny          !decomposes domain based on proc layout
  public :: decompose_npes           !decomposes domain based on proc elements
!EOP

contains
!BOP
! !ROUTINE: LIS_domain_init
! \label{LIS_domain_init}
!
! !INTERFACE: 
  subroutine LIS_domain_init
! !USES:
    !NONE
!
! !DESCRIPTION: 
!  This routine invokes the registry that defines the domain implementations, 
!  first. This is followed by calling the routines to read the runtime specific 
!  parameters for each domain instance. A call to create the 
!  domain is issued next, which is expected to generate the grid, tile, patch
!  and ensemble spaces for each domain. This routine is also expected to perform
!  any domain decomposition needed for parallel processing. Finally, the domain
!  decomposition information is disseminated to individual processors. 
!
! The calling sequence is:
! \begin{description}
!  \item[LIS\_domain\_plugin] (\ref{LIS_domain_plugin}) \newline
!    sets up function table registries for implemented domains
!  \item[readDomainInput] (\ref{readDomainInput}) \newline
!    invokes the generic method in the registry to read domain specific
!    input options from the configuration file
!  \item[LIS\_LMLC\_init] (\ref{LIS_LMLC_init}) \newline
!    initialize landmask and landcover data structures
!  \item[LIS\_topo\_init] (\ref{LIS_topo_init}) \newline
!    initialize topography data structures
!   \item[LIS\_soils\_init](\ref{LIS_soils_init}) \newline
!    initialize soils data structures
!  \item[make\_domain] (\ref{make_domain}) \newline
!    issues the call to generate the tile, patch, and grid dimensions, 
!    and associated datastructures. 
!  \end{description}
!EOP
    integer :: ierr
    integer :: i
    integer :: n
    integer :: k
    integer :: m
    integer, allocatable :: deblklist(:,:,:)
    integer :: stid, enid
    integer :: status
    integer :: gindex, ntiles
    
    type(ESMF_DistGrid) :: tileDG
    type(ESMF_DistGrid) :: gridDG
    type(ESMF_DistGrid) :: patchDG(LIS_rc%max_model_types)

    integer             :: ntiless
    integer             :: ngrids
    integer             :: tdeltas
    integer             :: gdeltas
    integer             :: npatches(LIS_rc%max_model_types)
    integer             :: patch_deltas(LIS_rc%max_model_types)

    TRACE_ENTER("dom_init")
!    call LIS_domain_plugin
!    call readinput(trim(LIS_rc%lis_map_proj)//char(0))
    call readDomainInput()
    call LIS_LMLC_init()
    call LIS_topo_init()
    call LIS_soils_init()
    call make_domain()

    do n=1, LIS_rc%nnest
!       deallocate(LIS_LMLC(n)%landmask)
!       deallocate(LIS_LMLC(n)%landcover) ! needed by VIC.4.1.1
       deallocate(LIS_LMLC(n)%surfacetype)
       if(LIS_rc%usetexturemap(n).ne."none") then 
          deallocate(LIS_soils(n)%texture)
       endif
       if(LIS_rc%usesoilfractionmap(n).ne."none") then 
          deallocate(LIS_soils(n)%sand)
          deallocate(LIS_soils(n)%clay)
          deallocate(LIS_soils(n)%silt)
          deallocate(LIS_soils(n)%soilffgrd)
       endif
       if(LIS_rc%useelevationmap(n).ne."none") then 
!          deallocate(LIS_topo(n)%elevation)
          deallocate(LIS_topo(n)%elevfgrd)
       endif
       if(LIS_rc%useslopemap(n).ne."none") then 
          deallocate(LIS_topo(n)%slope)
          deallocate(LIS_topo(n)%slopefgrd)
       endif
       if(LIS_rc%useaspectmap(n).ne."none") then 
          deallocate(LIS_topo(n)%aspect)
          deallocate(LIS_topo(n)%aspectfgrd)
       endif          
       if(LIS_rc%usecurvaturemap(n).ne."none") then
          deallocate(LIS_topo(n)%curvature)
          deallocate(LIS_topo(n)%curvfgrd)
       endif
    enddo

    do n=1,LIS_rc%nnest
       do k = 1, LIS_rc%ngrid(n)
          allocate(LIS_domain(n)%grid(k)%subgrid_tiles(&
               LIS_rc%surface_maxt*&
               LIS_rc%soilt_maxt*&
               LIS_rc%elev_maxt*&
               LIS_rc%slope_maxt*&
               LIS_rc%aspect_maxt*&
               LIS_rc%nensem(n)))
          LIS_domain(n)%grid(k)%ntiles = 0
       enddo
       do k = 1, LIS_rc%ntiles(n)
          gindex = LIS_domain(n)%tile(k)%index
          LIS_domain(n)%grid(gindex)%ntiles = &
             LIS_domain(n)%grid(gindex)%ntiles + 1
          ntiles = LIS_domain(n)%grid(gindex)%ntiles
          LIS_domain(n)%grid(gindex)%subgrid_tiles(ntiles) = k
       enddo
    enddo

    do n=1,LIS_rc%nnest
       ntiless = LIS_rc%ntiles(n)
       tdeltas = LIS_rc%ntiles(n)

       do m=1,LIS_rc%max_model_types
          npatches(m)     = LIS_rc%npatch(n,m)
          patch_deltas(m) = LIS_rc%npatch(n,m)
       enddo

       ngrids = LIS_rc%ngrid(n)
       gdeltas = LIS_rc%ngrid(n)

!-----------------------------------------------------------------------------
!  The grid, tile space sizes of the decomposed domains are gathered 
!  from each processor to compute the total sizes of the entire domain. 
!-----------------------------------------------------------------------------
#if (defined SPMD)
       call MPI_ALLGATHER(ntiless,1,MPI_INTEGER,&
            LIS_ntiless(n,:),1,MPI_INTEGER,&
            LIS_mpi_comm,ierr)
       call MPI_ALLGATHER(ngrids,1,MPI_INTEGER,&
            LIS_ngrids(n,:),1,MPI_INTEGER,&
            LIS_mpi_comm,ierr)
       call MPI_ALLGATHER(tdeltas,1,MPI_INTEGER,&
            LIS_tdeltas(n,:),1,MPI_INTEGER,&
            LIS_mpi_comm,ierr)
       call MPI_ALLGATHER(gdeltas,1,MPI_INTEGER,&
            LIS_gdeltas(n,:),1,MPI_INTEGER,&
            LIS_mpi_comm,ierr)

       do m=1,LIS_rc%max_model_types
          call MPI_ALLGATHER(npatches(m),1,MPI_INTEGER,&
               LIS_npatches(n,m,:),1,MPI_INTEGER, &
               LIS_mpi_comm,ierr)
       enddo
       do m=1,LIS_rc%max_model_types
          call MPI_ALLGATHER(patch_deltas(m),1,MPI_INTEGER,&
               LIS_patch_deltas(n,m,:),1,MPI_INTEGER, &
               LIS_mpi_comm,ierr)
       enddo
#else
       LIS_ntiless(n,:) = ntiless
       LIS_ngrids(n,:) = ngrids
       LIS_tdeltas(n,:) = tdeltas
       LIS_gdeltas(n,:) = gdeltas
       
       do m=1,LIS_rc%max_model_types
          LIS_npatches(n,m,:) = npatches(m)
       enddo
       do m=1,LIS_rc%max_model_types
          LIS_patch_deltas(n,m,:) = patch_deltas(m)
       enddo
#endif
       LIS_rc%glbntiles(n) = 0 
       do i=0,LIS_npes-1
          LIS_rc%glbntiles(n) = LIS_rc%glbntiles(n) + LIS_ntiless(n,i)
       enddo
       
       LIS_rc%glbngrid(n) = 0 
       do i=0,LIS_npes-1
          LIS_rc%glbngrid(n) = LIS_rc%glbngrid(n) + LIS_ngrids(n,i)
       enddo

       do m=1,LIS_rc%max_model_types
          LIS_rc%glbnpatch(n,m) = 0 
          do i=0,LIS_npes-1
             LIS_rc%glbnpatch(n,m) = LIS_rc%glbnpatch(n,m)+ LIS_npatches(n,m,i)
          enddo
       enddo

       if(LIS_masterproc) then 
          LIS_toffsets(n,0) = 0 
          do i=1,LIS_npes-1
             LIS_toffsets(n,i) = LIS_toffsets(n,i-1)+LIS_tdeltas(n,i-1)
          enddo
          LIS_goffsets(n,0) = 0 
          do i=1,LIS_npes-1
             LIS_goffsets(n,i) = LIS_goffsets(n,i-1)+LIS_gdeltas(n,i-1)
          enddo
          LIS_patch_offsets(n,:,0) =0
          do i=1,LIS_npes-1
             do m=1,LIS_rc%max_model_types
                LIS_patch_offsets(n,m,i) = LIS_patch_offsets(n,m,i-1)+&
                     LIS_patch_deltas(n,m,i-1)
             enddo
          enddo

       end if

#if (defined SPMD)
       call MPI_BCAST(LIS_toffsets(n,:), LIS_npes, MPI_INTEGER,0, &
            LIS_mpi_comm, ierr)
       call MPI_BCAST(LIS_goffsets(n,:), LIS_npes, MPI_INTEGER,0, &
            LIS_mpi_comm, ierr)
       do m=1,LIS_rc%max_model_types
          call MPI_BCAST(LIS_patch_offsets(n,m,:), LIS_npes, MPI_INTEGER,0, &
               LIS_mpi_comm, ierr)
       enddo
#endif

       allocate(deblklist(1,2,LIS_npes))      
       do i=0,LIS_npes-1
          stid = LIS_toffsets(n,i)+1
          enid = stid + LIS_ntiless(n,i)-1

          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)          
       enddo

       tileDG = ESMF_DistGridCreate(minIndex=(/1/), &
            maxIndex=(/LIS_rc%glbntiles(n)/),&
            deBlockList=deblklist,rc=status)
       call LIS_verify(status)

       do i=0,LIS_npes-1
          stid = LIS_goffsets(n,i)+1
          enid = stid + LIS_ngrids(n,i)-1

          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)
       enddo
        
       gridDG = ESMF_DistGridCreate(minIndex=(/1/), &
            maxIndex=(/LIS_rc%glbngrid(n)/),&
            deBlockList=deblklist,rc=status)
       call LIS_verify(status)

       do m=1,LIS_rc%max_model_types

          do i=0,LIS_npes-1
             stid = LIS_patch_offsets(n,m,i)+1
             enid = stid + LIS_npatches(n,m,i)-1
             
             deblklist(:,1,i+1) = (/stid/)
             deblklist(:,2,i+1) = (/enid/)
          enddo

          patchDG(m) = ESMF_DistGridCreate(minIndex=(/1/),&
               maxIndex=(/LIS_rc%glbnpatch(n,m)/),&
               deBlockList=deblklist,rc=status)
          call LIS_verify(status)
       enddo

       deallocate(deblklist)
       
       LIS_vecTile(n) = ESMF_GridCreate(name = "LIS Tile Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = tileDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LIS_verify(status)

       LIS_vecGrid(n) = ESMF_GridCreate(name = "LIS Grid Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = gridDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LIS_verify(status)

       do m=1,LIS_rc%max_model_types
          LIS_vecPatch(n,m) = ESMF_GridCreate(name="LIS Patch Space",&
               coordTypeKind=ESMF_TYPEKIND_R4, distGrid = patchDG(m),&
               gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
          call LIS_verify(status)
       enddo
       
       call LIS_histDataInit(n,LIS_rc%ntiles(n))

    enddo
    TRACE_EXIT("dom_init")

  end subroutine LIS_domain_init

!BOP
! !ROUTINE: LIS_domain_setup
! \label{LIS_domain_setup}
! 
! !INTERFACE: 
  subroutine LIS_domain_setup(n)
! !USES: 

! !ARGUMENTS: 
    integer,  intent(IN) :: n

! 
! !DESCRIPTION: 
!  This routines computes the variables used to map between the tile, grid, 
!  and patch spaces in LIS. The global indices including the halo 
!  regions are also computed. 
!EOP
    integer, allocatable :: ntiles_pergrid(:)
    integer, allocatable :: ntiles_pergrid_red(:)
    integer, allocatable :: npatch_pergrid(:,:)
    integer, allocatable :: npatch_pergrid_red(:,:)
    integer, allocatable :: gtmp(:,:)
    integer, allocatable :: gtmp1(:)
    integer              :: count, index
    integer              :: t, m, gid, gid_red, ierr
    integer              :: c,r, l, c1,c2, r1,r2
    integer              :: ntiles, npatch,stid

    TRACE_ENTER("dom_setup")
    allocate(LIS_domain(n)%ntiles_pergrid(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
    LIS_domain(n)%ntiles_pergrid = 0 ! EMK TEST
    allocate(LIS_domain(n)%str_tind(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
    LIS_domain(n)%str_tind = 0 ! EMK TEST
    allocate(ntiles_pergrid(LIS_rc%lnc(n)*LIS_rc%lnr(n)),stat=ierr)
    ntiles_pergrid = 0
    allocate(ntiles_pergrid_red(LIS_rc%lnc_red(n)*LIS_rc%lnr_red(n)),stat=ierr)
    ntiles_pergrid_red = 0
    allocate(npatch_pergrid(LIS_rc%lnc(n)*LIS_rc%lnr(n),LIS_rc%max_model_types))
    npatch_pergrid = 0
    allocate(npatch_pergrid_red(LIS_rc%lnc(n)*LIS_rc%lnr(n),LIS_rc%max_model_types))
    npatch_pergrid_red = 0
    do m=1,LIS_rc%max_model_types
       allocate(LIS_surface(n,m)%npatch_pergrid(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
       LIS_surface(n,m)%npatch_pergrid = 0
       allocate(LIS_surface(n,m)%str_patch_ind(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
       LIS_surface(n,m)%str_patch_ind = 0
    enddo

    do t=1,LIS_rc%ntiles(n)
       LIS_domain(n)%tile(t)%pens = 1.0/float(LIS_rc%nensem(n))      
    enddo
    do m=1,LIS_rc%max_model_types
       do t=1,LIS_rc%npatch(n,m)
          LIS_surface(n,m)%tile(t)%pens = 1.0/float(LIS_rc%nensem(n))
       enddo
    enddo

    do t=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
       ntiles_pergrid(t) = 0 
    enddo
    do t=1,LIS_rc%ntiles(n)
       index = LIS_domain(n)%gindex(LIS_domain(n)%tile(t)%col,&
            LIS_domain(n)%tile(t)%row)
       if(index.ne.-1) then 
          gid = LIS_domain(n)%tile(t)%col+&
               (LIS_domain(n)%tile(t)%row-1)*LIS_rc%lnc(n)
          ntiles_pergrid(gid) = ntiles_pergrid(gid)+1
       endif
    enddo

    do r=LIS_nss_halo_ind(n,LIS_localPet+1), LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1), LIS_ewe_halo_ind(n,LIS_localPet+1)
          c1 = c-LIS_ews_halo_ind(n,LIS_localPet+1)+1
          r1 = r-LIS_nss_halo_ind(n,LIS_localPet+1)+1
          gid = c1+(r1-1)*LIS_rc%lnc(n)
          
          if(r.ge.LIS_nss_ind(n,LIS_localPet+1).and. &
               r.le.LIS_nse_ind(n,LIS_localPet+1).and.&
               c.ge.LIS_ews_ind(n,LIS_localPet+1).and.&
               c.le.LIS_ewe_ind(n,LIS_localPet+1))then !points not in halo
             c2 = c-LIS_ews_ind(n,LIS_localPet+1)+1
             r2 = r-LIS_nss_ind(n,LIS_localPet+1)+1
             gid_red = c2+(r2-1)*LIS_rc%lnc_red(n)
             ntiles_pergrid_red(gid_red) = ntiles_pergrid(gid)
          endif
       enddo
    enddo


    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(gtmp1(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
    else
       allocate(gtmp1(1))
    endif

#if (defined SPMD)
    call MPI_GATHERV(ntiles_pergrid_red,LIS_deltas(n,LIS_localPet),MPI_INTEGER,gtmp1,&
         LIS_deltas(n,:),LIS_offsets(n,:),MPI_INTEGER,0,LIS_mpi_comm,ierr)
#else
    gtmp1 = ntiles_pergrid_red
#endif

    if(LIS_masterproc) then 
       count=1
       do l=1,LIS_npes
          do r=LIS_nss_ind(n,l), LIS_nse_ind(n,l)
             do c=LIS_ews_ind(n,l), LIS_ewe_ind(n,l)
                gtmp(c,r) = gtmp1(count)
                count = count+1
             enddo
          enddo
       enddo

       count=1
       do r=1,LIS_rc%gnr(n)
          do c=1,LIS_rc%gnc(n)
             LIS_domain(n)%ntiles_pergrid(count) = gtmp(c,r)
             count = count+1
          enddo
       enddo
       
       deallocate(gtmp)
       deallocate(gtmp1)
       
      
       if(LIS_domain(n)%ntiles_pergrid(1).ge.0) then 
          LIS_domain(n)%str_tind(1) = 1
       else
          LIS_domain(n)%str_tind(1) = 0 
       endif
       
       do t=2,LIS_rc%gnc(n)*LIS_rc%gnr(n)
          LIS_domain(n)%str_tind(t) = LIS_domain(n)%str_tind(t-1)+&
               LIS_domain(n)%ntiles_pergrid(t-1)
       enddo       
    else
       deallocate(gtmp1)
    endif

#if (defined SPMD)
    call MPI_Bcast(LIS_domain(n)%ntiles_pergrid,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
         MPI_INTEGER,0,LIS_mpi_comm,ierr)
    call MPI_Bcast(LIS_domain(n)%str_tind,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
         MPI_INTEGER,0,LIS_mpi_comm,ierr)
#endif   

    if(LIS_masterproc) then 
       LIS_rc%glbngrid_red(n)=0
       LIS_rc%glbntiles_red(n)=0
       do l=1,LIS_npes
          do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
             do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                if(r.ge.LIS_nss_ind(n,l).and.&
                     r.le.LIS_nse_ind(n,l).and.&
                     c.ge.LIS_ews_ind(n,l).and.&
                     c.le.LIS_ewe_ind(n,l))then !points not in halo
                   gid = c+(r-1)*LIS_rc%gnc(n)
                   ntiles = LIS_domain(n)%ntiles_pergrid(gid)
                   stid = LIS_domain(n)%str_tind(gid)
                   if ( ntiles .ne. 0 ) then
                      LIS_rc%glbngrid_red(n) = LIS_rc%glbngrid_red(n) + 1
                   endif
                   do t=1,ntiles
                      LIS_rc%glbntiles_red(n) = LIS_rc%glbntiles_red(n) + 1
                   enddo
                endif
             enddo
          enddo
       enddo
    endif    
#if (defined SPMD)
    call MPI_Bcast(LIS_rc%glbntiles_red(n),1,&
         MPI_INTEGER,0,LIS_mpi_comm,ierr)
    call MPI_Bcast(LIS_rc%glbngrid_red(n),1,&
         MPI_INTEGER,0,LIS_mpi_comm,ierr)
#endif

    do t=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
       npatch_pergrid(t,:) = 0
    enddo

    do m=1,LIS_rc%max_model_types
       do t=1,LIS_rc%npatch(n,m)
          index =  LIS_domain(n)%gindex(LIS_surface(n,m)%tile(t)%col,&
               LIS_surface(n,m)%tile(t)%row)
          if(index.ne.-1) then 
             gid = LIS_surface(n,m)%tile(t)%col + & 
                  (LIS_surface(n,m)%tile(t)%row-1)*LIS_rc%lnc(n)
             npatch_pergrid(gid,m) = npatch_pergrid(gid,m) + 1
          endif
       enddo
    enddo

    do r=LIS_nss_halo_ind(n,LIS_localPet+1), LIS_nse_halo_ind(n,LIS_localPet+1)
       do c=LIS_ews_halo_ind(n,LIS_localPet+1), LIS_ewe_halo_ind(n,LIS_localPet+1)
          c1 = c-LIS_ews_halo_ind(n,LIS_localPet+1)+1
          r1 = r-LIS_nss_halo_ind(n,LIS_localPet+1)+1
          gid = c1+(r1-1)*LIS_rc%lnc(n)
          
          if(r.ge.LIS_nss_ind(n,LIS_localPet+1).and. &
               r.le.LIS_nse_ind(n,LIS_localPet+1).and.&
               c.ge.LIS_ews_ind(n,LIS_localPet+1).and.&
               c.le.LIS_ewe_ind(n,LIS_localPet+1))then !points not in halo
             c2 = c-LIS_ews_ind(n,LIS_localPet+1)+1
             r2 = r-LIS_nss_ind(n,LIS_localPet+1)+1
             gid_red = c2+(r2-1)*LIS_rc%lnc_red(n)
             npatch_pergrid_red(gid_red,:) = npatch_pergrid(gid,:)
          endif
       enddo
    enddo

    do m=1,LIS_rc%max_model_types
       if(LIS_masterproc) then 
          allocate(gtmp(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(gtmp1(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
       else
          allocate(gtmp1(1))
       endif
       
#if (defined SPMD)
       call MPI_GATHERV(npatch_pergrid_red(:,m),LIS_deltas(n,LIS_localPet),&
            MPI_INTEGER,gtmp1,&
            LIS_deltas(n,:),LIS_offsets(n,:),MPI_INTEGER,0,LIS_mpi_comm,ierr)
#else
       gtmp1 = npatch_pergrid_red(:,m)
#endif

       if(LIS_masterproc) then 
          count=1
          do l=1,LIS_npes
             do r=LIS_nss_ind(n,l), LIS_nse_ind(n,l)
                do c=LIS_ews_ind(n,l), LIS_ewe_ind(n,l)
                   gtmp(c,r) = gtmp1(count)
                   count = count+1
                enddo
             enddo
          enddo
          
          count=1
          do r=1,LIS_rc%gnr(n)
             do c=1,LIS_rc%gnc(n)
                LIS_surface(n,m)%npatch_pergrid(count) = gtmp(c,r)
                count = count+1
             enddo
          enddo
          
          deallocate(gtmp)
          deallocate(gtmp1)
       
          
          if(LIS_surface(n,m)%npatch_pergrid(1).ge.0) then 
             LIS_surface(n,m)%str_patch_ind(1) = 1
          else
             LIS_surface(n,m)%str_patch_ind(1) = 0 
          endif
          
          do t=2,LIS_rc%gnc(n)*LIS_rc%gnr(n)
             LIS_surface(n,m)%str_patch_ind(t) = LIS_surface(n,m)%str_patch_ind(t-1)+&
                  LIS_surface(n,m)%npatch_pergrid(t-1)
          enddo
       else
          deallocate(gtmp1)
       endif
#if (defined SPMD)
       call MPI_Bcast(LIS_surface(n,m)%npatch_pergrid,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
            MPI_INTEGER,0,LIS_mpi_comm,ierr)
       call MPI_Bcast(LIS_surface(n,m)%str_patch_ind,LIS_rc%gnc(n)*LIS_rc%gnr(n),&
            MPI_INTEGER,0,LIS_mpi_comm,ierr)
#endif   

       if(LIS_masterproc) then 
          LIS_rc%glbnpatch_red(n,m)=0
          do l=1,LIS_npes
             do r=LIS_nss_halo_ind(n,l),LIS_nse_halo_ind(n,l)
                do c=LIS_ews_halo_ind(n,l),LIS_ewe_halo_ind(n,l)
                   if(r.ge.LIS_nss_ind(n,l).and.&
                        r.le.LIS_nse_ind(n,l).and.&
                        c.ge.LIS_ews_ind(n,l).and.&
                        c.le.LIS_ewe_ind(n,l))then !points not in halo
                      gid = c+(r-1)*LIS_rc%gnc(n)
                      npatch = LIS_surface(n,m)%npatch_pergrid(gid)
                      stid = LIS_surface(n,m)%str_patch_ind(gid)
                      do t=1,npatch
                         LIS_rc%glbnpatch_red(n,m) = &
                              LIS_rc%glbnpatch_red(n,m) + 1
                      enddo
                   endif
                enddo
             enddo
          enddo
       endif
#if (defined SPMD)
       call MPI_Bcast(LIS_rc%glbnpatch_red(n,m),1,&
            MPI_INTEGER,0,LIS_mpi_comm,ierr)
#endif

    enddo
    deallocate(ntiles_pergrid)
    deallocate(ntiles_pergrid_red)
    deallocate(npatch_pergrid)
    deallocate(npatch_pergrid_red)
    TRACE_EXIT("dom_setup")

  end subroutine LIS_domain_setup

!BOP
! !ROUTINE: LIS_quilt_domain
! \label{LIS_quilt_domain}
! 
! !INTERFACE: 
subroutine LIS_quilt_domain(n, nc, nr )
! !USES:     

! !ARGUMENTS: 
   integer, intent(in) :: n
   integer, intent(in) :: nc
   integer, intent(in) :: nr
! 
! !DESCRIPTION: 
! This routine generates the quilted domain extents and sizes based on the 
! processor layout and the size of the specified halos. 
! 
!EOP

   integer :: ips, ipe, jps, jpe
   integer :: i, j,ierr
   integer :: deltas
   integer :: ews_halo_ind
   integer :: ewe_halo_ind
   integer :: nss_halo_ind
   integer :: nse_halo_ind
   integer :: ews_ind
   integer :: ewe_ind
   integer :: nss_ind
   integer :: nse_ind

   TRACE_ENTER("dom_quilt")
   if ( LIS_rc%decompose_by_processes ) then
      call decompose_npes(n, nc, nr, ips, ipe, jps, jpe)
   else
      call decompose_nx_ny(nc, nr, ips, ipe, jps, jpe)
   endif

   ews_halo_ind = max(ips-LIS_rc%halox, 1)
   ewe_halo_ind = min(ipe+LIS_rc%halox, nc)
   nss_halo_ind = max(jps-LIS_rc%haloy, 1)
   nse_halo_ind = min(jpe+LIS_rc%haloy, nr)

   ews_ind = ips
   ewe_ind = min(ipe, nc)
   nss_ind = jps
   nse_ind = min(jpe, nr)

   LIS_rc%lnc(n) = ewe_halo_ind - ews_halo_ind + 1
   LIS_rc%lnr(n) = nse_halo_ind - nss_halo_ind + 1

   LIS_rc%lnc_red(n)= ewe_ind - ews_ind + 1
   LIS_rc%lnr_red(n)= nse_ind - nss_ind + 1

   write(unit=LIS_logunit,fmt=*)'[INFO] local domain',':(', &
                                LIS_rc%lnc(n),LIS_rc%lnr(n),')'
   write(unit=LIS_logunit,fmt=*)'[INFO] local domain without halo',':(', &
                                LIS_rc%lnc_red(n),LIS_rc%lnr_red(n),')'
!   write(unit=LIS_logunit,fmt=*)'parameter domain',':(', &
!                                LIS_rc%pnc(n),LIS_rc%pnr(n),')'
   write(unit=LIS_logunit,fmt=*)'[INFO] running domain',':(',nc,nr,')'

   deltas = LIS_rc%lnc_red(n)*LIS_rc%lnr_red(n)

#if (defined SPMD)
   call MPI_ALLGATHER(deltas,1,MPI_INTEGER,LIS_deltas(n,:),1,MPI_INTEGER,   &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(ews_ind,1,MPI_INTEGER,LIS_ews_ind(n,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nss_ind,1,MPI_INTEGER,LIS_nss_ind(n,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(ewe_ind,1,MPI_INTEGER,LIS_ewe_ind(n,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nse_ind,1,MPI_INTEGER,LIS_nse_ind(n,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)

   call MPI_ALLGATHER(ews_halo_ind,1,MPI_INTEGER,LIS_ews_halo_ind(n,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nss_halo_ind,1,MPI_INTEGER,LIS_nss_halo_ind(n,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(ewe_halo_ind,1,MPI_INTEGER,LIS_ewe_halo_ind(n,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nse_halo_ind,1,MPI_INTEGER,LIS_nse_halo_ind(n,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
#else
   LIS_deltas(n,:)  = deltas
   LIS_ews_ind(n,:) = ews_ind
   LIS_nss_ind(n,:) = nss_ind
   LIS_ewe_ind(n,:) = ewe_ind
   LIS_nse_ind(n,:) = nse_ind

   LIS_ews_halo_ind(n,:) = ews_halo_ind
   LIS_nss_halo_ind(n,:) = nss_halo_ind
   LIS_ewe_halo_ind(n,:) = ewe_halo_ind
   LIS_nse_halo_ind(n,:) = nse_halo_ind
#endif   
   if(LIS_masterproc) then
      LIS_offsets(n,0) = 0
      do i=1,LIS_npes-1
         LIS_offsets(n,i) = LIS_offsets(n,i-1)+LIS_deltas(n,i-1)
      enddo
   end if
#if (defined SPMD)
   call MPI_BCAST(LIS_offsets(n,:), LIS_npes, MPI_INTEGER,0, LIS_mpi_comm, ierr)
#endif
TRACE_EXIT("dom_quilt")
end subroutine LIS_quilt_domain


!BOP
! !ROUTINE: LIS_quilt_b_domain
! \label{LIS_quilt_b_domain}
! 
! !INTERFACE: 
subroutine LIS_quilt_b_domain(n, nc, nr, nc_raw, nr_raw)
! !USES:     

! !ARGUMENTS: 
   integer, intent(in)  :: n
   integer, intent(in)  :: nc
   integer, intent(in)  :: nr
   integer, intent(in)  :: nc_raw
   integer, intent(in)  :: nr_raw
! 
! !DESCRIPTION: 
! This routine generates the quilted domain extents and sizes based on the 
! processor layout and the size of the specified halos. 
! 
!EOP

   integer :: ips, ipe, jps, jpe
   integer :: Px, Py, P
   integer :: mytask_x, mytask_y
   integer :: i, j,ierr
   integer :: deltas
   integer :: ews_b_ind
   integer :: ewe_b_ind
   integer :: nss_b_ind
   integer :: nse_b_ind
   integer :: b_nc
   integer :: b_nr

   if ( LIS_rc%decompose_by_processes ) then
      call decompose_npes(n, nc_raw, nr_raw, ips, ipe, jps, jpe)
   else
      call decompose_nx_ny(nc_raw, nr_raw, ips, ipe, jps, jpe)
   endif

   b_nc = (nc - nc_raw)/2
   b_nr = (nr - nr_raw)/2

   ews_b_ind = max(ips-LIS_rc%halox, 1) - b_nc
   ewe_b_ind = min(ipe+LIS_rc%halox, nc_raw) + b_nc
   nss_b_ind = max(jps-LIS_rc%haloy, 1) - b_nr
   nse_b_ind = min(jpe+LIS_rc%haloy, nr_raw) + b_nr

   LIS_rc%lnc_b(n) = ewe_b_ind-ews_b_ind+1
   LIS_rc%lnr_b(n) = nse_b_ind-nss_b_ind+1

#if (defined SPMD)
   call MPI_ALLGATHER(ews_b_ind,1,MPI_INTEGER,&
        LIS_ews_b_ind(n,:),1,MPI_INTEGER,&
        LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nss_b_ind,1,MPI_INTEGER,&
        LIS_nss_b_ind(n,:),1,MPI_INTEGER,&
        LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(ewe_b_ind,1,MPI_INTEGER,&
        LIS_ewe_b_ind(n,:),1,MPI_INTEGER,&
        LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nse_b_ind,1,MPI_INTEGER,&
        LIS_nse_b_ind(n,:),1,MPI_INTEGER,&
        LIS_mpi_comm,ierr)
#else
   LIS_ews_b_ind(n,:) = ews_b_ind
   LIS_nss_b_ind(n,:) = nss_b_ind
   LIS_ewe_b_ind(n,:) = ewe_b_ind
   LIS_nse_b_ind(n,:) = nse_b_ind
#endif   
end subroutine LIS_quilt_b_domain


!BOP
! !ROUTINE: LIS_domain_finalize
! \label{LIS_domain_finalize}
!
! !INTERFACE: 
  subroutine LIS_domain_finalize
!
! !DESCRIPTION: 
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures in the specific instance of the 
!  domain implementation. 
! 
!EOP
  end subroutine LIS_domain_finalize

!BOP
! 
! !ROUTINE: make_domain
! \label{make_domain}
!
! !INTERFACE: 
  subroutine make_domain()
!
! !DESCRIPTION: 
!  This subroutine generates the tile, patch and grid data structures
!  based on the vegetation, soils, and topography datasets. The
!  routine first computes the dominant distributions for vegetation
!  and soils and then generates the multi-dimensional tile structure.  
!
!EOP
    implicit none

    integer :: c,r
    integer :: n
    logical :: soilt_selected
    logical :: soilf_selected
    logical :: elev_selected
    logical :: slope_selected
    logical :: aspect_selected
    logical :: curvature_selected
    
    do n=1,LIS_rc%nnest
       if(LIS_rc%usetexturemap(n).ne."none") then 
          soilt_selected = .true. 
       else
          soilt_selected = .false.
       endif      

       if(LIS_rc%usesoilfractionmap(n).ne."none") then 
          soilf_selected = .true. 
       else
          soilf_selected = .false.
       endif
       
       if(soilt_selected.and.soilf_selected) then 
          soilt_selected = .true. 
          soilf_selected = .false. 
       endif

       if(LIS_rc%useelevationmap(n).ne."none") then 
          elev_selected = .true. 
       else
          elev_selected = .false. 
       endif

       if(LIS_rc%useslopemap(n).ne."none") then 
          slope_selected = .true. 
       else
          slope_selected = .false. 
       endif

       if(LIS_rc%useaspectmap(n).ne."none") then 
          aspect_selected = .true. 
       else
          aspect_selected = .false. 
       endif

       if(LIS_rc%usecurvaturemap(n).ne."none") then
          curvature_selected = .true.
       else
          curvature_selected = .false.
       endif

!-----------------------------------------------------------------------
! normalize the parameter data distributions
!-----------------------------------------------------------------------
       ! EMK...Special handling of single surface tile case.  First find
       ! dominant surface model type in each grid point (LSM, Openwater, etc),
       ! then find the dominant tile for that model type.
       if (LIS_rc%surface_maxt == 1) then          
          call calculate_dominant_sfctile(n, LIS_LMLC(n)%dommask, &
               LIS_rc%nsurfacetypes, &
               LIS_rc%surface_minp, LIS_rc%surface_maxt, &
               LIS_LMLC(n)%surfacetype, LIS_LMLC(n)%landcover)
       else
          call calculate_domdistribution(n, LIS_rc%nsurfacetypes, &
               LIS_rc%surface_minp, LIS_rc%surface_maxt, LIS_LMLC(n)%landcover)
       end if
       if(soilt_selected) then 
          call calculate_domdistribution(n, LIS_rc%nsoiltypes,&
               LIS_rc%soilt_minp, LIS_rc%soilt_maxt, LIS_soils(n)%texture)
       endif

       if(soilf_selected) then 
          call calculate_domdistribution(n, LIS_rc%nsoilfbands,&
               LIS_rc%soilf_minp, LIS_rc%soilf_maxt, LIS_soils(n)%soilffgrd)
       endif

       if(elev_selected) then 
          call calculate_domdistribution(n, LIS_rc%nelevbands, &
               LIS_rc%elev_minp, LIS_rc%elev_maxt, LIS_topo(n)%elevfgrd)
       endif

       if(slope_selected) then 
          call calculate_domdistribution(n, LIS_rc%nslopebands, &
               LIS_rc%slope_minp, LIS_rc%slope_maxt, LIS_topo(n)%slopefgrd)
       endif

       if(aspect_selected) then 
          call calculate_domdistribution(n, LIS_rc%naspectbands, &
               LIS_rc%aspect_minp, LIS_rc%aspect_maxt, LIS_topo(n)%aspectfgrd)
       endif

       call create_tilespace(n, soilt_selected, soilf_selected, elev_selected, &
            slope_selected, aspect_selected, curvature_selected)
    enddo

  end subroutine make_domain

  ! EMK...Special subroutine for finding single dominant surface tile in each
  ! grid point. This first considers the single dominant surface model type 
  ! (LSM, Openwater, etc), and then identifies the dominant tile for that 
  ! model type.
  !
  ! The multiple surface tile case is more complicated and is not addressed
  ! here.  
  subroutine calculate_dominant_sfctile(n, dommask, &
       ntypes, minp, maxt, surfacetypes, fgrd)

     ! Defaults
     implicit none
     
     ! Arguments
     integer, intent(in)  :: n
     real, intent(in)     :: dommask(LIS_rc%lnc(n), LIS_rc%lnr(n))
     integer, intent(in)  :: ntypes
     real, intent(in)     :: minp
     integer, intent(in)  :: maxt
     real, intent(in)     :: surfacetypes(LIS_rc%lnc(n), LIS_rc%lnr(n), ntypes)
     real, intent(inout)  :: fgrd(LIS_rc%lnc(n), LIS_rc%lnr(n), ntypes)
     
     ! Locals
     integer :: c, r, t, m, mm, tt
     real    :: maxv, maxtilefrac
     real    :: model_type_fractions(LIS_rc%max_model_types)

     ! Sanity check
     if (maxt > 1) then
        write(LIS_logunit,*) &
             '[ERR] Calculate_dominant_sfctile only works for maxt = 1!'
        call LIS_endrun()
     end if
     
     ! Loop through each patch, figure out the dominant surface model type,
     ! and then find the dominant land cover category for that model type.
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)

           ! Skip points that are outside of the domain--these have no
           ! land cover values defined.
           if (dommask(c,r) .le. 0) cycle

           ! Consolidate the tile fractions into fractions by surface model 
           ! type (LSM, Openwater, etc)
           model_type_fractions(:) = 0
           do t=1,ntypes
              m = nint(surfacetypes(c,r,t))
              if (m .eq. 0) cycle ! Empty tile
              model_type_fractions(m) = &
                   model_type_fractions(m) + fgrd(c,r,t)
           end do ! t

           ! Now find which surface model type has the highest fraction.
           maxv = 0
           mm = 0
           do m = 1, LIS_rc%max_model_types
              ! Skip model type if not used
              if (.not. isSurfaceTypeSelected(LIS_rc%sf_model_type(m))) cycle
              if (model_type_fractions(m) > maxv) then
                 maxv = model_type_fractions(m)
                 mm = m
              end if
           end do ! m
           
           ! Sanity check
           if (mm .eq. 0) then
              write(LIS_logunit,*) &
                   '[ERR] No dominant surface model type found!'
              write(LIS_logunit,*) &
                   'c,r,fgrd: ',c,r,fgrd(c,r,:)
              flush(LIS_logunit)
              call LIS_endrun()
           end if

           ! We now know which surface model type is dominant.
           ! Identify and exclude tiles not used by that model type.
           do t = 1,ntypes
              if (nint(surfacetypes(c,r,t)) .eq. mm) cycle
              fgrd(c,r,t) = 0.0
           end do ! t

           ! At this point, remaining tiles are associated with only a single 
           ! surface model type. Exclude tiles below minimum tile grid area). 
           maxtilefrac=maxval(fgrd(c,r,:))
           do t=1,ntypes
              if (fgrd(c,r,t).lt.minp) then
                 fgrd(c,r,t)=0.0 
              endif
           enddo ! t

           ! Now find the tile with the largest fraction.
           maxv = 0
           tt = 0
           do t=1,ntypes
              if (fgrd(c,r,t) > maxv) then
                 maxv = fgrd(c,r,t)
                 tt = t
              endif
           enddo ! t
           
           ! Sanity check: Make sure we still have one tile left!
           if (tt .eq. 0) then
              write(LIS_logunit,*) &
                "[ERR] No surface model tile fraction >= minp,",minp
              write(LIS_logunit,*) &
                "  which is the 'Minimum cutoff percentage (surface type tiles):'"
              write(LIS_logunit,*) &
                "  value set in your lis.config file. Highest tile fraction value is:",maxtilefrac
              write(LIS_logunit,*) &
                "  Thus, surface tiles set to '0' for gridpoint:"
              write(LIS_logunit,*) &
                   'c,r,fgrd: ',c,r,fgrd(c,r,:)
              flush(LIS_logunit)
              call LIS_endrun()
           end if

           ! Exclude the smaller tiles, and set the dominant tile to 100%
           do t = 1,ntypes
              if (t .eq. tt) then
                 fgrd(c,r,t) = 1.0
              else
                 fgrd(c,r,t) = 0.0
              end if
           end do ! t

        end do ! c
     end do ! r
     
  end subroutine calculate_dominant_sfctile

!BOP
!
! !ROUTINE: calculate_domdistribution
! \label{calculate_domdistribution}
!
! !REVISION HISTORY:
!  09 Sept 2004: Sujay Kumar ; Initial version
!
! !INTERFACE:
  subroutine calculate_domdistribution(n, ntypes, minp, maxt, fgrd)
! !USES:
    implicit none

! !ARGUMENTS: 
    integer, intent(in)  :: n
    integer              :: ntypes
    real                 :: minp
    integer              :: maxt
    real                 :: fgrd(LIS_rc%lnc(n), LIS_rc%lnr(n), ntypes)
!
! !DESCRIPTION:
!  This routine determines the percentages of 
!  dominant types for subgrid tiling.
!  The routine uses the distribution of data within a 
!  grid cell and imposes the specified cutoff percentage and 
!  the maximum number of data types allowed to create the 
!  number of tiles. 
!
!  The arguments are: 
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ntypes]
!    number of categories or bins in the input data. 
!   \item[minp]
!    minimum cutoff percentage to be used in the renormalization
!   \item[maxt]
!    maximum number of tiles per grid cell to be used in the renormalization
!   \item[fgrd]
!    fraction of grid covered by a data type. 
!  \end{description}
!EOP

    real,    allocatable :: tsum(:,:)  !Temporary processing variable
    integer, allocatable :: pfrac(:,:,:)
    integer :: c, r, t, ierr, i, j
    real    :: rsum         
    real    :: fvt(ntypes)  
    real    :: maxv

    allocate(tsum(LIS_rc%lnc(n), LIS_rc%lnr(n)), stat=ierr)
    call LIS_verify(ierr,'Error allocating tsum.')
    tsum = 0.0
!----------------------------------------------------------------------      
! Exclude tiles with (minimum tile grid area),  
! normalize remaining tiles to 100%
!----------------------------------------------------------------------      
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n)            
          rsum=0.0
          do t=1,ntypes
             if(fgrd(c,r,t).lt.minp)then
                fgrd(c,r,t)=0.0 
             endif
             rsum=rsum+fgrd(c,r,t)
          enddo
!----------------------------------------------------------------------      
! renormalize veg fractions within a grid to 1
!----------------------------------------------------------------------      
          if(rsum.gt.0.0) then  
             do t=1,ntypes
                fgrd(c,r,t)=fgrd(c,r,t)/rsum
             enddo
             
             rsum=0.0
             do t=1,ntypes
                rsum=rsum+fgrd(c,r,t)
             enddo
             
             if(rsum.lt.0.9999.or.rsum.gt.1.0001)then 
                write(LIS_logunit,*) '[ERR] The tile distribution do not sum to 100%'
                call LIS_endrun()
             endif
          endif
       enddo
    enddo

    allocate(pfrac(LIS_rc%lnc(n),LIS_rc%lnr(n),ntypes))
!----------------------------------------------------------------------      
! Exclude tiles with SURFACE_MAXT (Maximum Tiles per grid), 
!   normalize remaining tiles to 100%
! Determine the grid predominance order of the tiles
!  PFRAC(NT) will contain the predominance order of tiles
!----------------------------------------------------------------------      
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n) 
          do t=1,ntypes
             fvt(t)=fgrd(c,r,t)
             pfrac(c,r,t)=0
          enddo
          do i=1,ntypes
             maxv=0.0
             t=0
             do j=1,ntypes
                if(fvt(j).gt.maxv)then
                   if(fgrd(c,r,j).gt.0) then
                      maxv=fvt(j)
                      t=j
                   endif
                endif
             enddo
             if(t.gt.0) then
                pfrac(c,r,t)=i
                fvt(t)=-999.0       
             endif
          enddo
       enddo
    enddo
!----------------------------------------------------------------------      
! Impose MAXT Cutoff
!----------------------------------------------------------------------
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n)
          rsum=0.0
          do t=1,ntypes
             if(pfrac(c,r,t).lt.1) then
                fgrd(c,r,t)=0.0    
                pfrac(c,r,t)=0  
             endif
             if(pfrac(c,r,t).gt.maxt) then
                fgrd(c,r,t)=0.0            
                pfrac(c,r,t)=0  
             endif
             rsum=rsum+fgrd(c,r,t)
          enddo
!----------------------------------------------------------------------
! renormalize veg fractions within a grid to 1
!----------------------------------------------------------------------
          if(rsum.gt.0.0) then  
             do t=1,ntypes
                fgrd(c,r,t)= fgrd(c,r,t)/rsum
             enddo
             
             rsum=0.0
             do t=1,ntypes
                rsum=rsum+ fgrd(c,r,t)  !recalculate rsum to check 
             enddo
             tsum(c,r)=rsum
           
             if(rsum.lt.0.9999.or.rsum.gt.1.0001)then  !check renormalization
                write(LIS_logunit,*) &
                     '[ERR] Renormalization failed in calculate_domdistribution'
                print*, c,r,fgrd(c,r,:)
                call LIS_endrun()
             endif
          endif
       enddo
    enddo
    deallocate(pfrac)
    deallocate(tsum)

  end subroutine calculate_domdistribution


!BOP
!
! !ROUTINE: create_tilespace
!  \label{create_tilespace}
!
! !REVISION HISTORY:
!  09 Sept 2004: Sujay Kumar ; Initial version
!
! !INTERFACE:
  subroutine create_tilespace(n,soilt_selected, soilf_selected, elev_selected,&
       slope_selected, aspect_selected, curvature_selected)

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: n
    logical, intent(in)   :: soilt_selected
    logical, intent(in)   :: soilf_selected
    logical, intent(in)   :: elev_selected
    logical, intent(in)   :: slope_selected
    logical, intent(in)   :: aspect_selected
    logical, intent(in)   :: curvature_selected

! !DESCRIPTION: 
!  This routine creates the tilespace based on the specified dimensions
!  of tiling. Vegetation is used by default. Soil texture, soil fractions, 
!  elevation, slope, aspect can be specified as additional dimensions for
!  tiling. The routine also generates the patch data structures which are
!  the tile data structures in each surface model class. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[soilt\_selected]
!    flag to indicate if soil texture based tiling is used
!   \item[soilf\_selected]
!    flag to indicate if soil fraction based tiling is used
!   \item[elev\_selected]
!    flag to indicate if elevation based tiling is used
!   \item[slope\_selected]
!    flag to indicate if slope based tiling is used
!   \item[aspect\_selected]
!    flag to indicate if aspect based tiling is used
!   \item[curvature\_selected]
!    flag to indicate if curvature based tiling is used
!  \end{description}
!EOP
    real          :: locallat, locallon
    integer       :: c, r, t, m
    integer       :: iv, it, ie, is, ia, ic
    integer       :: vegt
    integer       :: soilt
    real          :: sand, clay, silt
    real          :: elev, slope, aspect, curvature
    integer       :: sf_index
    integer       :: gnc, gnr
    integer       :: kk
    integer       :: kk_sf(LIS_rc%max_model_types)
    real, allocatable :: sumv(:)
    integer       :: ntiles_surface
    integer       :: npatch_surface
    integer       :: ntiles_soilt
    integer       :: ntiles_soilf
    integer       :: ntiles_soil
    integer       :: ntiles_elev
    integer       :: ntiles_slope
    integer       :: ntiles_aspect
    integer       :: ntiles_curvature
    integer       :: ntiles_land_landmask
    integer       :: gid
    real          :: temp
    integer       :: soilf_index
    integer       :: elev_index
    integer       :: slope_index
    integer       :: aspect_index
    integer       :: curvature_index
    

    gnc = LIS_rc%lnc(n)
    gnr = LIS_rc%lnr(n)
    
    LIS_rc%ntiles(n)=0
    LIS_rc%npatch(n,:) = 0 
    do r=1,LIS_rc%lnr(n)     
       do c=1,LIS_rc%lnc(n)             
          do m=1,LIS_rc%nensem(n)             
             if(LIS_LMLC(n)%dommask(c,r).gt.0) then

                ntiles_surface = compute_ntiles_surface(n,c,r)
                ntiles_soilt = compute_ntiles_soilt(n,c,r,soilt_selected)
                ntiles_soilf = compute_ntiles_soilf(n,c,r,soilf_selected)
                ntiles_elev = compute_ntiles_elev(n,c,r,elev_selected)
                ntiles_slope = compute_ntiles_slope(n,c,r,slope_selected)
                ntiles_aspect = compute_ntiles_aspect(n,c,r,aspect_selected)
                ntiles_curvature = compute_ntiles_curvature(n,c,r,curvature_selected)

                if(soilf_selected) then 
                   LIS_rc%ntiles(n) = LIS_rc%ntiles(n)+&
                        ntiles_surface* & 
                        ntiles_soilf* &
                        ntiles_elev*&
                        ntiles_slope*& 
                        ntiles_aspect
                else
                   LIS_rc%ntiles(n) = LIS_rc%ntiles(n)+&
                        ntiles_surface* & 
                        ntiles_soilt* &
                        ntiles_elev*&
                        ntiles_slope*& 
                        ntiles_aspect
                endif
                   
                do t=1,LIS_rc%nsurfacetypes
                   sf_index = nint(LIS_LMLC(n)%surfacetype(c,r,t))
                   npatch_surface = compute_npatches_surface(n,c,r,t,sf_index)
                   if(sf_index.ne.0) then
                      if(soilf_selected) then 
                         LIS_rc%npatch(n,sf_index) = & 
                              LIS_rc%npatch(n,sf_index) + &
                              npatch_surface*& 
                              ntiles_soilf* & 
                              ntiles_elev*& 
                              ntiles_slope*& 
                              ntiles_aspect                         
                      else
                         LIS_rc%npatch(n,sf_index) = & 
                              LIS_rc%npatch(n,sf_index) + &
                              npatch_surface*& 
                              ntiles_soilt* & 
                              ntiles_elev*& 
                              ntiles_slope*& 
                              ntiles_aspect
                      endif
                   endif
                enddo

             endif
          enddo
       enddo
    enddo

    ntiles_land_landmask = 0
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          if (LIS_LMLC(n)%landmask(c,r) .eq. 1) then
             ntiles_land_landmask = ntiles_land_landmask + 1
          end if
       end do
    end do
    write(LIS_logunit,*)'[INFO] Landmask LAND points = ', ntiles_land_landmask

    write(unit=LIS_logunit,fmt=*) &
         '[INFO] Total number of surface tiles',LIS_rc%ntiles(n)
    do m=1,LIS_rc%max_model_types
       if(isSurfaceTypeSelected(LIS_rc%sf_model_type(m))) then 
          write(unit=LIS_logunit,fmt=*) &
               '[INFO] Total number of tiles for ',&
               trim(LIS_rc%sf_model_type_name(m)),&
               LIS_rc%npatch(n,m)
       endif
    enddo
    LIS_rc%ngrid(n)=0
    
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n)
          if(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
               isValidSurfacePoint(LIS_LMLC(n)%landcover(c,r,:),&
               LIS_LMLC(n)%surfacetype(c,r,:))) then 
             LIS_rc%ngrid(n) = LIS_rc%ngrid(n)+1             
          elseif(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
               isValidSurfacePoint(LIS_LMLC(n)%landcover(c,r,:),&
               LIS_LMLC(n)%surfacetype(c,r,:))) then 
             write(LIS_logunit,*) '[ERR] The distribution of surface types in '
             write(LIS_logunit,*) '[ERR] the LANDCOVER field does not sum to 100%'
             write(LIS_logunit,*) '[ERR] mask ',LIS_LMLC(n)%dommask(c,r)
             write(LIS_logunit,*) '[ERR] veg ',LIS_LMLC(n)%landcover(c,r,:)
             call LIS_endrun()
          endif
       enddo
    enddo
!old implementation 
#if 0 
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n)
          if(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
               abs(sum(LIS_LMLC(n)%landcover(c,r,:))-1.0).lt.0.001.and.&
               ifContainsSelectedSurfaceType(&
               nint(LIS_LMLC(n)%surfacetype(c,r,:)))) then 
             LIS_rc%ngrid(n) = LIS_rc%ngrid(n)+1             
          elseif(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
               sum(LIS_LMLC(n)%landcover(c,r,:)).ne.0.and.&
               ifContainsSelectedSurfaceType(&
               nint(LIS_LMLC(n)%surfacetype(c,r,:)))) then 
             write(LIS_logunit,*) '[ERR] The distribution of surface types in '
             write(LIS_logunit,*) '[ERR] the LANDCOVER field does not sum to 100%'
             write(LIS_logunit,*) '[ERR] mask ',LIS_LMLC(n)%dommask(c,r)
             write(LIS_logunit,*) '[ERR] veg ',LIS_LMLC(n)%landcover(c,r,:)
             call LIS_endrun()
          endif
       enddo
    enddo
#endif

    write(LIS_logunit,*) '[INFO] Total number of grid points',LIS_rc%ngrid(n)


    allocate(LIS_domain(n)%tile(LIS_rc%ntiles(n)))
    do m=1,LIS_rc%max_model_types
       allocate(LIS_surface(n,m)%tile(LIS_rc%npatch(n,m)))
    enddo
    allocate(LIS_domain(n)%grid(LIS_rc%ngrid(n)))
    allocate(LIS_domain(n)%gindex(LIS_rc%lnc(n), LIS_rc%lnr(n)))
    
    LIS_domain(n)%minLat = 200.0
    LIS_domain(n)%minLon = 200.00
    LIS_domain(n)%maxLat = -200.0
    LIS_domain(n)%maxLon = -200.00
    
    kk = 1
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n)
          LIS_domain(n)%gindex(c,r) = -1
          
!          locallat = LIS_rc%gridDesc(n,4)+(r-1)*LIS_rc%gridDesc(n,10)
!          locallon = LIS_rc%gridDesc(n,5)+(c-1)*LIS_rc%gridDesc(n,9)

!          call ij_to_latlon(LIS_domain(n)%lisproj,float(c),float(r),&
!               locallat,locallon)
          locallat = LIS_domain(n)%lat(c+(r-1)*LIS_rc%lnc(n))
          locallon = LIS_domain(n)%lon(c+(r-1)*LIS_rc%lnc(n))

          if(localLat.lt.LIS_domain(n)%minLat) LIS_domain(n)%minLat = localLat
          if(localLat.gt.LIS_domain(n)%maxLat) LIS_domain(n)%maxLat = localLat
          if(localLon.lt.LIS_domain(n)%minLon) LIS_domain(n)%minLon = localLon
          if(localLon.gt.LIS_domain(n)%maxLon) LIS_domain(n)%maxLon = localLon

          if(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
               isValidSurfacePoint(LIS_LMLC(n)%landcover(c,r,:),&
               LIS_LMLC(n)%surfacetype(c,r,:))) then 
             LIS_domain(n)%grid(kk)%lat = locallat
             LIS_domain(n)%grid(kk)%lon = locallon
             LIS_domain(n)%grid(kk)%col = c
             LIS_domain(n)%grid(kk)%row = r
             LIS_domain(n)%gindex(c,r) = kk
             kk = kk+1
          endif
!old implementation
#if 0 
          if(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
               abs(sum(LIS_LMLC(n)%landcover(c,r,:))-1.0).lt.0.001.and.&
               ifContainsSelectedSurfaceType(&
               nint(LIS_LMLC(n)%surfacetype(c,r,:)))) then 
             LIS_domain(n)%grid(kk)%lat = locallat
             LIS_domain(n)%grid(kk)%lon = locallon
             LIS_domain(n)%grid(kk)%col = c
             LIS_domain(n)%grid(kk)%row = r
             LIS_domain(n)%gindex(c,r) = kk
             kk = kk+1
          endif
#endif

       enddo
    enddo
    kk = 0
    kk_sf = 0 
    do r=1,LIS_rc%lnr(n)
       do c=1,LIS_rc%lnc(n)
          do m=1,LIS_rc%nensem(n)
             if(LIS_LMLC(n)%dommask(c,r).gt.0) then
                
                ntiles_surface = compute_ntiles_surface(n,c,r)
                ntiles_soilt = compute_ntiles_soilt(n,c,r,soilt_selected)
                ntiles_soilf = compute_ntiles_soilf(n,c,r,soilf_selected)
                ntiles_elev = compute_ntiles_elev(n,c,r,elev_selected)
                ntiles_slope = compute_ntiles_slope(n,c,r,slope_selected)
                ntiles_aspect = compute_ntiles_aspect(n,c,r,aspect_selected)
                ntiles_curvature = compute_ntiles_curvature(n,c,r,curvature_selected)

                if(soilf_selected) then 
                   ntiles_soil = ntiles_soilf
                else
                   ntiles_soil = ntiles_soilt
                endif

                do iv = 1, ntiles_surface
                   call get_vegt_value(n,c,r,iv,vegt)
                   call get_surface_value(n,c,r,iv,sf_index)
                   do it = 1, ntiles_soil

                      sand = -1.0
                      clay = -1.0
                      silt = -1.0
                      soilt = -1.0

                      if(soilf_selected) then 
                         call get_soilf_value(n,c,r,it,soilf_selected,&
                              sand,clay,silt,soilf_index)
                      else
                         call get_soilt_value(n,c,r,it,soilt_selected,soilt)
                      endif

                      do ie = 1, ntiles_elev
                         call get_elev_value(n,c,r,ie,elev_selected,&
                              elev, elev_index)
                         do is=1, ntiles_slope
                            call get_slope_value(n,c,r,is,slope_selected,&
                                 slope, slope_index)
                            do ia=1, ntiles_aspect
                               call get_aspect_value(n,c,r,ia,aspect_selected,&
                                    aspect, aspect_index)
                               call get_curvature_value(n,c,r,ia,curvature_selected,&
                                    curvature, curvature_index)
                            
                               kk = kk+1

                               LIS_domain(n)%tile(kk)%ensem = m
                               LIS_domain(n)%tile(kk)%row=r    
                               LIS_domain(n)%tile(kk)%col=c    
                               LIS_domain(n)%tile(kk)%index = &
                                    LIS_domain(n)%gindex(c,r)
                               LIS_domain(n)%tile(kk)%vegt=vegt
                               LIS_domain(n)%tile(kk)%soilt=soilt
                               LIS_domain(n)%tile(kk)%sand=sand
                               LIS_domain(n)%tile(kk)%clay=clay
                               LIS_domain(n)%tile(kk)%silt=silt
                               LIS_domain(n)%tile(kk)%sftype = &
                                    sf_index
                               LIS_domain(n)%tile(kk)%elev = &
                                    elev
                               LIS_domain(n)%tile(kk)%slope = &
                                    slope
                               LIS_domain(n)%tile(kk)%aspect = &
                                    aspect
                               LIS_domain(n)%tile(kk)%curvature = &
                                    curvature
                               LIS_domain(n)%tile(kk)%fgrd = &
                                    LIS_LMLC(n)%landcover(c,r,vegt)
                               
                               if(soilt_selected) then 
                                  LIS_domain(n)%tile(kk)%fgrd = & 
                                       LIS_domain(n)%tile(kk)%fgrd * &
                                       LIS_soils(n)%texture(c,r,soilt)
                               elseif(soilf_selected) then 
                                  LIS_domain(n)%tile(kk)%fgrd = & 
                                       LIS_domain(n)%tile(kk)%fgrd * &
                                       LIS_soils(n)%soilffgrd(c,r,soilf_index)
                               endif
                               
                               if(elev_selected) then 
                                  LIS_domain(n)%tile(kk)%fgrd=& 
                                       LIS_domain(n)%tile(kk)%fgrd*&
                                       LIS_topo(n)%elevfgrd(c,r,elev_index)
                               endif

                               if(slope_selected) then 
                                  LIS_domain(n)%tile(kk)%fgrd=& 
                                       LIS_domain(n)%tile(kk)%fgrd*&
                                       LIS_topo(n)%slopefgrd(c,r,slope_index)
                               endif
                            
                               if(aspect_selected) then 
                                  LIS_domain(n)%tile(kk)%fgrd=& 
                                       LIS_domain(n)%tile(kk)%fgrd*&
                                       LIS_topo(n)%aspectfgrd(c,r,aspect_index)
                               endif

!                               if(curvature_selected) then
!                                  LIS_domain(n)%tile(kk)%fgrd=&
!                                       LIS_domain(n)%tile(kk)%fgrd*&
!                                       LIS_topo(n)%curvfgrd(c,r,curvature_index)
!                               endif
                               
                               if(sf_index.gt.0) then 
                                  kk_sf(sf_index) = kk_sf(sf_index) + 1
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%ensem = m
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%tile_id = kk
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%row=r    
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%col=c    
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%index = &
                                       LIS_domain(n)%gindex(c,r)
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%vegt=vegt
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%soilt=soilt
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%sand=sand
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%clay=clay
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%silt=silt

                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%elev=elev
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%slope=slope
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%aspect=aspect
                                  LIS_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%fgrd=&
                                       LIS_LMLC(n)%landcover(c,r,vegt)

                                  if(soilt_selected) then 
                                     LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LIS_soils(n)%texture(c,r,soilt)
                                  elseif(soilf_selected) then 
                                     LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LIS_soils(n)%soilffgrd(c,r,soilf_index)
                                  endif

                                  if(elev_selected) then 
                                     LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LIS_topo(n)%elevfgrd(c,r,elev_index)
                                  endif
                                  if(slope_selected) then 
                                     LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LIS_topo(n)%slopefgrd(c,r,slope_index)
                                  endif
                                  if(aspect_selected) then 
                                     LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LIS_topo(n)%aspectfgrd(c,r,aspect_index)
                                  endif                                  
                                  if(curvature_selected) then
                                     LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = &
                                          LIS_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * &
                                          LIS_topo(n)%curvfgrd(c,r,curvature_index)
                                  endif

                               endif

                            enddo
                         end do
                      end do
                   end do
                enddo
             end if
          end do
       end do
    enddo

#if 0 
!normalization check 
    allocate(sumv(LIS_rc%ngrid(n)))
    sumv = 0 
    do t=1,LIS_rc%ntiles(n)
       c = LIS_domain(n)%tile(t)%col
       r = LIS_domain(n)%tile(t)%row
       gid = LIS_domain(n)%gindex(c,r)
       
       sumv(gid) = sumv(gid) + LIS_domain(n)%tile(t)%fgrd
    enddo
    
    do t=1,LIS_rc%ngrid(n)
       print*, t, sumv(gid)
    enddo
    stop
#endif
    call LIS_domain_setup(n)

  end subroutine create_tilespace

!BOP
! 
! !ROUTINE: isSurfaceTypeSelected
! \label{isSurfaceTypeSelected}
! 
! !INTERFACE: 
  function isSurfaceTypeSelected(surface_type)
! !ARGUMENTS:     
    integer       :: surface_type
! 
! !DESCRIPTION:     
! This function determines if the input surface type is among 
! the surface model types that are currently chosen in the 
! LIS run
!
!EOP
    logical       :: isSurfaceTypeSelected
    integer       :: i 

    isSurfaceTypeSelected = .false. 
    
    do i=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(i).eq.surface_type) then 
          isSurfaceTypeSelected = .true. 
          exit
       endif
    enddo

  end function isSurfaceTypeSelected

!BOP
! 
! !ROUTINE: ifContainsSelectedSurfaceType
! \label{ifContainsSelectedSurfaceType}
! 
! !INTERFACE: 
  function ifContainsSelectedSurfaceType(surface_types)
! !ARGUMENTS:     
    integer       :: surface_types(LIS_rc%nsurfacetypes)
! 
! !DESCRIPTION: 
!  This function determines if any of the input list of surface types
!  is among the surface model types that are currently chosen in the 
!  LIS run
!EOP
    logical       :: ifContainsSelectedSurfaceType
    integer       :: i,k

    ifContainsSelectedSurfaceType = .false. 
    
    do k=1,LIS_rc%nsurfacetypes       
       do i=1,LIS_rc%nsf_model_types
          if(LIS_rc%sf_model_type_select(i).eq.surface_types(k)) then 
             ifContainsSelectedSurfaceType = .true. 
             exit
          endif
       enddo
    enddo

  end function ifContainsSelectedSurfaceType


!BOP
! 
! !ROUTINE: isValidSurfacePoint
! \label{isValidSurfacePoint}
! 
! !INTERFACE: 
  function isValidSurfacePoint(landcover,surface_types)
! !ARGUMENTS:     
    real       :: surface_types(LIS_rc%nsurfacetypes)
    real       :: landcover(LIS_rc%nsurfacetypes)
! 
! !DESCRIPTION: 
!  This function determines if any of the input list of surface types
!  is among the surface model types that are currently chosen in the 
!  LIS run
!EOP
    logical       :: isValidSurfacePoint,selected
    integer       :: i,k

    isValidSurfacePoint = .false. 
    
    do k=1,LIS_rc%nsurfacetypes       
       selected = .false.
       do i=1,LIS_rc%nsf_model_types
          if(LIS_rc%sf_model_type_select(i).eq.surface_types(k)) then 
             selected = .true. 
             exit
          endif
       enddo
       if(selected.and.landcover(k).gt.0) then 
          isValidSurfacePoint = .true. 
       endif
    enddo

  end function isValidSurfacePoint

!BOP
! 
! !ROUTINE: compute_ntiles_surface
! \label{compute_ntiles_surface}
! 
! !INTERFACE: 
  function compute_ntiles_surface(n,c,r)
! !ARGUMENTS: 
    integer :: n 
    integer :: c
    integer :: r
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of landcover types in a given grid cell
!EOP

    integer :: compute_ntiles_surface
    integer :: t
    
    compute_ntiles_surface = 0 
    do t=1,LIS_rc%nsurfacetypes        
       if(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
            LIS_LMLC(n)%landcover(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(&
            nint(LIS_LMLC(n)%surfacetype(c,r,t)))) then
          compute_ntiles_surface = compute_ntiles_surface + 1
       endif
    enddo
  end function compute_ntiles_surface

!BOP
! !ROUTINE: compute_npatches_surface
! \label{compute_npatches_surface}
! 
! !INTERFACE: 
  function compute_npatches_surface(n,c,r,t,sf_index)
! !ARGUMENTS: 
    integer :: n 
    integer :: c
    integer :: r
    integer :: t
    integer :: sf_index
!
! !DESCRIPTION: 
!  This function computes the number of patches (for the given 
!  surface type) based on the distribution of landcover types 
!  in a given grid cell
!EOP
    integer :: compute_npatches_surface
    
    compute_npatches_surface = 0 
    if(LIS_LMLC(n)%dommask(c,r).gt.0.and.&
         LIS_LMLC(n)%landcover(c,r,t).gt.0.0.and.&
         isSurfaceTypeSelected(&
         nint(LIS_LMLC(n)%surfacetype(c,r,t)))) then
       compute_npatches_surface = compute_npatches_surface + 1
    endif
  end function compute_npatches_surface

!BOP
! !ROUTINE: compute_ntiles_soilt
! \label{compute_ntiles_soilt}
! 
! !INTERFACE: 
  function compute_ntiles_soilt(n,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of soil texture types in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_soilt
    integer :: t

    kk = 0 
    if(flag) then 
       do t=1,LIS_rc%nsoiltypes        
          if(LIS_soils(n)%texture(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_soilt = max(1,kk)

  end function compute_ntiles_soilt

!BOP
! !ROUTINE: compute_ntiles_soilf
! \label{compute_ntiles_soilf}
! 
! !INTERFACE: 
  function compute_ntiles_soilf(n,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of soil fraction data in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_soilf
    integer :: t

    kk = 0 
    if(flag) then 
       do t=1,LIS_rc%nsoilfbands
          if(LIS_soils(n)%soilffgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_soilf = max(1,kk)

  end function compute_ntiles_soilf

!BOP
!
! !ROUTINE: compute_ntiles_elev
! \label{compute_ntiles_elev}
! 
! !INTERFACE: 
  function compute_ntiles_elev(n,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of elevation bands in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_elev
    integer :: t

    kk = 0 
    
    if(flag) then 
       do t=1,LIS_rc%nelevbands        
          if(LIS_topo(n)%elevfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_elev = max(1,kk)

  end function compute_ntiles_elev

!BOP
!
! !ROUTINE: compute_ntiles_slope
! \label{compute_ntiles_slope}
! 
! !INTERFACE: 
  function compute_ntiles_slope(n,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: c
    integer :: r
    logical :: flag
!
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of slope bands in a given grid cell
!EOP 
    integer :: kk 
    integer :: compute_ntiles_slope
    integer :: t

    kk = 0 
    
    if(flag) then 
       do t=1,LIS_rc%nslopebands        
          if(LIS_topo(n)%slopefgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_slope = max(1,kk)

  end function compute_ntiles_slope

!BOP
!
! !ROUTINE: compute_ntiles_aspect
! \label{compute_ntiles_aspect}
! 
! !INTERFACE: 
  function compute_ntiles_aspect(n,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of aspect bands in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_aspect
    integer :: t

    kk = 0 
    
    if(flag) then 
       do t=1,LIS_rc%naspectbands        
          if(LIS_topo(n)%aspectfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_aspect = max(1,kk)

  end function compute_ntiles_aspect

!BOP
!
! !ROUTINE: compute_ntiles_curvature
! \label{compute_ntiles_curvature}
! 
! !INTERFACE: 
  function compute_ntiles_curvature(n,c,r,flag)
! !ARGUMENTS: 
    integer :: n
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of curvature tiles in a given grid cell
!EOP
    integer :: kk
    integer :: compute_ntiles_curvature
    integer :: t

    kk = 0

    if(flag) then
!       do t=1,LIS_rc%ncurvbands
       do t=1,1
          if(LIS_topo(n)%curvfgrd(c,r,t).gt.0.0) then
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_curvature = max(1,kk)

  end function compute_ntiles_curvature


!BOP
! 
! !ROUTINE: get_vegt_value
! \label{get_vegt_value}
! 
! !INTERFACE: 
  subroutine get_vegt_value(n,c,r,i,vegt)
! !ARGUMENTS: 
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    integer  :: vegt
! 
! !DESCRIPTION: 
!  This routine computes the vegetation type value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk

    kk = 0 
    do t=1, LIS_rc%nsurfacetypes
       if(LIS_LMLC(n)%landcover(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(&
            nint(LIS_LMLC(n)%surfacetype(c,r,t)))) then 
          kk = kk + 1
          if(kk.eq.i) then 
             vegt = t
             return
          endif
       endif
    enddo
  end subroutine get_vegt_value

!BOP
! !ROUTINE: get_surface_value
! \label{get_surface_value}
! 
! !INTERFACE:
  subroutine get_surface_value(n,c,r,i,sf_index)
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    integer  :: sf_index
! 
! !DESCRIPTION: 
!  This routine computes the surface type value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk

    kk = 0 
    do t=1, LIS_rc%nsurfacetypes
       if(LIS_LMLC(n)%landcover(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(&
            nint(LIS_LMLC(n)%surfacetype(c,r,t)))) then 
          kk = kk + 1
          if(kk.eq.i) then 
             sf_index = nint(LIS_LMLC(n)%surfacetype(c,r,t))
             return
          endif
       endif
    end do
  end subroutine get_surface_value

!BOP
! !ROUTINE: get_soilt_value
! \label{get_soilt_value}
!
! !INTERFACE:
  subroutine get_soilt_value(n,c,r,i,soilt_selected,soilt)
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: soilt_selected
    integer  :: soilt
! 
! !DESCRIPTION: 
!  This routine computes the soil texture value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    soilt = -1
    if(soilt_selected) then 
       kk = 0 
       do t=1, LIS_rc%nsoiltypes
          if(LIS_soils(n)%texture(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                soilt = t
                return
             endif
          end if
       enddo
       if(soilt.eq.-1) then 
          print*, c,r,LIS_soils(n)%texture(c,r,:)
       endif
    endif

  end subroutine get_soilt_value

!BOP
! !ROUTINE: get_soilf_value
! \label{get_soilf_value}
!
! !INTERFACE:
  subroutine get_soilf_value(n,c,r,i,soilf_selected,sand,clay,silt,soilf_index)
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: soilf_selected
    real     :: sand
    real     :: clay
    real     :: silt
    integer  :: soilf_index
! 
! !DESCRIPTION: 
!  This routine computes the elevation value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    sand = -1
    clay = -1
    silt = -1
    soilf_index = -1
    if(soilf_selected) then 
       kk = 0 
       do t=1, LIS_rc%nsoilfbands
          if(LIS_soils(n)%soilffgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                sand = LIS_soils(n)%sand(c,r,t)
                clay = LIS_soils(n)%clay(c,r,t)
                silt = LIS_soils(n)%silt(c,r,t)
                soilf_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_soilf_value

!BOP
! !ROUTINE: get_elev_value
! \label{get_elev_value}
!
! !INTERFACE:
  subroutine get_elev_value(n,c,r,i,elev_selected,elev,elev_index)
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: elev_selected
    real     :: elev
    integer  :: elev_index
! 
! !DESCRIPTION: 
!  This routine computes the elevation value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    elev = -1
    elev_index = -1
    if(elev_selected) then 
       kk = 0 
       do t=1, LIS_rc%nelevbands
          if(LIS_topo(n)%elevfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                elev = LIS_topo(n)%elevation(c,r,t)
                elev_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_elev_value

!BOP
! 
! !ROUTINE: get_slope_value
! \label{get_slope_value}
!
! !INTERFACE:
  subroutine get_slope_value(n,c,r,i,slope_selected,slope,slope_index)
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: slope_selected
    real     :: slope
    integer  :: slope_index
! 
! !DESCRIPTION: 
!  This routine computes the slope value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    slope = -1
    slope_index = -1
    if(slope_selected) then 
       kk = 0 
       do t=1, LIS_rc%nslopebands
          if(LIS_topo(n)%slopefgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                slope = LIS_topo(n)%slope(c,r,t)
                slope_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_slope_value

!BOP
! !ROUTINE: get_aspect_value
! \label{get_aspect_value}
! 
! !INTERFACE: 
  subroutine get_aspect_value(n,c,r,i,aspect_selected,aspect,aspect_index)
! !ARGUMENTS: 
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: aspect_selected
    real     :: aspect
    integer  :: aspect_index
! 
! !DESCRIPTION: 
!  This routine computes the aspect value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    aspect = -1
    aspect_index = -1
    if(aspect_selected) then 
       kk = 0 
       do t=1, LIS_rc%naspectbands
          if(LIS_topo(n)%aspectfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                aspect = LIS_topo(n)%aspect(c,r,t)
                aspect_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_aspect_value

!BOP
! !ROUTINE: get_curvature_value
! \label{get_curvature_value}
! 
! !INTERFACE: 
  subroutine get_curvature_value(n,c,r,i,curvature_selected,curvature,curvature_index)
! !ARGUMENTS: 
    integer  :: n
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: curvature_selected
    real     :: curvature 
    integer  :: curvature_index
! 
! !DESCRIPTION: 
!  This routine computes the curvature value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk

    curvature = -1
    curvature = -1
    if(curvature_selected) then
       kk = 0
!       do t=1, LIS_rc%naspectbands
       do t=1, 1
          if(LIS_topo(n)%curvfgrd(c,r,t).gt.0.0) then
             kk = kk + 1
             if(kk.eq.i) then
                curvature = LIS_topo(n)%curvature(c,r,t)
                curvature_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_curvature_value

!BOP
! !ROUTINE: readDomainInput
! \label{readDomainInput}
! 
! !INTERFACE: 
  subroutine readDomainInput()
! 
! !DESCRIPTION: 
!  This routine reads the running domain information from the
!  LIS domain and parameter file.
!EOP

    integer                    :: n 
    integer                    :: ios
    integer                    :: ncId, nrId
    integer                    :: ncbId, nrbId
    integer                    :: latId,lonId
    integer                    :: latbId,lonbId
    integer                    :: gnc,gnr
    integer                    :: gnc_b, gnr_b
    integer                    :: ftn
    integer                    :: k, gindex
    logical                    :: file_exists
    character*50               :: map_proj
    integer                    :: c,r
    real,      allocatable     :: lat(:,:)
    real,      allocatable     :: lon(:,:)
    real,      allocatable     :: lat_b(:,:)
    real,      allocatable     :: lon_b(:,:)
    real,      allocatable     :: lat_b1(:,:)
    real,      allocatable     :: lon_b1(:,:)
    real,      allocatable     :: locallat(:,:)
    real,      allocatable     :: locallon(:,:)
    real,      allocatable     :: locallat_b(:,:)
    real,      allocatable     :: locallon_b(:,:)
    integer                    :: b_nc, b_nr
    real                       :: maxLat, minlat, maxLon, minLon
! _______________________________________________________

    write(LIS_logunit,fmt=*)'[INFO] DOMAIN details:'
    
    if ( LIS_rc%decompose_by_processes ) then
        write(LIS_logunit,*) '[INFO] Decomposing by processes. ' // &
                             '   Ignoring both'
        write(LIS_logunit,*) '[INFO]    "Number of processors along x:"'
        write(LIS_logunit,*) '[INFO]  and'
        write(LIS_logunit,*) '[INFO]    "Number of processors along y:"'
    else
       if(LIS_rc%npesx*LIS_rc%npesy.ne.LIS_npes) then
          write(unit=LIS_logunit,fmt=*) '[ERR] Layout does not match ' // &
                                        'the number of processors...'
          write(unit=LIS_logunit,fmt=*) '[ERR] npex, npey, ', &
                                        LIS_rc%npesx, 'x',    &
                                        LIS_rc%npesy, '!=', LIS_npes
          write(unit=LIS_logunit,fmt=*) '[ERR] Stopping program ..'
          call LIS_endrun()
       endif
    endif
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    do n=1,LIS_rc%nnest
       inquire(file=LIS_rc%paramfile(n), exist=file_exists)
       if(file_exists) then 
          
          ios = nf90_open(path=LIS_rc%paramfile(n),&
               mode=NF90_NOWRITE,ncid=ftn)
          call LIS_verify(ios,'Error in nf90_open in readDomainInput:paramfile')
          
          ios = nf90_inq_dimid(ftn,"east_west",ncId)
          call LIS_verify(ios,'Error in nf90_inq_dimid in readDomainInput:east_west')
        
          ios = nf90_inq_dimid(ftn,"north_south",nrId)
          call LIS_verify(ios,'Error in nf90_inq_dimid in readDomainInput:north_south')
          
          ios = nf90_inquire_dimension(ftn,ncId, len=gnc)
          call LIS_verify(ios,'Error in nf90_inquire_dimension in readDomainInput:ncId')
          
          ios = nf90_inquire_dimension(ftn,nrId, len=gnr)
          call LIS_verify(ios,'Error in nf90_inquire_dimension in readDomainInput:nrId')
          
          allocate(lat(gnc,gnr))
          allocate(lon(gnc,gnr))
        
          ios = nf90_inq_varid(ftn,'lat',latid)
          call LIS_verify(ios,'lat field not found in the LIS param file')
          
          ios = nf90_inq_varid(ftn,'lon',lonid)
          call LIS_verify(ios,'lon field not found in the LIS param file')
          
          ios = nf90_get_var(ftn,latid,lat)
          call LIS_verify(ios,'Error in nf90_get_var for latid in readDomainInput')
          
          ios = nf90_get_var(ftn,lonid,lon)
          call LIS_verify(ios,'Error in nf90_get_var for lonid in readDomainInput')
                  
          call LIS_quilt_domain(n, gnc,gnr)

          LIS_rc%gnc(n) = gnc
          LIS_rc%gnr(n) = gnr

          maxLat = -100000
          maxLon = -100000
          minLat = 100000
          minLon = 100000

          do r=1,gnr
             do c=1,gnc
                if(lat(c,r).gt.maxLat) then 
                   maxlat = lat(c,r)
                endif
                if(lon(c,r).gt.maxLon) then 
                   maxLon = lon(c,r)
                endif
                if(lat(c,r).lt.minLat) then 
                   minlat = lat(c,r)
                endif
                if(lon(c,r).lt.minLon) then 
                   minLon = lon(c,r)
                endif
             enddo
          enddo
          !default setting
          LIS_rc%minLat(n) = minLat
          LIS_rc%minLon(n) = minLon
          LIS_rc%maxLat(n) = maxLat
          LIS_rc%maxLon(n) = maxLon

          allocate(locallat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(locallon(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          
          locallat = lat(&
               LIS_ews_halo_ind(n,LIS_localPet+1):&         
               LIS_ewe_halo_ind(n,LIS_localPet+1), &
               LIS_nss_halo_ind(n,LIS_localPet+1): &
               LIS_nse_halo_ind(n,LIS_localPet+1))
          
          locallon = lon(&
               LIS_ews_halo_ind(n,LIS_localPet+1):&         
               LIS_ewe_halo_ind(n,LIS_localPet+1), &
               LIS_nss_halo_ind(n,LIS_localPet+1): &
               LIS_nse_halo_ind(n,LIS_localPet+1))
          
          allocate(LIS_domain(n)%lat(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(LIS_domain(n)%lon(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          do r=1,LIS_rc%lnr(n)
             do c=1,LIS_rc%lnc(n)
                gindex = c+(r-1)*LIS_rc%lnc(n)
                LIS_domain(n)%lat(gindex) = locallat(c,r)
                LIS_domain(n)%lon(gindex) = locallon(c,r)
             enddo
          enddo

          ! EMK NEW
          allocate(LIS_domain(n)%glat(LIS_rc%gnc(n)*LIS_rc%gnr(n))) ! EMK
          allocate(LIS_domain(n)%glon(LIS_rc%gnc(n)*LIS_rc%gnr(n))) ! EMK
          do r = 1, LIS_rc%gnr(n)
             do c = 1, LIS_rc%gnc(n)
                gindex = c+(r-1)*LIS_rc%gnc(n)
                LIS_domain(n)%glat(gindex) = lat(c,r)
                LIS_domain(n)%glon(gindex) = lon(c,r)
             end do
          end do
! b grid
          ios = nf90_inq_dimid(ftn,"east_west_b",ncbId)
          call LIS_verify(ios,'Error in nf90_inq_dimid in readDomainInput:east_west_b')
        
          ios = nf90_inq_dimid(ftn,"north_south_b",nrbId)
          call LIS_verify(ios,'Error in nf90_inq_dimid in readDomainInput:north_south_b')
          
          ios = nf90_inquire_dimension(ftn,ncbId, len=gnc_b)
          call LIS_verify(ios,'Error in nf90_inquire_dimension in readDomainInput:ncbId')
          
          ios = nf90_inquire_dimension(ftn,nrbId, len=gnr_b)
          call LIS_verify(ios,'Error in nf90_inquire_dimension in readDomainInput:nrbId')
          
          b_nc = (gnc_b - gnc)/2
          b_nr = (gnr_b - gnr)/2

          allocate(lat_b1(gnc_b,gnr_b))
          allocate(lon_b1(gnc_b,gnr_b))
          
          allocate(lat_b(-b_nc:gnc+b_nc,-b_nr:gnr+b_nr))
          allocate(lon_b(-b_nc:gnc+b_nc,-b_nr:gnr+b_nr))
        
          ios = nf90_inq_varid(ftn,'lat_b',latbid)
          call LIS_verify(ios,'lat_b field not found in the LIS param file')
          
          ios = nf90_inq_varid(ftn,'lon_b',lonbid)
          call LIS_verify(ios,'lon_b field not found in the LIS param file')
          
          ios = nf90_get_var(ftn,latbid,lat_b1)
          call LIS_verify(ios,&
               'Error in nf90_get_var for latb in readDomainInput:latbid')
          
          ios = nf90_get_var(ftn,lonbid,lon_b1)
          call LIS_verify(ios,&
               'Error in nf90_get_var for lonb in readDomainInput:lonbid')
          
          do r=-b_nr+1, gnr+b_nr
             do c=-b_nc+1, gnc+b_nc
                lat_b(c,r) = lat_b1(c+b_nc,r+b_nr)
                lon_b(c,r) = lon_b1(c+b_nc,r+b_nr)
             enddo
          enddo
                
          deallocate(lat_b1)
          deallocate(lon_b1)

          call LIS_quilt_b_domain(n, gnc_b,gnr_b, gnc, gnr)

          LIS_rc%gnc_b(n) = gnc_b
          LIS_rc%gnr_b(n) = gnr_b

          allocate(locallat_b(LIS_rc%lnc_b(n),LIS_rc%lnr_b(n)))
          allocate(locallon_b(LIS_rc%lnc_b(n),LIS_rc%lnr_b(n)))
          
          locallat_b = lat_b(&
               LIS_ews_b_ind(n,LIS_localPet+1):&         
               LIS_ewe_b_ind(n,LIS_localPet+1), &
               LIS_nss_b_ind(n,LIS_localPet+1): &
               LIS_nse_b_ind(n,LIS_localPet+1))
          
          locallon_b = lon_b(&
               LIS_ews_b_ind(n,LIS_localPet+1):&         
               LIS_ewe_b_ind(n,LIS_localPet+1), &
               LIS_nss_b_ind(n,LIS_localPet+1): &
               LIS_nse_b_ind(n,LIS_localPet+1))
          
          deallocate(lat_b)
          deallocate(lon_b)

          allocate(LIS_domain(n)%lat_b(LIS_rc%lnc_b(n)*LIS_rc%lnr_b(n)))
          allocate(LIS_domain(n)%lon_b(LIS_rc%lnc_b(n)*LIS_rc%lnr_b(n)))
          
          do r=1,LIS_rc%lnr_b(n)
             do c=1,LIS_rc%lnc_b(n)
                gindex = c+(r-1)*LIS_rc%lnc_b(n)
                LIS_domain(n)%lat_b(gindex) = locallat_b(c,r)
                LIS_domain(n)%lon_b(gindex) = locallon_b(c,r)
             enddo
          enddo

          deallocate(locallat_b)
          deallocate(locallon_b)

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',map_proj)
          call LIS_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')
          
          ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
               LIS_domain(n)%stlat)
          call LIS_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
          
          ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
               LIS_domain(n)%stlon)
          call LIS_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')
          
          ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',LIS_domain(n)%dx)
          call LIS_verify(ios, 'Error in nf90_get_att: DX')

          ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',LIS_domain(n)%dy)
          call LIS_warning(ios, 'Error in nf90_get_att: DY')

          if(ios.ne.0) then 
             LIS_domain(n)%dy = 0.0
          endif

          LIS_rc%gridDesc(n,32) = gnc
          LIS_rc%gridDesc(n,33) = gnr
          LIS_rc%gridDesc(n,34) = lat(1,1)
          LIS_rc%gridDesc(n,35) = lon(1,1)
          LIS_rc%gridDesc(n,37) = lat(gnc,gnr)
          LIS_rc%gridDesc(n,38) = lon(gnc,gnr)
          LIS_rc%gridDesc(n,39) = LIS_domain(n)%dx
          LIS_rc%gridDesc(n,40) = LIS_domain(n)%dy


          if(map_proj.eq."EQUIDISTANT CYLINDRICAL") then              
             LIS_rc%lis_map_proj = "latlon"

             LIS_rc%gridDesc(n,4) = LIS_domain(n)%stlat + &
                  (LIS_nss_halo_ind(n,LIS_localPet+1)-1)*LIS_domain(n)%dy
             LIS_rc%gridDesc(n,5) = LIS_domain(n)%stlon + &
                  (LIS_ews_halo_ind(n,LIS_localPet+1)-1)*LIS_domain(n)%dx
             LIS_rc%gridDesc(n,7) = LIS_domain(n)%stlat + & 
                  (LIS_nse_halo_ind(n,LIS_localPet+1)-1)*LIS_domain(n)%dy
             LIS_rc%gridDesc(n,8) = LIS_domain(n)%stlon + &
                  (LIS_ewe_halo_ind(n,LIS_localPet+1)-1)*LIS_domain(n)%dx
             
             LIS_rc%minLat(n) = LIS_rc%gridDesc(n,4) 
             LIS_rc%minLon(n) = LIS_rc%gridDesc(n,5) 
             LIS_rc%maxLat(n) = LIS_rc%gridDesc(n,7) 
             LIS_rc%maxLon(n) = LIS_rc%gridDesc(n,8) 

             write(unit=LIS_logunit,fmt=*) '[INFO] local domain ',&
                  LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,7),&
                  LIS_rc%gridDesc(n,5),LIS_rc%gridDesc(n,8), 'd: ',n
             
             LIS_rc%gridDesc(n,1) = 0
             LIS_rc%gridDesc(n,9) = LIS_domain(n)%dx
       
             if(LIS_rc%gridDesc(n,1).eq.0) then 
                LIS_rc%gridDesc(n,10) = LIS_domain(n)%dy
                LIS_rc%gridDesc(n,6) = 128
                LIS_rc%gridDesc(n,11) = 64
                LIS_rc%gridDesc(n,20) = 64
             endif
             if(LIS_rc%gridDesc(n,7).lt.LIS_rc%gridDesc(n,4)) then
                write(unit=LIS_logunit,fmt=*) '[ERR] lat2 must be greater than lat1'
                write(LIS_logunit,*) '[ERR] ',LIS_rc%gridDesc(n,7),LIS_rc%gridDesc(n,4)
                write(unit=LIS_logunit,fmt=*) '[ERR] Stopping run...'
                call LIS_endrun
             endif
             if(LIS_rc%gridDesc(n,8).lt.LIS_rc%gridDesc(n,5)) then
                write(unit=LIS_logunit,fmt=*) '[ERR] lon2 must be greater than lon1'
                write(LIS_logunit,*) '[ERR] ',LIS_rc%gridDesc(n,8),LIS_rc%gridDesc(n,5)
                write(unit=LIS_logunit,fmt=*) '[ERR] Stopping run...'
                call LIS_endrun
             endif
             
             LIS_rc%gridDesc(n,2) = nint((LIS_rc%gridDesc(n,8)-LIS_rc%gridDesc(n,5))&
                  /LIS_rc%gridDesc(n,9))+ 1
             LIS_rc%gridDesc(n,3) = nint((LIS_rc%gridDesc(n,7)-LIS_rc%gridDesc(n,4))&
                  /LIS_rc%gridDesc(n,10)) + 1  
             
                          
             !local domain of each processor
             call map_set(PROJ_LATLON, LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
                  0.0, LIS_rc%gridDesc(n,9),LIS_rc%gridDesc(n,10), 0.0,&
                  LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)
                          
             
          elseif(map_proj.eq."LAMBERT CONFORMAL") then      
             LIS_rc%lis_map_proj = "lambert"
             
             ! CM Ensure that the number of lat/lon dimensions is 2D for this projection
             if(LIS_rc%nlatlon_dimensions == '1D') then
               write(LIS_logunit,*) &
                  '[ERR] The lambert map projection cannot be written with 1D lat/lon fields.'
               write(LIS_logunit,*) &
                  '[WARN] The lat/lon fields will be written in 2D'
               LIS_rc%nlatlon_dimensions = '2D'
             end if

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',LIS_domain(n)%truelat1)
             call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT1')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT2',LIS_domain(n)%truelat2)
             call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT2')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',LIS_domain(n)%truelon)
             call LIS_verify(ios, 'Error in nf90_get_att: STANDARD_LON')

             LIS_rc%gridDesc(n,1) = 3
             LIS_rc%gridDesc(n,2) = LIS_rc%lnc(n)
             LIS_rc%gridDesc(n,3) = LIS_rc%lnr(n)
             LIS_rc%gridDesc(n,4) = lat(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,5) = lon(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,6) = 8
             LIS_rc%gridDesc(n,7) = LIS_domain(n)%truelat2
             LIS_rc%gridDesc(n,8) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,9) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,10) = LIS_domain(n)%truelat1
             LIS_rc%gridDesc(n,11) = LIS_domain(n)%truelon

             call map_set(PROJ_LC,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
                  LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),&
                  LIS_rc%gridDesc(n,10),LIS_domain(n)%truelat2,&
                  LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)


          elseif(map_proj.eq."MERCATOR") then              

             ! CM Ensure that the number of lat/lon dimensions is 2D for this projection
             if(LIS_rc%nlatlon_dimensions == '1D') then
               write(LIS_logunit,*) &
                  '[ERR] The MERCATOR map projection cannot be written with 1D lat/lon fields.'
               write(LIS_logunit,*) &
                  '[WARN] The lat/lon fields will be written in 2D'
               LIS_rc%nlatlon_dimensions = '2D'
             end if

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',LIS_domain(n)%truelat1)
             call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT1')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',LIS_domain(n)%truelon)
             call LIS_verify(ios, 'Error in nf90_get_att: STANDARD_LON')

             LIS_rc%gridDesc(n,1) = 1
             LIS_rc%gridDesc(n,2) = LIS_rc%lnc(n)
             LIS_rc%gridDesc(n,3) = LIS_rc%lnr(n)
             LIS_rc%gridDesc(n,4) = lat(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,5) = lon(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,6) = 8
             LIS_rc%gridDesc(n,7) = 0.0
             LIS_rc%gridDesc(n,8) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,9) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,10) = LIS_domain(n)%truelat1
             LIS_rc%gridDesc(n,11) = LIS_domain(n)%truelon

             call map_set(PROJ_MERC,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
                  LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),&
                  LIS_rc%gridDesc(n,10),&
                  0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)
             
          elseif(map_proj.eq."POLAR STEREOGRAPHIC") then              

             ! CM Ensure that the number of lat/lon dimensions is 2D for this projection 
             if(LIS_rc%nlatlon_dimensions == '1D') then
               write(LIS_logunit,*) &
                  '[ERR] The polar stereographic map projection cannot be written with 1D lat/lon fields.'
               write(LIS_logunit,*) &
                  '[WARN] The lat/lon fields will be written in 2D'
               LIS_rc%nlatlon_dimensions = '2D'
             end if

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',LIS_domain(n)%truelat1)
             call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT1')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'ORIENT',LIS_domain(n)%orient)
             call LIS_verify(ios, 'Error in nf90_get_att: ORIENT')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',LIS_domain(n)%truelon)
             call LIS_verify(ios, 'Error in nf90_get_att: STANDARD_LON')

             LIS_rc%gridDesc(n,1) = 5
             LIS_rc%gridDesc(n,2) = LIS_rc%lnc(n)
             LIS_rc%gridDesc(n,3) = LIS_rc%lnr(n)
             LIS_rc%gridDesc(n,4) = lat(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,5) = lon(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,6) = 8
             LIS_rc%gridDesc(n,7) = LIS_domain(n)%orient
             LIS_rc%gridDesc(n,8) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,9) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,10) = LIS_domain(n)%truelat1
             LIS_rc%gridDesc(n,11) = LIS_domain(n)%truelon

             call map_set(PROJ_PS,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
                  LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),&
                  0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)
             
          elseif(map_proj.eq."GAUSSIAN") then 

             ! CM Ensure that the number of lat/lon dimensions is 2D for this projection
             if(LIS_rc%nlatlon_dimensions == '1D') then
               write(LIS_logunit,*) &
                  '[ERR] The gaussian map projection cannot be written with 1D lat/lon fields.'
               write(LIS_logunit,*) &
                  '[WARN] The lat/lon fields will be written in 2D'
               LIS_rc%nlatlon_dimensions = '2D'
             end if

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'NUMBER OF LAT CIRCLES',&
                  LIS_domain(n)%nlatcircles)
             call LIS_verify(ios, 'Error in nf90_get_att: NUMBER_OF_LAT_CIRCLES')

             LIS_rc%gridDesc(n,1) = 4
             LIS_rc%gridDesc(n,2) = LIS_rc%lnc(n)
             LIS_rc%gridDesc(n,3) = LIS_rc%lnr(n)
             LIS_rc%gridDesc(n,4) = lat(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,5) = lon(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,6) = 128
             LIS_rc%gridDesc(n,7) = lat(LIS_ewe_halo_ind(n,LIS_localPet+1),&
                  LIS_nse_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,8) = lon(LIS_ewe_halo_ind(n,LIS_localPet+1),&
                  LIS_nse_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,9) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,10) = LIS_domain(n)%nlatcircles
             LIS_rc%gridDesc(n,11) = 64
             LIS_rc%gridDesc(n,20) = 64

          elseif(map_proj.eq."HRAP") then 
   
             ! CM Ensure that the number of lat/lon dimensions is 2D for this projection
             if(LIS_rc%nlatlon_dimensions == '1D') then
               write(LIS_logunit,*) &
                  '[ERR] The HRAP map projection cannot be written with 1D lat/lon fields.'
               write(LIS_logunit,*) &
                  '[WARN] The lat/lon fields will be written in 2D'
               LIS_rc%nlatlon_dimensions = '2D'
             end if

             LIS_rc%gridDesc(n,1) = 8
             LIS_rc%gridDesc(n,2) = LIS_rc%lnc(n)
             LIS_rc%gridDesc(n,3) = LIS_rc%lnr(n)
             LIS_rc%gridDesc(n,4) = lat(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,5) = lon(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,6) = 8
             LIS_rc%gridDesc(n,7) = LIS_domain(n)%orient
             LIS_rc%gridDesc(n,8) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,9) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,10) = LIS_domain(n)%truelat1
             LIS_rc%gridDesc(n,11) = LIS_domain(n)%truelon
             
             call map_set(PROJ_HRAP,LIS_rc%gridDesc(n,4),LIS_rc%gridDesc(n,5),&
                  LIS_rc%gridDesc(n,8)*1000.0,LIS_rc%gridDesc(n,11),LIS_rc%gridDesc(n,10),&
                  0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_domain(n)%lisproj)
             
          elseif(map_proj.eq."UTM") then 
        
             ! CM Ensure that the number of lat/lon dimensions is 2D for this projection
             if(LIS_rc%nlatlon_dimensions == '1D') then
               write(LIS_logunit,*) &
                  '[ERR] The UTM map projection cannot be written with 1D lat/lon fields.'
               write(LIS_logunit,*) &
                  '[WARN] The lat/lon fields will be written in 2D'
               LIS_rc%nlatlon_dimensions = '2D'
             end if 


             LIS_rc%gridDesc(n,1) = 7
             LIS_rc%gridDesc(n,2) = LIS_rc%lnc(n)
             LIS_rc%gridDesc(n,3) = LIS_rc%lnr(n)
             LIS_rc%gridDesc(n,4) = lat(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,5) = lon(LIS_ews_halo_ind(n,LIS_localPet+1),&
                  LIS_nss_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,6) = 128
             LIS_rc%gridDesc(n,7) = lat(LIS_ewe_halo_ind(n,LIS_localPet+1),&
                  LIS_nse_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,8) = lon(LIS_ewe_halo_ind(n,LIS_localPet+1),&
                  LIS_nse_halo_ind(n,LIS_localPet+1))
             LIS_rc%gridDesc(n,9) = LIS_domain(n)%dx
             LIS_rc%gridDesc(n,10) = LIS_domain(n)%dy
             LIS_rc%gridDesc(n,11) = 64
             LIS_rc%gridDesc(n,20) = 64
             
          endif

          deallocate(lat)
          deallocate(lon)
          deallocate(locallat)
          deallocate(locallon)

       else
          write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), ' does not exist'
          write(LIS_logunit,*) '[ERR] program stopping ...'
          call LIS_endrun
       endif

       write(LIS_logunit,*)'[INFO] --------------------Domain ',n,'----------------------'
       do k=1,13
          write(unit=LIS_logunit,fmt=*) '[INFO] (',k,',',LIS_rc%gridDesc(n,k),')'
       enddo
       write(LIS_logunit,*)'[INFO] --------------------------------------------------------------'
       
    enddo
#endif

  end subroutine readDomainInput


!BOP
! !ROUTINE: decompose_nx_ny
! \label{decompose_nx_ny}
!
! !INTERFACE:
subroutine decompose_nx_ny(nc, nr, ips, ipe, jps, jpe)
! !USES:
!  none

! !ARGUMENTS:
   integer, intent(in)  :: nc, nr
   integer, intent(out) :: ips, ipe, jps, jpe
!
! !DESCRIPTION:
! This routine decomposes a domain based on the specified
! processor layout (for example 4x3).
!
! The arguments are:
! \begin{description}
! \item[nc] number of columns of the domain
! \item[nr] number of columns of the domain
! \item[ips] starting column for the process
! \item[ipe] ending column for the process
! \item[jps] starting row for the process
! \item[jpe] ending row for the process
! \end{description}
!
!EOP

   integer :: mytask_x, mytask_y
   integer :: Px, Py, P
   integer :: i, j

#ifdef MPDECOMP2
    mytask_x = mod( LIS_localPet , LIS_rc%npesx ) + 1
    mytask_y = ( LIS_localPet / LIS_rc%npesx ) + 1

    j = 1
    ips = -1
    do i=1, nc
      call LIS_mpDecomp_2(i,j,1,nc,1,nr,LIS_rc%npesx, LIS_rc%npesy,Px,Py,P)
      if(Px.eq. mytask_x) then
        ipe = i
        if(ips.eq.-1) ips = i
      endif
    enddo

    i = 1
    jps = -1
    do j=1, nr
      call LIS_mpDecomp_2(i,j,1,nc,1,nr,LIS_rc%npesx,LIS_rc%npesy,Px,Py,P)
      if(Py.eq.mytask_y) then
        jpe = j
        if(jps.eq.-1) jps = j
      endif
    enddo
#else
   mytask_x = mod(LIS_localPet, LIS_rc%npesx)
   mytask_y = LIS_localPet / LIS_rc%npesx

   j = 1
   ips = -1
   do i = 1, nc+1
      call LIS_mpDecomp(i,j,1,nc+1,1,nr+1,LIS_rc%npesx,LIS_rc%npesy,Px,Py,P)

      if ( Px .eq. mytask_x ) then
         ipe = i
         if ( ips .eq. -1 ) then
            ips = i
         endif
      endif
   enddo

   i = 1
   jps = -1
   do j = 1, nr+1
      call LIS_mpDecomp(i,j,1,nc+1,1,nr+1,LIS_rc%npesx,LIS_rc%npesy,Px,Py,P)
      if ( Py .eq. mytask_y ) then
         jpe = j
         if ( jps .eq. -1 ) then
            jps = j
         endif
      endif
   enddo
#endif
end subroutine decompose_nx_ny


!BOP
! !ROUTINE: decompose_npes
! \label{decompose_npes}
!
! !INTERFACE:
subroutine decompose_npes(n, gnc, gnr, ips, ipe, jps, jpe)
! !USES:
!  none

! !ARGUMENTS:
   integer, intent(in)  :: n, gnc, gnr
   integer, intent(out) :: ips, ipe, jps, jpe

! !DESCRIPTION:
! This routine decomposes a domain based on the number of
! processing elements.
!
! The arguments are:
! \begin{description}
! \item[n] nest index
! \item[gnc] number of columns of the domain
! \item[gnr] number of columns of the domain
! \item[ips] starting column for the process
! \item[ipe] ending column for the process
! \item[jps] starting row for the process
! \item[jpe] ending row for the process
! \end{description}
!
! The local variables are:
! \item[file_exists]
!    logical indicating whether the LIS domain and parameter file
!    containing the dommask exists
! \item[ios]
!    return code from various routine calls
! \item[nid]
!    NetCDF handle to the LIS domain and parameter file
! \item[dommaskid]
!    NetCDF handle to the dommask field in the LIS domain and parameter file
! \item[sc, ec, sr, er]
!    starting column, ending column, starting row, and ending row
!    indicies used to sub-set into the landmask for counting the number
!    of land-based grid-cells
! \item[f1, f2]
!    the factors of LIS_npes used to create the layout of the running domain
! \item[Pc, Pr]
!    the number of processor columns and of processor rows, respectively
! \item[firstP, secondP]
!    The running domain is decomposed into a Pc-by-Pr layout in two passes,
!    either by Pc then by Pr or by Pr the by Pc, depending on the shape
!    of the running domain.  firstP is Pc (or Pr) and secondP is Pr (or Pc).
! \item[first\_size, second\_size]
!    first\_size is either gnc or gnr and second\_size is either gnr or gnc
!    depending on the order in which the processor layout is being processed.
! \item[pe]
!    a processor element
! \item[LP]
!    the total number of land-based grid-cells
! \item[T]
!    the per-process target for land-based grid-cells
! \item[gross\_target]
!    the target for land-based grid-cells for the first pass
! \item[i, j]
!    looping indicies
! \item[columnwise]
!    logical indicating whether a pass is along the Pc direction or the Pr.
!    columnwise == .true. is in the Pc (or column) direction.
! \item[lcount]
!    array containing either per-column or per-row counts of the land-based
!    grid-cells
! \item[s\_array, e\_array]
!    arrays containing the starting and ending indicies into the running
!    domain where processor cuts are to be made
! \item[ips\_array, ipe\_array]
!    arrays contains the starting and ending indicies into the running
!    domain in the column direction where processor cuts are to be made.
!    used to set ips and ipe.
! \item[jps\_array, jpe\_array]
!    arrays contains the starting and ending indicies into the running
!    domain in the row direction where processor cuts are to be made.
!    used to set jps and jpe.
! \item[istat]
!    return code from various MPI routine calls
!
!EOP


! == Overview of the load balancing decomposition scheme
!
! We know the number of requested processes, LIS_npes.  So find the factors
! of LIS_npes and choose the middle ones to create a squarish layout.
! Call these factors Pc and Pr, where Pc is the number of processors in
! the x or column direction and where Pr is the number of processors in
! the y or row direction.  We will create a Pc by Pr layout.
!
! Note that this decomposition scheme is *not* recommended for running
! domains where the ratio of gnc to gnr is very large (likewise where gnr to
! gnc is very large).  For running domains of this shape, use a specified
! layout; for example,
!
!    Decompose by processes:          .false.
!    Number of processors along x:    64
!    Number of processors along y:    1
!
! Note that this decomposition scheme balances with respect to the number
! of land-based grid-cells per process.  It does not take tiling or other
! considerations into account.
!
! Note that this decomposition scheme does *not* yet support nesting.
!
! Let's assume that we will decompose the running domain first in the
! column direction and then in the row direction.
!
!
! == First gross cut
!
! Determine number of land-based grid-cells per column, L_c_i, of the
! running domain.  Summing these values gives the total number of land-based
! grid-cells, LP.
!
! The target number of land-based grid-cells per process, T, is given by
!    T = LP / LIS_npes.
!
! Since we are decomposing in the column direction first, determine how
! many columns from the land mask satisfies
!    Sum L_c_i = Pr * T
!
! This gives the starting and ending indicies for cutting the running domain
! along the x or column direction.
!
!
! == Second finer cut
!
! There are Pc processor columns in the layout.  Refer to each processor
! column as Pc_i.
!
! For each processor column, Pc_i, determine the number of land-based
! grid-cells per row, L_r_j, given the number of actual columns in Pc_i.
! Determine how many rows from the land mask satisfies
!    Sum L_r_j = T
!
! This gives the starting and ending indicies for cutting the running domain
! along the y or row direction for each Pc_i.
!
! == Example
!
! For this example, we will assume that we will decompose in the x or column
! direction first.  Say that 12 processes are requested.  This will be
! factored into a Pc=4 by Pr=3 layout.
!
! The figure below represents both the gridded landmask, with dimensions
! gnc by gnr, and the load-balanced decomposition.  Dots ('.') represent
! the grid-cells (no distinction between land or water).  Dashes ('-'
! and '|') represent the boundaries or cuts between the processing
! elements, denoted P0, P1, ..., P11.  c_N and r_M denote columns and
! rows of the landmask.
!
! Since we are decomposing in the column direction first, we want to
! count the number of land-points in each column of the landmask for
! columns c_1 through c_{gnc}.  With these per-column counts we now know
! the total number of land-points in the landmask.  Call this total LP.
! Thus we want a per-process target, T, of T = LP / 12.
!
! Now let's consider the first gross cut.  Since we have a 4x3 layout,
! there are four processor columns, denoted Pc_1, Pc_2, Pc_3, and Pc_4.
! Pc_1 contains processes P0, P4, and P8.  Pc_2 contains P1, P5, and P9.
! Pc_3 contains P2, P6, and P10.  Pc_4 contains P3, P7, and P11.  Each
! Pc_i contains three processes, so when we start cutting in the column
! direction, we want each Pc_i to have 3 times the per-process target;
! i.e., we want Pr*T.  So we start counting how many columns into the
! landmask we must go to get Pr*T land-points.  Here that is column c_m.
! Then staring from c_{m+1}, we count how many columns into the landmask
! we must go to get another Pr*T land-points.  Here that is c_n.  Then
! staring from c_{n+1}, we count through column c_o.  And the last
! processor column is simply c_{o+1} through c_gnc.  Thus, Pc_1 spans
! landmask columns c_1 through c_m; Pc_2 spans columns c_{m+1} through
! c_n; Pc_3 spans columns c_{n+1} through c_o; Pc_4 spans columns
! c_{o+1} through c_{gnc}.
!
! The second finer cut loops over the processor columns, cutting each
! processor column into three processor rows.  So for each processor
! column, we count the number of land-points in each row of the landmask
! for rows r_1 through r_gnr constrained to the columns of that
! processor column.  For processor column 1, Pc_1, we count the number
! of land-points in each row of the landmask constrained to columns c_1
! through c_m.  With these values,  we count how many rows into the
! landmask we must go to get T land-points, where T is the per-process
! target.  Here that is row r_j.  Then starting from r_{j+1} we count
! how many rows into the landmask we must go to get another T
! land-points.  Here that is row r_k.  The final cut is simply from
! r_{k+1} through r_gnr.  This is repeated for Pc_2 through Pc_4.
!
!
!       <--   Pc_1  --><--    Pc_2    --><-- Pc_3 --><--    Pc_4    -->
!
! r_gnr +-------------++----------------++----------++----------------+
!       |.............||................||..........||................|
!       |.............||................||..........||................|
!       |.... P8 .....||...... P9 ......||.. P10 ...||................|
!       |.............|+----------------+|..........||................|
!       |.............|+----------------+|..........||.... P11 .......|
!       +-------------+|................||..........||................|
! r_k   +-------------+|................|+----------+|................|
!       |.............||................|+----------+|................|
!       |.............||................||..........|+----------------+
!       |.............||................||..........|+----------------+
!       |.... P4 .....||...... P5 ......||... P6 ...||................|
!       |.............||................||..........||................|
!       |.............|+----------------+|..........||................|
!       |.............|+----------------+|..........||..... P7 .......|
!       |.............||................||..........||................|
!       |.............||................|+----------+|................|
!       |.............||................|+----------+|................|
! r_j+1 +-------------+|................||..........||................|
! r_j   +-------------+|................||..........||................|
!       |.............||................||..........|+----------------+
!       |.............||................||..........|+----------------+
! .     |.............||................||..........||................|
! .     |..... P0 ....||...... P1 ......||... P2 ...||..... P3 .......|
! .     |.............||................||..........||................|
!       |.............||................||..........||................|
! r_2   |.............||................||..........||................|
! r_1   +-------------++----------------++----------++----------------+
!
!
!       cc  ...       cc                c           c                 c
!       12            mm                n           o               gnc
!                      +
!                      1
!
!-----------------------------------------------------------------------


   logical :: file_exists
   integer :: ios, nid, landmaskid
   integer :: sc, ec, sr, er
   integer :: f1, f2
   integer :: Pc, Pr, firstP, secondP, first_size, second_size, pe
   integer :: LP, T, gross_target, i, j, l1, l2
   logical :: columnwise
   !integer :: ncId, nrId, nc, nr
   integer, allocatable, dimension(:) :: lcount, s_array, e_array
   integer, allocatable, dimension(:) :: ips_array, ipe_array, &
                                         jps_array, jpe_array
#if ( defined SPMD )
   integer :: istat(MPI_STATUS_SIZE)
#endif

   ips = -9999
   ipe = -9999
   jps = -9999
   jpe = -9999

   if ( LIS_masterproc ) then

      inquire(file=trim(LIS_rc%paramfile(n)), exist=file_exists)

      if ( .not. file_exists ) then
         write(LIS_logunit,*) '[ERR] '//trim(LIS_rc%paramfile(n)), &
            ' does not exist'
         write(LIS_logunit,*) '[ERR] program stopping'
         call LIS_endrun
      endif

      ios = nf90_open(path=trim(LIS_rc%paramfile(n)), &
                      mode=NF90_NOWRITE,ncid=nid)
      call LIS_verify(ios,'Error in nf90_open in decompose_npes')

      !ios = nf90_inq_dimid(nid,"east_west",ncId)
      !call LIS_verify(ios,'Error in nf90_inq_dimid in decompose_npes')

      !ios = nf90_inq_dimid(nid,"north_south",nrId)
      !call LIS_verify(ios,'Error in nf90_inq_dimid in decompose_npes')

      !ios = nf90_inquire_dimension(nid,ncId,len=nc)
      !call LIS_verify(ios,'Error in nf90_inquire_dimension in decompose_npes')

      !ios = nf90_inquire_dimension(nid,nrId,len=nr)
      !call LIS_verify(ios,'Error in nf90_inquire_dimension in decompose_npes')

      ios = nf90_inq_varid(nid,'LANDMASK',landmaskid)
      call LIS_verify(ios,'landmask field not found in the LIS param file')

      allocate(ips_array(0:LIS_npes-1))
      allocate(ipe_array(0:LIS_npes-1))
      allocate(jps_array(0:LIS_npes-1))
      allocate(jpe_array(0:LIS_npes-1))

      ! Initialize for use by the first gross cut
      ips_array = 1
      ipe_array = gnc
      jps_array = 1
      jpe_array = gnr

      ! Determine layout
      call inner_factors(LIS_npes, f1, f2)
      if ( gnc >= gnr ) then
         Pc = f1
         Pr = f2
         firstP = Pc
         secondP = Pr
         first_size = gnc
         second_size = gnr
         columnwise = .true.
      else
         Pc = f2
         Pr = f1
         firstP = Pr
         secondP = Pc
         first_size = gnr
         second_size = gnc
         columnwise = .false.
      endif

      write(LIS_logunit,*) '[INFO] producing a ', Pc, ' by ', Pr, ' layout'

      ! Count global number of land-based grid-cells
      allocate(lcount(first_size))
      do i = 1, first_size
         if ( columnwise ) then
            call set_count_indicies(c=i, r=1)
         else
            call set_count_indicies(c=1, r=i)
         endif
         lcount(i) = count_mask(nid, landmaskid, sc, ec, sr, er)
      enddo
      LP = sum(lcount)
      T = nint(real(LP) / real(LIS_npes))
      if ( LP < 1 ) then
         write(LIS_logunit, *) '[ERR] cannot divide number of land ' // &
            'grid-cells over requested number of processes.'
         call LIS_endrun
      endif

      ! Perform first gross cut of the domain
      gross_target = secondP*T
      allocate(s_array(firstP))
      allocate(e_array(firstP))
      call cut(firstP, first_size, s_array, e_array, lcount, gross_target)
      do j = 1, Pr
         do i = 1, Pc
            pe = proc_element(c=i, r=j)
            if ( columnwise ) then
               ips_array(pe) = s_array(i)
               ipe_array(pe) = e_array(i)
            else
               jps_array(pe) = s_array(j)
               jpe_array(pe) = e_array(j)
            endif
         enddo
      enddo
      deallocate(s_array)
      deallocate(e_array)
      deallocate(lcount)

      ! Perform second finer cut of the domain
      allocate(lcount(second_size))
      allocate(s_array(secondP))
      allocate(e_array(secondP))

      ! Toggle value of columnwise
      columnwise = .not. columnwise

      do l1 = 1, firstP         ! layout loop
         do i = 1, second_size  ! domain loop
            if ( columnwise ) then
               call set_count_indicies(c=i, r=l1)
            else
               call set_count_indicies(c=l1, r=i)
            endif
            lcount(i) = count_mask(nid, landmaskid, sc, ec, sr, er)
         enddo
         call cut(secondP, second_size, s_array, e_array, lcount, T)
         do l2 = 1, secondP
            if ( columnwise ) then
               ! now count columnwise / l2 represents column index
               pe = proc_element(c=l2, r=l1)
               ips_array(pe) = s_array(l2)
               ipe_array(pe) = e_array(l2)
            else
               ! now count rowwise / l2 represents row index
               pe = proc_element(c=l1, r=l2)
               jps_array(pe) = s_array(l2)
               jpe_array(pe) = e_array(l2)
            endif
         enddo
      enddo
      deallocate(s_array)
      deallocate(e_array)
      deallocate(lcount)

#if ( defined SPMD )
      do i = 1, LIS_npes-1
         call MPI_SEND(ips_array(i), 1, MPI_INTEGER, i, 0, LIS_mpi_comm, ios)
         call MPI_SEND(ipe_array(i), 1, MPI_INTEGER, i, 0, LIS_mpi_comm, ios)
         call MPI_SEND(jps_array(i), 1, MPI_INTEGER, i, 0, LIS_mpi_comm, ios)
         call MPI_SEND(jpe_array(i), 1, MPI_INTEGER, i, 0, LIS_mpi_comm, ios)
      enddo
#endif
      ips = ips_array(0)
      ipe = ipe_array(0)
      jps = jps_array(0)
      jpe = jpe_array(0)
#if ( defined SPMD )
   else
      allocate(ips_array(1))
      allocate(ipe_array(1))
      allocate(jps_array(1))
      allocate(jpe_array(1))

      call MPI_RECV(ips_array, 1, MPI_INTEGER, 0, 0, LIS_mpi_comm, istat, ios)
      call MPI_RECV(ipe_array, 1, MPI_INTEGER, 0, 0, LIS_mpi_comm, istat, ios)
      call MPI_RECV(jps_array, 1, MPI_INTEGER, 0, 0, LIS_mpi_comm, istat, ios)
      call MPI_RECV(jpe_array, 1, MPI_INTEGER, 0, 0, LIS_mpi_comm, istat, ios)

      ips = ips_array(1)
      ipe = ipe_array(1)
      jps = jps_array(1)
      jpe = jpe_array(1)

      deallocate(ips_array)
      deallocate(ipe_array)
      deallocate(jps_array)
      deallocate(jpe_array)
#endif
   endif

   if ( LIS_masterproc ) then
      deallocate(ips_array)
      deallocate(ipe_array)
      deallocate(jps_array)
      deallocate(jpe_array)

      ios = nf90_close(nid)
      call LIS_verify(ios,'Error in nf90_close in decompose_npes')
   endif

   if ( ips == ipe .or. jps == jpe ) then
      write(LIS_logunit, *) '[ERR] decomposition for process ', LIS_localPet, &
                            ' is not at least 2x2'
      write(LIS_logunit, *) 'Stopping'
      call LIS_endrun
   endif

   contains

   integer function proc_element(c, r)
      implicit none
      integer, intent(in) :: c, r
      proc_element = Pc*(r-1)+c-1
   end function proc_element

   subroutine set_count_indicies(c, r)
      implicit none
      integer, intent(in) :: c, r
      integer :: pe

      ! When processing columnwise, c represents an index into the domain
      ! grid (1 through gnc); r represents an index into the processor
      ! layout (1 through Pr).
      ! When processing rowwise, c represents an index into the processor
      ! layout (1 through Pc); r represnets an index into the domain
      ! grid (1 through gnr).

      ! pe represents the processing element.  Usually you need both
      ! the column and row indices in the processor layout to determine
      ! the processing element, but due to the way the values in the
      ! ips_array, ipe_array, jps_array, and jpe_array arrays repeat,
      ! you only need a representative processing element not the
      ! actual one.  Hence the unusual equation for determining pe.

      if ( columnwise ) then
         pe = Pc*(r-1)
         sc = c; ec = c; sr = jps_array(pe); er = jpe_array(pe)
      else
         pe = c-1
         sc = ips_array(pe); ec = ipe_array(pe); sr = r; er = r
      endif
   end subroutine set_count_indicies
end subroutine decompose_npes


!BOP
! !ROUTINE: cut
! \label{cut}
!
! !INTERFACE:
subroutine cut(NP, asize, s_array, e_array, lcount, cut_target)
! !USES:
!  none
   implicit none

! !ARGUMENTS:
   integer, intent(in) :: NP, asize
   integer, intent(in), dimension(asize) :: lcount
   integer, intent(in) :: cut_target
   integer, intent(out), dimension(NP) :: s_array, e_array

! !DESCRIPTION:
! This routine determines the starting and ending indices for each process.
!
! This routine will be first called to perform a gross cut of the domain,
! where it determines the starting and ending indices for each processor
! column (or row depending on how called).  Then this routine will be called
! to determine the starting and ending row (column) indices for each process
! within a given processor column (row).
!
! The arguments are:
! \begin{description}
! \item[NP] number of processors, either Pc or Pr, where the domain
!    is laid out as Pc processes in x-direction by Pr processes
!    in the y-direction
! \item[asize] size of dimension, either gnc or gnr, where gnc is
!    the number of columns in the domain and gnr is the number of rows.
! \item[lcount] array of land-based grid-cell counts
! \item[cut\_target] target number of land-based grid-cells per cut.
! \item[s\_array] array of starting indices
! \item[e\_array] array of ending indices
! \end{description}
!
!EOP
   integer :: p, i
   integer :: slice_count, before, after

   s_array(1) = 1
   e_array(NP) = asize
   do p = 1, NP-1
      slice_count = 0
      do i = s_array(p), asize
         slice_count = slice_count + lcount(i)
         if ( slice_count >= cut_target ) then
            before = (slice_count - lcount(i))
            after = slice_count
            if ( abs(before-cut_target) < abs(after-cut_target) ) then
               e_array(p) = i - 1
               s_array(p+1) = i
            else
               e_array(p) = i
               s_array(p+1) = i + 1
            endif
            exit
         endif
      enddo
   enddo
end subroutine cut


!BOP
! !ROUTINE: inner_factors
! \label{inner_factors}
!
! !INTERFACE:
subroutine inner_factors(N, f1, f2)
! !USES:
!  none

   implicit none

! !ARGUMENTS:
   integer, intent(in)  :: N
   integer, intent(out) :: f1, f2

! !DESCRIPTION:
! This routine finds the inner-most factors of the input N.
! For example, the factors of N=255 are 1, 3, 5, 15, 17, 51, 85, 255.
! This routine finds f1=17 and f2=15.  This routine always returns f1 >= f2.
!
! The arguments are:
! \begin{description}
! \item[N] number to factor
! \item[f1] one factor
! \item[f2] corresponding other factor
! \end{description}
!
!EOP

   real :: rN
   real :: x
   integer :: ix
   integer :: i

   rN = real(N)
    x = sqrt(rN)
   ix = int(x)

   if ( ix**2 == N ) then
      ! N is square
      f1 = ix
      f2 = ix
   else
      do i = ix+1, N
         if ( mod(N, i) == 0 ) then
            ! found a factor
            f1 = i
            f2 = N / f1
            exit
         endif
      enddo
   endif
end subroutine inner_factors


!BOP
!
! !ROUTINE: count_mask
! \label{count_mask}
!
! !INTERFACE:
integer function count_mask(nid, landmaskid, sc, ec, sr, er) result(l_count)
! !USES:
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   use netcdf
#endif

   implicit none

! !ARGUMENTS:
   integer, intent(in)  :: nid, landmaskid, sc, ec, sr, er
!
! !DESCRIPTION:
!  Reads a column or row from the landmask and counts the number of land
!  grid-cells.
!
!  The arguments are:
!  \begin{description}
!  \item[nid]
!     handle of the LIS domain and parameter file (file containing the
!     landmask)
!  \item[landmaskid]
!     handle to the landmask field
!  \item[sc]
!     starting column
!  \item[ec]
!     ending column
!  \item[sr]
!     starting row
!  \item[er]
!     ending row
!  \item[l_count]
!     number of land grid-cells in the requested slice
!  \end{description}
!
!EOP

   integer :: ios, c_count, r_count, asize, i
   real, allocatable, dimension(:) :: array

   c_count = ec - sc + 1
   r_count = er - sr + 1

   if ( sc == ec ) then
      asize = r_count
   else
      asize = c_count
   endif
   allocate(array(asize))

   ios = nf90_get_var(nid, landmaskid, array,     &
                      start=(/sc,sr/),            &
                      count=(/c_count,r_count/))
   call LIS_verify(ios,'Error in nf90_get_var in count_mask')

   l_count = 0
   do i = 1, asize
      if ( array(i) > 0 ) then
         l_count = l_count + 1
      endif
   enddo
   deallocate(array)
end function count_mask

end module LIS_domainMod
