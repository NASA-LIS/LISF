!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_domainMod
!BOP
!
! !MODULE: LDT_domainMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage the creation of 
!   the domain (grid, tile, patch and ensemble spaces) used in LDT.
!
!   \subsubsection{Overview}
!   The `domain' is defined in LDT as a cartesian grid in two dimensions
!   with an ordering from the lower left corner to the upper right corner. 
!   This rule is followed for the running domain of LDT that may use different
!   map projections. 
!  
!   When the running domain is created, LDT generates a 1-d vector 
!   representation of the 2-d input grid based on the input landmask (where
!   open water points are typically excluded, unless a lake or open water surface
!   types are included). The basic computational unit in LDT is called 
!   a ``tile''. If options for subgrid variability are employed, LDT will
!   create multiple number of tiles per each grid cell. 
!   The subgrid tiles can be defined based on the distribution of 
!   vegetation types, soil characteristics or topography. Each tile can 
!   further define a specified number of ensembles.
!
!   LDT also supports the use of multiple surface types. For example, 
!   the domain could consist of `land' surfaces and `lake' surfaces, 
!   which can be executed over land points and inland water points, 
!   respectively.  Within LDT, these surface models are defined to 
!   run over ``patches''. The patches across different surface models 
!   sum to the tile dimension. 
!    
! 
! !REVISION HISTORY: 
!  17 Feb 2004    Sujay Kumar  Initial Specification
!  24 Aug 2008    Sujay Kumar  Implemented halo support 
!   3 Aug 2012    Sujay Kumar  Added support for flexible tiling
!  25 Aug 2015    KR Arsenault Added in support for parallel processing
! 
  use ESMF
  use LDT_coreMod
  use LDT_paramDataMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_mpiMod
  use map_utils

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_domainInit               ! initialize specified domains
  public :: LDT_timeInit
  public :: LDT_domain_setup
  public :: LDT_quilt_domain             ! generate quilted domains
!  public :: LDT_domain_finalize          ! cleanup allocated structures
  public :: LDT_setParameterDomainSpecs  ! sets map projection information for parameters
  public :: LDT_setDomainSpecs
  public :: isSurfaceTypeSelected
!EOP

contains

!BOP
! !ROUTINE: LDT_domainInit
! \label{LDT_domainInit}
!
! !INTERFACE: 
  subroutine LDT_domainInit()

! !USES:
   !NONE
!
! !DESCRIPTION: 
!
!  This routine invokes the registry that defines the domain implementations, 
!  first. This is followed by calling the routines to read the runtime 
!  specific parameters for each domain instance. A call to create the domain
!  is issued next, which is expected to generate the grid, tile, patch
!  and ensemble spaces for each domain. This routine is also expected to perform
!  any domain decomposition needed for parallel processing. Finally, the domain
!  decomposition information is disseminated to individual processors. 
!
! The calling sequence is:
! \begin{description}
!  \item[LDT\_domain\_plugin] (\ref{LDT_domain_plugin}) \newline
!    sets up function table registries for implemented domains
!  \item[readinput] (\ref{readinput}) \newline
!    invokes the generic method in the registry to read domain specific
!    input options from the configuration file
!  \item[make\_domain] (\ref{make_domain}) \newline
!    issues the call to generate the tile, patch, and grid dimensions, 
!    and associated data structures. 
!  \end{description}
!
!EOP
    integer :: ierr
    integer :: i 
    integer :: n
    integer :: k
    integer :: m
    integer :: gindex, ntiles
    integer, allocatable :: deblklist(:,:,:)
    integer :: stid, enid
    integer :: status
    
    type(ESMF_DistGrid) :: tileDG
    type(ESMF_DistGrid) :: gridDG
    type(ESMF_DistGrid) :: patchDG(LDT_rc%max_model_types)
    type(ESMF_DistGrid) :: gridEnsDG

    integer             :: ntiless
    integer             :: ngrids
    integer             :: tdeltas
    integer             :: gdeltas
    integer             :: npatches(LDT_rc%max_model_types)
    integer             :: patch_deltas(LDT_rc%max_model_types)
!    integer             :: odeltas(LDT_rc%ndas)

! _______________________________________________________
    
    allocate(LDT_rc%nensem(LDT_rc%nnest))
    LDT_rc%nensem(:) = 1
    call ESMF_ConfigFindLabel(LDT_config,"Number of ensembles per tile:",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%nensem(n),rc=status)
       call LDT_verify(status,'Number of ensembles per tile: not defined')
    enddo
       
    allocate(LDT_rc%npatch(LDT_rc%nnest, LDT_rc%max_model_types))
    
    allocate(LDT_rc%surface_maxt(LDT_rc%nnest))
    allocate(LDT_rc%surface_minp(LDT_rc%nnest))
    allocate(LDT_rc%soilt_maxt(LDT_rc%nnest))
    allocate(LDT_rc%soilt_minp(LDT_rc%nnest))
    allocate(LDT_rc%soilf_maxt(LDT_rc%nnest))
    allocate(LDT_rc%soilf_minp(LDT_rc%nnest))
    allocate(LDT_rc%elev_maxt(LDT_rc%nnest))
    allocate(LDT_rc%elev_minp(LDT_rc%nnest))
    allocate(LDT_rc%slope_maxt(LDT_rc%nnest))
    allocate(LDT_rc%slope_minp(LDT_rc%nnest))
    allocate(LDT_rc%aspect_maxt(LDT_rc%nnest))
    allocate(LDT_rc%aspect_minp(LDT_rc%nnest))
    
    call ESMF_ConfigFindLabel(LDT_config,"Maximum number of surface type tiles per grid:",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%surface_maxt(n),rc=status)
       call LDT_verify(status,'Maximum number of surface type tiles per grid: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Minimum cutoff percentage (surface type tiles):",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%surface_minp(n),rc=status)
       call LDT_verify(status,'Minimum cutoff percentage (surface type tiles): not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Maximum number of soil texture tiles per grid:",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soilt_maxt(n),rc=status)
       call LDT_verify(status,'Maximum number of soil texture tiles per grid: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Minimum cutoff percentage (soil texture tiles):",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soilt_minp(n),rc=status)
       call LDT_verify(status,'Minimum cutoff percentage (soil texture tiles): not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Maximum number of soil fraction tiles per grid:",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soilf_maxt(n),rc=status)
       call LDT_verify(status,'Maximum number of soil fraction tiles per grid: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LDT_config, "Minimum cutoff percentage (soil fraction tiles):",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soilf_minp(n),rc=status)
       call LDT_verify(status,'Minimum cutoff percentage (soil fraction tiles): not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Maximum number of elevation bands per grid:",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elev_maxt(n),rc=status)
       call LDT_verify(status,'Maximum number of elevation bands per grid: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Minimum cutoff percentage (elevation bands):",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elev_minp(n),rc=status)
       call LDT_verify(status,'Minimum cutoff percentage (elevation bands): not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Maximum number of slope bands per grid:",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%slope_maxt(n),rc=status)
       call LDT_verify(status,'Maximum number of slope bands per grid: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "Minimum cutoff percentage (slope bands):",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%slope_minp(n),rc=status)
       call LDT_verify(status,'Minimum cutoff percentage (slope bands): not defined')
    enddo


    call ESMF_ConfigFindLabel(LDT_config, "Maximum number of aspect bands per grid:",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%aspect_maxt(n),rc=status)
       call LDT_verify(status,'Maximum number of aspect bands per grid: not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LDT_config, "Minimum cutoff percentage (aspect bands):",rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%aspect_minp(n),rc=status)
       call LDT_verify(status,'Minimum cutoff percentage (aspect bands): not defined')
    enddo

    
!-- Create vegetation tile and grid space:
    call make_domain()


!-- Temporary addition needed for writing output diagnostics(LDT_diagnoseOutputVar):
    do n=1,LDT_rc%nnest
       do k = 1, LDT_rc%ngrid(n)
          allocate(LDT_domain(n)%grid(k)%subgrid_tiles(&
               LDT_rc%surface_maxt(n)*&
               LDT_rc%soilt_maxt(n)*&
               LDT_rc%elev_maxt(n)*&
               LDT_rc%slope_maxt(n)*&
               LDT_rc%aspect_maxt(n)*&
               LDT_rc%nensem(n)))
          LDT_domain(n)%grid(k)%ntiles = 0
       enddo
       do k = 1, LDT_rc%ntiles(n)
          gindex = LDT_domain(n)%tile(k)%index
          LDT_domain(n)%grid(gindex)%ntiles = &
             LDT_domain(n)%grid(gindex)%ntiles + 1
          ntiles = LDT_domain(n)%grid(gindex)%ntiles
          LDT_domain(n)%grid(gindex)%subgrid_tiles(ntiles) = k
       enddo
    enddo

    
!-- Begin decomposition of tiles and gridcells for MPI calls:
    do n=1,LDT_rc%nnest
       LDT_ntiless(n,LDT_localPet) = LDT_rc%ntiles(n)
       LDT_tdeltas(n,LDT_localPet) = LDT_rc%ntiles(n)

       do m=1,LDT_rc%max_model_types
          LDT_npatches(n,m,LDT_localPet) = LDT_rc%npatch(n,m)
          LDT_patch_deltas(n,m,LDT_localPet) = LDT_rc%npatch(n,m)
       enddo

       LDT_ngrids(n,LDT_localPet)  = LDT_rc%ngrid(n)
       LDT_gdeltas(n,LDT_localPet) = LDT_rc%ngrid(n)

!-----------------------------------------------------------------------------
!  The grid, tile space sizes of the decomposed domains are gathered 
!  from each processor to compute the total sizes of the entire domain. 
!-----------------------------------------------------------------------------
#if (defined SPMD)
       call MPI_ALLGATHER(LDT_ntiless(n,LDT_localPet),1,MPI_INTEGER,&
            LDT_ntiless(n,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(LDT_ngrids(n,LDT_localPet),1,MPI_INTEGER,&
            LDT_ngrids(n,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(LDT_tdeltas(n,LDT_localPet),1,MPI_INTEGER,&
            LDT_tdeltas(n,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(LDT_gdeltas(n,LDT_localPet),1,MPI_INTEGER,&
            LDT_gdeltas(n,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)

       do m=1,LDT_rc%max_model_types
          call MPI_ALLGATHER(LDT_npatches(n,m,LDT_localPet),1,MPI_INTEGER,&
               LDT_npatches(n,m,:),1,MPI_INTEGER, &
               MPI_COMM_WORLD,ierr)
       enddo
       do m=1,LDT_rc%max_model_types
          call MPI_ALLGATHER(LDT_patch_deltas(n,m,LDT_localPet),1,MPI_INTEGER,&
               LDT_patch_deltas(n,m,:),1,MPI_INTEGER, &
               MPI_COMM_WORLD,ierr)
       enddo
#endif

       LDT_rc%glbntiles(n) = 0 
       do i=0,LDT_npes-1
          LDT_rc%glbntiles = LDT_rc%glbntiles(n) + LDT_ntiless(n,i)
       enddo
       
       LDT_rc%glbngrid(n) = 0 
       do i=0,LDT_npes-1
          LDT_rc%glbngrid(n) = LDT_rc%glbngrid(n) + LDT_ngrids(n,i)
       enddo

       do m=1,LDT_rc%max_model_types
          LDT_rc%glbnpatch(n,m) = 0 
          do i=0,LDT_npes-1
             LDT_rc%glbnpatch(n,m) = LDT_rc%glbnpatch(n,m)+ LDT_npatches(n,m,i)
          enddo
       enddo
       
       if(LDT_masterproc) then 
          LDT_toffsets(n,0) = 0 
          do i=1,LDT_npes-1
             LDT_toffsets(n,i) = LDT_toffsets(n,i-1)+LDT_tdeltas(n,i-1)
          enddo
          LDT_goffsets(n,0) = 0 
          do i=1,LDT_npes-1
             LDT_goffsets(n,i) = LDT_goffsets(n,i-1)+LDT_gdeltas(n,i-1)
          enddo
          
          LDT_patch_offsets(n,:,0) =0
          do i=1,LDT_npes-1
             do m=1,LDT_rc%max_model_types
                LDT_patch_offsets(n,m,i) = LDT_patch_offsets(n,m,i-1)+&
                     LDT_patch_deltas(n,m,i-1)
             enddo
          enddo
       end if
       
#if (defined SPMD)
       call MPI_BCAST(LDT_toffsets(n,:), LDT_npes, MPI_INTEGER,0, &
            MPI_COMM_WORLD, ierr)
       call MPI_BCAST(LDT_goffsets(n,:), LDT_npes, MPI_INTEGER,0, &   ! updated for parallel
            MPI_COMM_WORLD, ierr)
       do m=1,LDT_rc%max_model_types
          call MPI_BCAST(LDT_patch_offsets(n,m,:), LDT_npes, MPI_INTEGER,0, &
               MPI_COMM_WORLD, ierr)
       enddo
#endif
       
       allocate(deblklist(1,2,LDT_npes))      
       do i=0,LDT_npes-1
          stid = LDT_toffsets(n,i)+1
          enid = stid + LDT_ntiless(n,i)-1

          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)          
       enddo

       tileDG = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/LDT_rc%glbntiles(n)/),&
            deBlockList=deblklist,rc=status)
       call LDT_verify(status)

       do i=0,LDT_npes-1
          stid = LDT_goffsets(n,i)+1
          enid = stid + LDT_ngrids(n,i)-1

          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)
       enddo

       gridDG = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/LDT_rc%glbngrid(n)/),&
            deBlockList=deblklist,rc=status)
       call LDT_verify(status)

       do m=1,LDT_rc%max_model_types

          do i=0,LDT_npes-1
             stid = LDT_patch_offsets(n,m,i)+1
             enid = stid + LDT_npatches(n,m,i)-1
             
             deblklist(:,1,i+1) = (/stid/)
             deblklist(:,2,i+1) = (/enid/)
          enddo

          patchDG(m) = ESMF_DistGridCreate(minIndex=(/1/),&
               maxIndex=(/LDT_rc%glbnpatch(n,m)/),&
               deBlockList=deblklist,rc=status)
          call LDT_verify(status)
       enddo

       deallocate(deblklist)

       allocate(deblklist(2,2,LDT_npes))
       do i=0,LDT_npes-1
          stid = LDT_goffsets(n,i)+1
          enid = stid+LDT_ngrids(n,i)-1
          deblklist(:,1,i+1) = (/stid,1/)
          deblklist(:,2,i+1) = (/enid,LDT_rc%nensem(n)/)
       enddo

       gridEnsDG = ESMF_DistGridCreate(minIndex=(/1,1/),&
            maxIndex=(/LDT_rc%glbngrid(n),LDT_rc%nensem(n)/),&
            deBlockList=deblklist,rc=status)
       call LDT_verify(status)

       LDT_vecTile(n) = ESMF_GridCreate(name = "LDT Tile Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = tileDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LDT_verify(status)

       LDT_vecGrid(n) = ESMF_GridCreate(name = "LDT Grid Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = gridDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LDT_verify(status)


       LDT_ensOngrid(n) = ESMF_GridCreate(name = "LDT Grid Ensemble Space",&
            coordTypeKind = ESMF_TYPEKIND_R4, distGrid = gridEnsDG,&
            gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,0/),rc=status)
       call LDT_verify(status)

       do m=1,LDT_rc%max_model_types
          LDT_vecPatch(n,m) = ESMF_GridCreate(name="LDT Patch Space",&
               coordTypeKind=ESMF_TYPEKIND_R4, distGrid = patchDG(m),&
               gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
          call LDT_verify(status)
       enddo
       deallocate(deblklist)

       allocate(LDT_domain(n)%datamask(LDT_rc%lnc(n), LDT_rc%lnr(n)))
       LDT_domain(n)%datamask = 1

    enddo

end subroutine LDT_domainInit

!BOP
! !ROUTINE: LDT_timeInit
! \label{LDT_timeInit}
!
! !INTERFACE: 
  subroutine LDT_timeInit()

! !USES:
!
! !DESCRIPTION: 
!
!  This routine reads in LDT config starting and ending 
!   date and time entries, and also the LIS output timestep
!   and temporal averaging interval entries.
!
!EOP
    integer :: i 
    integer :: n
    integer :: status
    character*20 :: stime
! _______________________________________________________

    ! Read in ldt.config Starting date and time:
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%syr,label="Starting year:",&
         rc=status)
    call LDT_verify(status,'Starting year: not defined')
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%smo,label="Starting month:",&
         rc=status)
    call LDT_verify(status,'Starting month: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%sda,label="Starting day:",&
         rc=status)
    call LDT_verify(status,'Starting day: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%shr,label="Starting hour:",&
         rc=status)
    call LDT_verify(status,'Starting hour: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%smn,label="Starting minute:",&
         rc=status)
    call LDT_verify(status,'Starting minute: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%sss,label="Starting second:",&
         rc=status)
    call LDT_verify(status,'Starting second: not defined')

    ! Read in ldt.config Ending date and time:
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%eyr1,label="Ending year:",&
         rc=status)
    call LDT_verify(status,'Ending year: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%emo1,label="Ending month:",&
         rc=status)
    call LDT_verify(status,'Ending month: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%eda1,label="Ending day:",&
         rc=status)
    call LDT_verify(status,'Ending day: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%ehr1,label="Ending hour:",&
         rc=status)
    call LDT_verify(status,'Ending hour: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%emn1,label="Ending minute:",&
         rc=status)
    call LDT_verify(status,'Ending minute: not defined')
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%ess1,label="Ending second:",&
         rc=status)
    call LDT_verify(status,'Ending second: not defined')
    
    LDT_rc%eyr = LDT_rc%eyr1
    LDT_rc%emo = LDT_rc%emo1
    LDT_rc%eda = LDT_rc%eda1
    LDT_rc%ehr = LDT_rc%ehr1
    LDT_rc%emn = LDT_rc%emn1
    LDT_rc%ess = LDT_rc%ess1
    
    LDT_rc%YR= LDT_rc%SYR
    LDT_rc%MO= LDT_rc%SMO
    LDT_rc%DA= LDT_rc%SDA
    LDT_rc%HR= LDT_rc%SHR
    LDT_rc%MN= LDT_rc%SMN
    LDT_rc%SS= LDT_rc%SSS
    
    call ESMF_ConfigFindLabel(LDT_config,"LIS output timestep:",rc=status)
    call LDT_warning(status,'*')
    call LDT_warning(status,'* WARNING:: LIS output timestep: not defined')
    call LDT_warning(status,'*')
    call ESMF_ConfigGetAttribute(LDT_config,stime,rc=status,default="1da")

    call LDT_parseTimeString(stime, LDT_rc%ts)
    LDT_rc%nts(:) = LDT_rc%ts

    ! Call LDT time manager init routine:
    call LDT_timemgr_init(LDT_rc)

    call ESMF_ConfigFindLabel(LDT_config,"Temporal averaging interval:",rc=status)
    call LDT_warning(status,'*')
    call LDT_warning(status,'* WARNING:: Temporal averaging interval: not defined')
    call LDT_warning(status,'*')
    call ESMF_ConfigGetAttribute(LDT_config,stime,&
         rc=status, default="1da")

    call LDT_parseTimeString(stime, LDT_rc%tavgInterval)

end subroutine LDT_timeInit


!BOP
! !ROUTINE: LDT_setDomainSpecs
! \label{LDT_setDomainSpecs}
!
! !INTERFACE: 
  subroutine LDT_setDomainSpecs()

! !USES:
   use LDT_domain_pluginMod, only : LDT_domain_plugin
!
! !DESCRIPTION: 
!
!  This routine reads in LDT config run domain and projection
!   entries.
!
!EOP
    integer   :: n, c,r,gindex
    integer   :: bufsize


    call LDT_domain_plugin

    do n=1,LDT_rc%nnest

       ! Read in the domain's config file entries:
       call readinput(trim(LDT_rc%lis_map_proj(n))//char(0),n)

       allocate(LDT_domain(n)%lat(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LDT_domain(n)%lon(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       ! Assign the LDT run domain lat and lon point values:
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             gindex = c+(r-1)*LDT_rc%lnc(n)
             call ij_to_latlon(LDT_domain(n)%ldtproj,&
                  real(c), real(r), LDT_domain(n)%lat(gindex),&
                  LDT_domain(n)%lon(gindex))     
          enddo
       enddo

       ! Set a buffer window around main run-domain;
       !  (used for budget-bilinear interpolation processes) 
       bufsize = 4
       
       LDT_rc%gnc_buf(n) = LDT_rc%gnc(n) + bufsize
       LDT_rc%gnr_buf(n) = LDT_rc%gnr(n) + bufsize

       LDT_rc%lnc_b(n) = LDT_rc%gnc_buf(n) 
       LDT_rc%lnr_b(n) = LDT_rc%gnr_buf(n) 
       
       allocate(LDT_domain(n)%lat_b(&
            LDT_rc%gnc_buf(n)*LDT_rc%gnr_buf(n)))
       allocate(LDT_domain(n)%lon_b(&
            LDT_rc%gnc_buf(n)*LDT_rc%gnr_buf(n)))
       
       do r=-1,LDT_rc%gnr(n)+2
          do c=-1,LDT_rc%gnc(n)+2
             gindex = (c+2) + (r+2 - 1)*LDT_rc%gnc_buf(n)
             call ij_to_latlon(LDT_domain(n)%ldtproj,&
                  real(c), real(r), LDT_domain(n)%lat_b(gindex),&
                  LDT_domain(n)%lon_b(gindex))
          enddo
       enddo
       
    enddo

  end subroutine LDT_setDomainSpecs


!BOP
! !ROUTINE: LDT_domain_setup
! \label{LDT_domain_setup}
! 
! !INTERFACE: 
  subroutine LDT_domain_setup(n)

! !USES: 

! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine computes the variables used to map between the tile, grid, 
!  and patch spaces in LDT. The global indices including the halo regions 
!  are also computed. 
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
    integer              :: ntiles, npatch, stid

! Gathering the number of tiles per grid to the master processor, in the
! right order:

    allocate(LDT_domain(n)%ntiles_pergrid(LDT_rc%gnc(n)*LDT_rc%gnr(n)))
    allocate(LDT_domain(n)%str_tind(LDT_rc%gnc(n)*LDT_rc%gnr(n)))
    allocate(ntiles_pergrid(LDT_rc%lnc(n)*LDT_rc%lnr(n)),stat=ierr)
    allocate(ntiles_pergrid_red(LDT_rc%lnc_red(n)*LDT_rc%lnr_red(n)),stat=ierr)

    allocate(npatch_pergrid(LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%max_model_types))
    allocate(npatch_pergrid_red(LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%max_model_types))
    do m=1,LDT_rc%max_model_types
       allocate(LDT_surface(n,m)%npatch_pergrid(LDT_rc%gnc(n)*LDT_rc%gnr(n)))
       allocate(LDT_surface(n,m)%str_patch_ind(LDT_rc%gnc(n)*LDT_rc%gnr(n)))
    enddo

    do t=1,LDT_rc%ntiles(n)
       LDT_domain(n)%tile(t)%pens = 1.0/float(LDT_rc%nensem(n))      
    enddo
    do m=1,LDT_rc%max_model_types
       do t=1,LDT_rc%npatch(n,m)
          LDT_surface(n,m)%tile(t)%pens = 1.0/float(LDT_rc%nensem(n))
       enddo
    enddo

    do t=1,LDT_rc%lnc(n)*LDT_rc%lnr(n)
       ntiles_pergrid(t) = 0 
    enddo
    do t=1,LDT_rc%ntiles(n)
       index = LDT_domain(n)%gindex(LDT_domain(n)%tile(t)%col,&
            LDT_domain(n)%tile(t)%row)
       if(index.ne.-1) then 
          gid = LDT_domain(n)%tile(t)%col+&
               (LDT_domain(n)%tile(t)%row-1)*LDT_rc%lnc(n)
          ntiles_pergrid(gid) = ntiles_pergrid(gid)+1
       endif
    enddo

    do r=LDT_nss_halo_ind(n,LDT_localPet+1), LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1), LDT_ewe_halo_ind(n,LDT_localPet+1)
          c1 = c-LDT_ews_halo_ind(n,LDT_localPet+1)+1
          r1 = r-LDT_nss_halo_ind(n,LDT_localPet+1)+1
          gid = c1+(r1-1)*LDT_rc%lnc(n)
          
          if(r.ge.LDT_nss_ind(n,LDT_localPet+1).and. &
               r.le.LDT_nse_ind(n,LDT_localPet+1).and.&
               c.ge.LDT_ews_ind(n,LDT_localPet+1).and.&
               c.le.LDT_ewe_ind(n,LDT_localPet+1))then !points not in halo
             c2 = c-LDT_ews_ind(n,LDT_localPet+1)+1
             r2 = r-LDT_nss_ind(n,LDT_localPet+1)+1
             gid_red = c2+(r2-1)*LDT_rc%lnc_red(n)
             ntiles_pergrid_red(gid_red) = ntiles_pergrid(gid)
          endif
       enddo
    enddo

    deallocate(ntiles_pergrid)


    if(LDT_masterproc) then 
       allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
       allocate(gtmp1(LDT_rc%gnc(n)*LDT_rc%gnr(n)))
    else
       allocate(gtmp1(1))
    endif

#if (defined SPMD)
    call MPI_GATHERV(ntiles_pergrid_red,&
         LDT_deltas(n,LDT_localPet),MPI_INTEGER,gtmp1,&
         LDT_deltas(n,:),LDT_offsets(n,:),&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
#else
    gtmp1 = ntiles_pergrid_red
#endif

    deallocate(ntiles_pergrid_red)

    if(LDT_masterproc) then 
       count=1
       do l=1,LDT_npes
          do r=LDT_nss_ind(n,l), LDT_nse_ind(n,l)
             do c=LDT_ews_ind(n,l), LDT_ewe_ind(n,l)
                gtmp(c,r) = gtmp1(count)
                count = count+1
             enddo
          enddo
       enddo

       count=1
       do r=1,LDT_rc%gnr(n)
          do c=1,LDT_rc%gnc(n)
             LDT_domain(n)%ntiles_pergrid(count) = gtmp(c,r)
             count = count+1
          enddo
       enddo
       
       deallocate(gtmp)
       deallocate(gtmp1)
       
      
       if(LDT_domain(n)%ntiles_pergrid(1).ge.0) then 
          LDT_domain(n)%str_tind(1) = 1
       else
          LDT_domain(n)%str_tind(1) = 0 
       endif
       
       do t=2,LDT_rc%gnc(n)*LDT_rc%gnr(n)
          LDT_domain(n)%str_tind(t) = LDT_domain(n)%str_tind(t-1)+&
               LDT_domain(n)%ntiles_pergrid(t-1)
       enddo       
    else
       deallocate(gtmp1)
    endif
#if (defined SPMD)
    call MPI_Bcast(LDT_domain(n)%ntiles_pergrid,LDT_rc%gnc(n)*LDT_rc%gnr(n),&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(LDT_domain(n)%str_tind,LDT_rc%gnc(n)*LDT_rc%gnr(n),&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif   

    if(LDT_masterproc) then 
       LDT_rc%glbngrid_red(n)=0
       LDT_rc%glbntiles_red(n)=0
       do l=1,LDT_npes
          do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
             do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                if(r.ge.LDT_nss_ind(n,l).and.&
                     r.le.LDT_nse_ind(n,l).and.&
                     c.ge.LDT_ews_ind(n,l).and.&
                     c.le.LDT_ewe_ind(n,l))then !points not in halo
                   gid = c+(r-1)*LDT_rc%gnc(n)
                   ntiles = LDT_domain(n)%ntiles_pergrid(gid)
                   stid = LDT_domain(n)%str_tind(gid)
                   if ( ntiles .ne. 0 ) then
                      LDT_rc%glbngrid_red(n) = LDT_rc%glbngrid_red(n) + 1
                   endif
                   do t=1,ntiles
                      LDT_rc%glbntiles_red(n) = LDT_rc%glbntiles_red(n) + 1
                   enddo
                endif
             enddo
          enddo
       enddo
    endif    
#if (defined SPMD)
    call MPI_Bcast(LDT_rc%glbntiles_red(n),1,&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(LDT_rc%glbngrid_red(n),1,&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

    do t=1,LDT_rc%lnc(n)*LDT_rc%lnr(n)
       npatch_pergrid(t,:) = 0
    enddo

    do m=1,LDT_rc%max_model_types
       do t=1,LDT_rc%npatch(n,m)
          index =  LDT_domain(n)%gindex(LDT_surface(n,m)%tile(t)%col,&
               LDT_surface(n,m)%tile(t)%row)
          if(index.ne.-1) then 
             gid = LDT_surface(n,m)%tile(t)%col + & 
                  (LDT_surface(n,m)%tile(t)%row-1)*LDT_rc%lnc(n)
             npatch_pergrid(gid,m) = npatch_pergrid(gid,m) + 1
          endif
       enddo
    enddo

    do r=LDT_nss_halo_ind(n,LDT_localPet+1), LDT_nse_halo_ind(n,LDT_localPet+1)
       do c=LDT_ews_halo_ind(n,LDT_localPet+1), LDT_ewe_halo_ind(n,LDT_localPet+1)
          c1 = c-LDT_ews_halo_ind(n,LDT_localPet+1)+1
          r1 = r-LDT_nss_halo_ind(n,LDT_localPet+1)+1
          gid = c1+(r1-1)*LDT_rc%lnc(n)
          
          if(r.ge.LDT_nss_ind(n,LDT_localPet+1).and. &
               r.le.LDT_nse_ind(n,LDT_localPet+1).and.&
               c.ge.LDT_ews_ind(n,LDT_localPet+1).and.&
               c.le.LDT_ewe_ind(n,LDT_localPet+1))then !points not in halo
             c2 = c-LDT_ews_ind(n,LDT_localPet+1)+1
             r2 = r-LDT_nss_ind(n,LDT_localPet+1)+1
             gid_red = c2+(r2-1)*LDT_rc%lnc_red(n)
             npatch_pergrid_red(gid_red,:) = npatch_pergrid(gid,:)
          endif
       enddo
    enddo

    do m=1,LDT_rc%max_model_types
       if(LDT_masterproc) then 
          allocate(gtmp(LDT_rc%gnc(n),LDT_rc%gnr(n)))
          allocate(gtmp1(LDT_rc%gnc(n)*LDT_rc%gnr(n)))
       else
          allocate(gtmp1(1))
       endif
       
#if (defined SPMD)
       call MPI_GATHERV(npatch_pergrid_red(:,m),LDT_deltas(n,LDT_localPet),MPI_INTEGER,gtmp1,&
            LDT_deltas(n,:),LDT_offsets(n,:),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
       gtmp1 = npatch_pergrid_red(:,m)
#endif

       if(LDT_masterproc) then 
          count=1
          do l=1,LDT_npes
             do r=LDT_nss_ind(n,l), LDT_nse_ind(n,l)
                do c=LDT_ews_ind(n,l), LDT_ewe_ind(n,l)
                   gtmp(c,r) = gtmp1(count)
                   count = count+1
                enddo
             enddo
          enddo
          
          count=1
          do r=1,LDT_rc%gnr(n)
             do c=1,LDT_rc%gnc(n)
                LDT_surface(n,m)%npatch_pergrid(count) = gtmp(c,r)
                count = count+1
             enddo
          enddo
          
          deallocate(gtmp)
          deallocate(gtmp1)
       
          
          if(LDT_surface(n,m)%npatch_pergrid(1).ge.0) then 
             LDT_surface(n,m)%str_patch_ind(1) = 1
          else
             LDT_surface(n,m)%str_patch_ind(1) = 0 
          endif
          
          do t=2,LDT_rc%gnc(n)*LDT_rc%gnr(n)
             LDT_surface(n,m)%str_patch_ind(t) = LDT_surface(n,m)%str_patch_ind(t-1)+&
                  LDT_surface(n,m)%npatch_pergrid(t-1)
          enddo
       else
          deallocate(gtmp1)
       endif
#if (defined SPMD)
       call MPI_Bcast(LDT_surface(n,m)%npatch_pergrid,LDT_rc%gnc(n)*LDT_rc%gnr(n),&
            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(LDT_surface(n,m)%str_patch_ind,LDT_rc%gnc(n)*LDT_rc%gnr(n),&
            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif   

       if(LDT_masterproc) then 
          LDT_rc%glbnpatch_red(n,m)=0
          do l=1,LDT_npes
             do r=LDT_nss_halo_ind(n,l),LDT_nse_halo_ind(n,l)
                do c=LDT_ews_halo_ind(n,l),LDT_ewe_halo_ind(n,l)
                   if(r.ge.LDT_nss_ind(n,l).and.&
                        r.le.LDT_nse_ind(n,l).and.&
                        c.ge.LDT_ews_ind(n,l).and.&
                        c.le.LDT_ewe_ind(n,l))then !points not in halo
                      gid = c+(r-1)*LDT_rc%gnc(n)
                      npatch = LDT_surface(n,m)%npatch_pergrid(gid)
                      stid = LDT_surface(n,m)%str_patch_ind(gid)
                      do t=1,npatch
                         LDT_rc%glbnpatch_red(n,m) = &
                              LDT_rc%glbnpatch_red(n,m) + 1
                      enddo
                   endif
                enddo
             enddo
          enddo
       endif
#if (defined SPMD)
       call MPI_Bcast(LDT_rc%glbnpatch_red(n,m),1,&
            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

    enddo
    deallocate(npatch_pergrid)
    deallocate(npatch_pergrid_red)

  end subroutine LDT_domain_setup

!BOP
! !ROUTINE: LDT_quilt_domain
! \label{LDT_quilt_domain}
! 
! !INTERFACE: 
  subroutine LDT_quilt_domain( n, nc, nr )

! !USES:     

! !ARGUMENTS: 
    integer, intent(in)  :: n 
    integer, intent(in)  :: nc
    integer, intent(in)  :: nr
! 
! !DESCRIPTION: 
! This routine generates the quilted domain extents and sizes based on the 
! processor layout and the size of the specified halos. 
! 
!EOP

    integer    :: ips, ipe, jps, jpe
    integer    :: Px, Py, P
    integer    :: mytask_x, mytask_y
    integer    :: i, j,ierr

    integer    :: deltas
    integer    :: ews_halo_ind
    integer    :: ewe_halo_ind
    integer    :: nss_halo_ind
    integer    :: nse_halo_ind
    integer    :: ews_ind
    integer    :: ewe_ind
    integer    :: nss_ind
    integer    :: nse_ind


    mytask_x = mod( LDT_localPet, LDT_rc%npesx )
    mytask_y = LDT_localPet / LDT_rc%npesx

    j = 1
    ips = -1
    do i=1, nc+1
       call LDT_mpDecomp(i,j,1,nc+1,1,nr+1,LDT_rc%npesx, LDT_rc%npesy,Px,Py,P)
       if(Px.eq.mytask_x) then 
          ipe = i 
          if(ips.eq.-1) then
            ips = i 
          endif
       endif
    enddo
     
    i = 1
    jps = -1
    do j=1, nr+1
       call LDT_mpDecomp(i,j,1,nc+1,1,nr+1,LDT_rc%npesx,LDT_rc%npesy, Px, Py, P)
       if(Py.eq.mytask_y) then 
          jpe = j
          if(jps.eq.-1) then 
            jps = j 
          endif
       endif
    enddo

! Original LDT code:
#if 0
     LDT_ews_halo_ind(n,LDT_localPet+1) = max(ips-LDT_rc%halox, 1)   ! Switched haloy and halox 
     LDT_ewe_halo_ind(n,LDT_localPet+1) = min(ipe+LDT_rc%halox, nc)
     LDT_nss_halo_ind(n,LDT_localPet+1) = max(jps-LDT_rc%haloy, 1)
     LDT_nse_halo_ind(n,LDT_localPet+1) = min(jpe+LDT_rc%haloy, nr)
     LDT_ews_ind(n,LDT_localPet+1) = ips
     LDT_ewe_ind(n,LDT_localPet+1) = min(ipe, nc)
     LDT_nss_ind(n,LDT_localPet+1) = jps
     LDT_nse_ind(n,LDT_localPet+1) = min(jpe, nr)

     LDT_rc%lnc(n) = LDT_ewe_halo_ind(n,LDT_localPet+1)-&
          LDT_ews_halo_ind(n,LDT_localPet+1)+1
     LDT_rc%lnr(n) = LDT_nse_halo_ind(n,LDT_localPet+1)-&
          LDT_nss_halo_ind(n,LDT_localPet+1)+1

     LDT_rc%lnc_red(n)= (LDT_ewe_ind(n,LDT_localPet+1)-&
          LDT_ews_ind(n,LDT_localPet+1))+1
     LDT_rc%lnr_red(n)= (LDT_nse_ind(n,LDT_localPet+1)-&
          LDT_nss_ind(n,LDT_localPet+1))+1
#endif

! LIS-7 code:
!#if 0
     ews_halo_ind = max(ips-LDT_rc%halox, 1)
     ewe_halo_ind = min(ipe+LDT_rc%halox, nc)
     nss_halo_ind = max(jps-LDT_rc%haloy, 1)
     nse_halo_ind = min(jpe+LDT_rc%haloy, nr)

     ews_ind = ips
     ewe_ind = min(ipe, nc)
     nss_ind = jps
     nse_ind = min(jpe, nr)

     LDT_rc%lnc(n) = ewe_halo_ind-&
          ews_halo_ind+1
     LDT_rc%lnr(n) = nse_halo_ind-&
          nss_halo_ind+1

     LDT_rc%lnc_red(n)= ewe_ind-&
          ews_ind+1
     LDT_rc%lnr_red(n)= nse_ind-&
          nss_ind+1
!#endif

     write(unit=LDT_logunit,fmt=*)'local domain',':(',LDT_rc%lnc(n),LDT_rc%lnr(n),')'
     write(unit=LDT_logunit,fmt=*)'local domain without halo',':(',&
          LDT_rc%lnc_red(n),LDT_rc%lnr_red(n),')'
     write(unit=LDT_logunit,fmt=*)'running domain',':(',nc,nr,')'

! Former LDT version:
!     LDT_deltas(n,LDT_localPet) = LDT_rc%lnc_red(n)*LDT_rc%lnr_red(n)
! LIS-7 version:
     deltas = LDT_rc%lnc_red(n)*LDT_rc%lnr_red(n)

#if (defined SPMD)
! Former LDT code:
#if 0
     call MPI_ALLGATHER(LDT_deltas(n,LDT_localPet),1,MPI_INTEGER,LDT_deltas(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_ews_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_ews_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_nss_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_nss_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_ewe_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_ewe_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_nse_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_nse_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_ews_halo_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_ews_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_nss_halo_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_nss_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_ewe_halo_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_ewe_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LDT_nse_halo_ind(n,LDT_localPet+1),1,MPI_INTEGER,LDT_nse_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
#endif

! LIS-7 code:
!#if 0
     call MPI_ALLGATHER(deltas,1,MPI_INTEGER,LDT_deltas(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(ews_ind,1,MPI_INTEGER,LDT_ews_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(nss_ind,1,MPI_INTEGER,LDT_nss_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(ewe_ind,1,MPI_INTEGER,LDT_ewe_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(nse_ind,1,MPI_INTEGER,LDT_nse_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)

     call MPI_ALLGATHER(ews_halo_ind,1,MPI_INTEGER,LDT_ews_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(nss_halo_ind,1,MPI_INTEGER,LDT_nss_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(ewe_halo_ind,1,MPI_INTEGER,LDT_ewe_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(nse_halo_ind,1,MPI_INTEGER,LDT_nse_halo_ind(n,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
!#endif

#else
     LDT_deltas(n,:) = deltas
     LDT_ews_ind(n,:) = ews_ind
     LDT_nss_ind(n,:) = nss_ind
     LDT_ewe_ind(n,:) = ewe_ind
     LDT_nse_ind(n,:) = nse_ind

     LDT_ews_halo_ind(n,:) = ews_halo_ind
     LDT_nss_halo_ind(n,:) = nss_halo_ind
     LDT_ewe_halo_ind(n,:) = ewe_halo_ind
     LDT_nse_halo_ind(n,:) = nse_halo_ind
#endif   

     if(LDT_masterproc) then 
        LDT_offsets(n,0) = 0 
        do i=1,LDT_npes-1
           LDT_offsets(n,i) = LDT_offsets(n,i-1)+LDT_deltas(n,i-1)
        enddo
     end if
#if (defined SPMD)
     call MPI_BCAST(LDT_offsets(n,:), LDT_npes, MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
#endif

   end subroutine LDT_quilt_domain


!BOP
! 
! !ROUTINE: LDT_setParameterDomainSpecs
! \label{LDT_setParameterDomainSpecs}
! 
! !INTERFACE: 
   subroutine LDT_setParameterDomainSpecs( proj, param_dd, gridDesc )

! !USES: 

! !ARGUMENTS: 
     integer,intent(in)    :: proj
     real,   intent(in)    :: param_dd(20)
     real,   intent(inout) :: gridDesc(20)
! 
! !DESCRIPTION: 
!  This subroutine sets the grid descriptor array section that
!  defines the parameter domain and projection information, based
!  on the map projection used for defining surface parameters
!
!  OBSOLETE CODE!
! 
!EOP      
#if 0
  !- Parameters specified in latlon:
     if(proj.eq.0) then 
        gridDesc(34) = param_dd(1)
        gridDesc(35) = param_dd(2)
        gridDesc(37) = param_dd(3)
        gridDesc(38) = param_dd(4)
        gridDesc(39) = param_dd(6)
        gridDesc(40) = param_dd(5)

        gridDesc(33) = nint((gridDesc(37)-gridDesc(34))&
             /gridDesc(39)) + 1
        gridDesc(32) = nint((gridDesc(38)-gridDesc(35))&
             /gridDesc(40)) + 1

  !- Parameters specified in gaussian:
     elseif(proj.eq.4) then  
        gridDesc(41)  = 4 ! param_dd(1)
        gridDesc(44)  = param_dd(1)
        gridDesc(45)  = param_dd(2)
        gridDesc(47)  = param_dd(3)
        gridDesc(48)  = param_dd(4)
        gridDesc(49)  = param_dd(5)
        gridDesc(50)  = param_dd(6)
        
        gridDesc(42) = nint(360/gridDesc(49))
        gridDesc(43) = 2*gridDesc(50)
        gridDesc(46) = 128

  !- UTM:
     elseif(proj.eq.7) then 
!-  THIS NEEDS TO BE UPDATED (KRA) ::
!        gridDesc(32) = param_dd(4)
!        gridDesc(33) = param_dd(5)
!        gridDesc(34) = param_dd(2)
!        gridDesc(35) = param_dd(3)
!        gridDesc(37) = param_dd(2)+(param_dd(5)-1)*param_dd(6)
!        gridDesc(38) = param_dd(3)+(param_dd(4)-1)*param_dd(6)
!        gridDesc(39) = param_dd(6) !resolution
!        gridDesc(40) = param_dd(1) !zone

        gridDesc(34) = param_dd(1)
        gridDesc(35) = param_dd(2)
        gridDesc(37) = param_dd(3)
        gridDesc(38) = param_dd(4)
        gridDesc(39) = param_dd(5)
        gridDesc(40) = param_dd(6)
! - KRA
 
  !- Projections not currently supported:
     else
        write(*,*)"This parameter projection is not currently supported", proj
        write(*,*)"Stopping program ...."
        call LDT_endrun 
     endif
#endif
   end subroutine LDT_setParameterDomainSpecs

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


    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then 
          soilt_selected = .true. 
       else
          soilt_selected = .false.
       endif      

       if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then 
          soilf_selected = .true. 
       else
          soilf_selected = .false.
       endif

       if(LDT_LSMparam_struc(n)%elevation%selectOpt.eq.1) then
          elev_selected = .true. 
       else
          elev_selected = .false. 
       endif

       if(LDT_LSMparam_struc(n)%slope%selectOpt.eq.1) then
          slope_selected = .true. 
       else
          slope_selected = .false. 
       endif

       if(LDT_LSMparam_struc(n)%aspect%selectOpt.eq.1) then
          aspect_selected = .true. 
       else
          aspect_selected = .false. 
       endif
!-----------------------------------------------------------------------
! normalize the vegetation distributions
!-----------------------------------------------------------------------
       ! EMK...Special handling of single surface tile case.  First find
       ! dominant surface model type in each grid point (LSM, Openwater, etc),
       ! then find the dominant tile for that model type.
       if (LDT_rc%surface_maxt(n) == 1) then          
          call calculate_dominant_sfctile(n, &
               LDT_LSMparam_struc(n)%dommask%value, &
               LDT_LSMparam_struc(n)%sfctype%num_bins, &
               LDT_rc%surface_minp(n), LDT_rc%surface_maxt(n), &
               LDT_LSMparam_struc(n)%sfctype%value, &
               LDT_LSMparam_struc(n)%landcover%value)
       else
          call calculate_domdistribution(n, &
               LDT_LSMparam_struc(n)%sfctype%num_bins,  &
               LDT_rc%surface_minp(n), LDT_rc%surface_maxt(n), &
               LDT_LSMparam_struc(n)%landcover%value)
       endif

       if(soilt_selected) then 
          call calculate_domdistribution(n, &
               LDT_LSMparam_struc(n)%texture%num_bins,&
               LDT_rc%soilt_minp(n), LDT_rc%soilt_maxt(n), &
               LDT_LSMparam_struc(n)%texture%value)
       endif

       if(soilf_selected) then 
          call calculate_domdistribution(n, &
               LDT_LSMparam_struc(n)%sand%num_bins,&
               LDT_rc%soilf_minp(n), LDT_rc%soilf_maxt(n), &
               LDT_LSMparam_struc(n)%sand%value )
       endif

       if(elev_selected) then 
          call calculate_domdistribution(n, &
               LDT_LSMparam_struc(n)%elevation%num_bins,&
               LDT_rc%elev_minp(n), LDT_rc%elev_maxt(n), &
               LDT_LSMparam_struc(n)%elevfgrd%value)
       endif

       if(slope_selected) then 
          call calculate_domdistribution(n, &
               LDT_LSMparam_struc(n)%slope%num_bins,&
               LDT_rc%slope_minp(n), LDT_rc%slope_maxt(n), &
               LDT_LSMparam_struc(n)%slopefgrd%value)
       endif

       if(aspect_selected) then 
          call calculate_domdistribution(n, &
               LDT_LSMparam_struc(n)%aspect%num_bins, &
               LDT_rc%aspect_minp(n), LDT_rc%aspect_maxt(n),&
               LDT_LSMparam_struc(n)%aspectfgrd%value)
       endif

       call create_tilespace(n,soilt_selected, soilf_selected, elev_selected, &
            slope_selected, aspect_selected)
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
     real, intent(in)     :: dommask(LDT_rc%lnc(n), LDT_rc%lnr(n), LDT_LSMparam_struc(n)%landmask%num_bins)
     integer, intent(in)  :: ntypes
     real, intent(in)     :: minp
     integer, intent(in)  :: maxt
     real, intent(in)     :: surfacetypes(LDT_rc%lnc(n), LDT_rc%lnr(n), ntypes)
     real, intent(inout)  :: fgrd(LDT_rc%lnc(n), LDT_rc%lnr(n), ntypes)
     
     ! Locals
     integer :: c, r, t, m, mm, tt
     real    :: maxv
     real    :: model_type_fractions(LDT_rc%max_model_types)

     ! Sanity check
     if (maxt > 1) then
        write(LDT_logunit,*) &
             '[ERR], calculate_dominant_sfctile only works for maxt = 1!'
        call LDT_endrun()
     end if
     
     ! Loop through each patch, figure out the dominant surface model type,
     ! and then find the dominant land cover category for that model type.
     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)

           ! Skip points that are outside of the domain--these have no
           ! land cover values defined.
           if (dommask(c,r,1) .le. 0) cycle

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
           do m = 1, LDT_rc%max_model_types
              ! Skip model type if not used
              if (.not. isSurfaceTypeSelected(LDT_rc%sf_model_type(m))) cycle
              if (model_type_fractions(m) > maxv) then
                 maxv = model_type_fractions(m)
                 mm = m
              end if
           end do ! m
           
           ! Sanity check
           if (mm .eq. 0) then
              write(LDT_logunit,*) &
                   '[ERR] No dominant surface model type found!'
              write(LDT_logunit,*) &
                   'c,r,fgrd: ',c,r,fgrd(c,r,:)
              flush(LDT_logunit)
              call LDT_endrun()
           end if

           ! We now know which surface model type is dominant.
           ! Identify and exclude tiles not used by that model type.
           do t = 1,ntypes
              if (nint(surfacetypes(c,r,t)) .eq. mm) cycle
              fgrd(c,r,t) = 0.0
           end do ! t

           ! At this point, remaining tiles are associated with only a single 
           ! surface model type. Exclude tiles below minimum tile grid area). 
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
              write(LDT_logunit,*) &
                   '[ERR] No surface tiles remain!'
              write(LDT_logunit,*) &
                   'c,r,fgrd: ',c,r,fgrd(c,r,:)
              flush(LDT_logunit)
              call LDT_endrun()
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
    real                 :: fgrd(LDT_rc%lnc(n), LDT_rc%lnr(n), ntypes)
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

    allocate(tsum(LDT_rc%lnc(n), LDT_rc%lnr(n)), stat=ierr)
    call LDT_verify(ierr,'Error allocating tsum.')
    tsum = 0.0
!----------------------------------------------------------------------      
! Exclude tiles with (minimum tile grid area),  
! normalize remaining tiles to 100%
!----------------------------------------------------------------------      
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)            
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
                if(rsum.gt.0.0)fgrd(c,r,t)=fgrd(c,r,t)/rsum
             enddo
             
             rsum=0.0
             do t=1,ntypes
                rsum=rsum+fgrd(c,r,t)
             enddo
             
             if(rsum.lt.0.9999.or.rsum.gt.1.0001)then 
                write(LDT_logunit,*) 'The tile distribution do not sum to 100%'
                call LDT_endrun()
             endif
          endif
       enddo
    enddo

    allocate(pfrac(LDT_rc%lnc(n),LDT_rc%lnr(n),ntypes))
!----------------------------------------------------------------------      
! Exclude tiles with SURFACE_MAXT (Maximum Tiles per grid), 
!   normalize remaining tiles to 100%
! Determine the grid predominance order of the tiles
!  PFRAC(NT) will contain the predominance order of tiles
!----------------------------------------------------------------------      
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n) 
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
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
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
                if(rsum.gt.0.0)fgrd(c,r,t)= fgrd(c,r,t)/rsum
             enddo
             
             rsum=0.0
             do t=1,ntypes
                rsum=rsum+ fgrd(c,r,t)  !recalculate rsum to check 
             enddo
             tsum(c,r)=rsum
           
             if(rsum.lt.0.9999.or.rsum.gt.1.0001)then  !check renormalization
                write(LDT_logunit,*) &
                     'Renormalization failed in calculate_domdistribution'
                call LDT_endrun()
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
       slope_selected, aspect_selected)

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: n
    logical, intent(in)   :: soilt_selected
    logical, intent(in)   :: soilf_selected
    logical, intent(in)   :: elev_selected
    logical, intent(in)   :: slope_selected
    logical, intent(in)   :: aspect_selected

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
!   \item[soilf\_selected]
!    flag to indicate if soil fraction based tiling is used
!   \item[elev\_selected]
!    flag to indicate if elevation based tiling is used
!   \item[slope\_selected]
!    flag to indicate if slope based tiling is used
!   \item[aspect\_selected]
!    flag to indicate if aspect based tiling is used
!  \end{description}
!EOP
    real          :: locallat, locallon
    integer       :: c, r, t, m
    integer       :: iv, it, ie, is, ia
    integer       :: vegt
    integer       :: soilt
    real          :: sand, clay, silt
    real          :: elev, slope, aspect
    integer       :: sf_index
    integer       :: gnc, gnr
    integer       :: kk
    integer       :: kk_sf(LDT_rc%max_model_types)
    real, allocatable :: sumv(:)
    integer       :: ntiles_surface
    integer       :: npatch_surface
    integer       :: ntiles_soilt
    integer       :: ntiles_soilf
    integer       :: ntiles_soil
    integer       :: ntiles_elev
    integer       :: ntiles_slope
    integer       :: ntiles_aspect
    integer       :: gid
    real          :: temp
    integer       :: soilf_index
    integer       :: elev_index
    integer       :: slope_index
    integer       :: aspect_index
! _____________________________________________
    
    gnc = LDT_rc%lnc(n)
    gnr = LDT_rc%lnr(n)
    
    LDT_rc%ntiles(n)=0
    LDT_rc%npatch(n,:) = 0 

    do r=1,LDT_rc%lnr(n)     
       do c=1,LDT_rc%lnc(n)             
          do m=1,LDT_rc%nensem(n)
             if(LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0) then
 
                ntiles_surface = compute_ntiles_surface(n,c,r)
                ntiles_soilt = compute_ntiles_soilt(n,c,r,soilt_selected)
                ntiles_soilf = compute_ntiles_soilf(n,c,r,soilf_selected)
                ntiles_elev = compute_ntiles_elev(n,c,r,elev_selected)
                ntiles_slope = compute_ntiles_slope(n,c,r,slope_selected)
                ntiles_aspect = compute_ntiles_aspect(n,c,r,aspect_selected)

                if(soilf_selected) then 
                   LDT_rc%ntiles(n) = LDT_rc%ntiles(n)+&
                        ntiles_surface* & 
                        ntiles_soilf* &
                        ntiles_elev*&
                        ntiles_slope*& 
                        ntiles_aspect
                else
                   LDT_rc%ntiles(n) = LDT_rc%ntiles(n)+&
                        ntiles_surface* & 
                        ntiles_soilt* &
                        ntiles_elev*&
                        ntiles_slope*& 
                        ntiles_aspect
                endif

                do t=1,LDT_LSMparam_struc(n)%sfctype%num_bins
                   sf_index = nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,t))
                   npatch_surface = compute_npatches_surface(n,c,r,t,sf_index)
                   if(sf_index.ne.0) then
                      if(soilf_selected) then 
                         LDT_rc%npatch(n,sf_index) = & 
                              LDT_rc%npatch(n,sf_index) + &
                              npatch_surface*& 
                              ntiles_soilf* & 
                              ntiles_elev*& 
                              ntiles_slope*& 
                              ntiles_aspect                         
                      else
                         LDT_rc%npatch(n,sf_index) = & 
                              LDT_rc%npatch(n,sf_index) + &
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

    write(unit=LDT_logunit,fmt=*) &
         'Total number of surface tiles',LDT_rc%ntiles(n)
    do m=1,LDT_rc%max_model_types
       if(isSurfaceTypeSelected(LDT_rc%sf_model_type(m))) then 
          write(unit=LDT_logunit,fmt=*) &
               'Total number of tiles for ',&
               trim(LDT_rc%sf_model_type_name(m)),&
               LDT_rc%npatch(n,m)
       endif
    enddo

    LDT_rc%ngrid(n)=0
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)

        ! Determine number of grid points, based on landmask/cover points present:
          if(LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0.and.&
               abs(sum(LDT_LSMparam_struc(n)%landcover%value(c,r,:))-1.0).lt.0.001.and.&
               ifContainsSelectedSurfaceType(n,&
               nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,:)))) then 
             LDT_rc%ngrid(n) = LDT_rc%ngrid(n)+1             

        ! Check if landcover fractions < 1.0 :
          elseif(LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0.and.&
               sum(LDT_LSMparam_struc(n)%landcover%value(c,r,:)).ne.0.and.&
               ifContainsSelectedSurfaceType(n,&
               nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,:)))) then 
             write(LDT_logunit,*) 'The distribution of surface types in '
             write(LDT_logunit,*) 'the LANDCOVER field does not sum to 100%'
             write(LDT_logunit,*) 'mask ',LDT_LSMparam_struc(n)%landmask%value(c,r,1)
             write(LDT_logunit,*) 'veg ',LDT_LSMparam_struc(n)%landcover%value(c,r,:)
             call LDT_endrun()
          endif
       enddo
    enddo
    write(LDT_logunit,*) "Total number of grid points:",LDT_rc%ngrid(n)

    allocate(LDT_domain(n)%tile(LDT_rc%ntiles(n)))
    do m=1,LDT_rc%max_model_types
       allocate(LDT_surface(n,m)%tile(LDT_rc%npatch(n,m)))
    enddo
    allocate(LDT_domain(n)%grid(LDT_rc%ngrid(n)))
    allocate(LDT_domain(n)%gindex(LDT_rc%lnc(n), LDT_rc%lnr(n)))
    
    LDT_domain(n)%minLat = 200.0
    LDT_domain(n)%minLon = 200.00
    LDT_domain(n)%maxLat = -200.0
    LDT_domain(n)%maxLon = -200.00
    
    kk = 1
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
          LDT_domain(n)%gindex(c,r) = -1

          call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
               locallat,locallon)
          
          if(localLat.lt.LDT_domain(n)%minLat) LDT_domain(n)%minLat = localLat
          if(localLat.gt.LDT_domain(n)%maxLat) LDT_domain(n)%maxLat = localLat
          if(localLon.lt.LDT_domain(n)%minLon) LDT_domain(n)%minLon = localLon
          if(localLon.gt.LDT_domain(n)%maxLon) LDT_domain(n)%maxLon = localLon

          if(LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0.and.&
               abs(sum(LDT_LSMparam_struc(n)%landcover%value(c,r,:))-1.0).lt.0.001.and.&
               ifContainsSelectedSurfaceType(n,&
               nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,:)))) then 
             LDT_domain(n)%grid(kk)%lat = locallat
             LDT_domain(n)%grid(kk)%lon = locallon
             LDT_domain(n)%grid(kk)%col = c
             LDT_domain(n)%grid(kk)%row = r
             LDT_domain(n)%gindex(c,r) = kk
             kk = kk+1
          endif
       enddo
    enddo

    kk = 0
    kk_sf = 0 
    do r=1,LDT_rc%lnr(n)
       do c=1,LDT_rc%lnc(n)
          do m=1,LDT_rc%nensem(n)
             if(LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0) then
                
                ntiles_surface = compute_ntiles_surface(n,c,r)
                ntiles_soilt = compute_ntiles_soilt(n,c,r,soilt_selected)
                ntiles_soilf = compute_ntiles_soilf(n,c,r,soilf_selected)
                ntiles_elev = compute_ntiles_elev(n,c,r,elev_selected)
                ntiles_slope = compute_ntiles_slope(n,c,r,slope_selected)
                ntiles_aspect = compute_ntiles_aspect(n,c,r,aspect_selected)

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
                            
                               kk = kk+1

                               LDT_domain(n)%tile(kk)%ensem = m
                               LDT_domain(n)%tile(kk)%row=r    
                               LDT_domain(n)%tile(kk)%col=c    
                               LDT_domain(n)%tile(kk)%index = &
                                    LDT_domain(n)%gindex(c,r)
                               LDT_domain(n)%tile(kk)%vegt=vegt
                               LDT_domain(n)%tile(kk)%soilt=soilt
                               LDT_domain(n)%tile(kk)%sand=sand
                               LDT_domain(n)%tile(kk)%clay=clay
                               LDT_domain(n)%tile(kk)%silt=silt
                               LDT_domain(n)%tile(kk)%sftype = &
                                    sf_index
                               LDT_domain(n)%tile(kk)%elev = &
                                    elev
                               LDT_domain(n)%tile(kk)%slope = &
                                    slope
                               LDT_domain(n)%tile(kk)%aspect = &
                                    aspect
                               LDT_domain(n)%tile(kk)%fgrd = &
                                    LDT_LSMparam_struc(n)%landcover%value(c,r,vegt)
                               
                               if(soilt_selected) then 
                                  LDT_domain(n)%tile(kk)%fgrd = & 
                                       LDT_domain(n)%tile(kk)%fgrd * &
                                       LDT_LSMparam_struc(n)%texture%value(c,r,soilt)
                               elseif(soilf_selected) then 
                                  LDT_domain(n)%tile(kk)%fgrd = & 
                                       LDT_domain(n)%tile(kk)%fgrd * &
                                       LDT_LSMparam_struc(n)%soilsfgrd%value(c,r,soilf_index)
                               endif
                               
                               if(elev_selected) then 
                                  LDT_domain(n)%tile(kk)%fgrd=& 
                                       LDT_domain(n)%tile(kk)%fgrd*&
                                       LDT_LSMparam_struc(n)%elevfgrd%value(c,r,elev_index)
                               endif

                               if(slope_selected) then 
                                  LDT_domain(n)%tile(kk)%fgrd=& 
                                       LDT_domain(n)%tile(kk)%fgrd*&
                                       LDT_LSMparam_struc(n)%slopefgrd%value(c,r,slope_index)
                               endif
                            
                               if(aspect_selected) then 
                                  LDT_domain(n)%tile(kk)%fgrd=& 
                                       LDT_domain(n)%tile(kk)%fgrd*&
                                       LDT_LSMparam_struc(n)%aspectfgrd%value(c,r,aspect_index)
                               endif
                               
                               if(sf_index.gt.0) then 
                                  kk_sf(sf_index) = kk_sf(sf_index) + 1
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%ensem = m
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%tile_id = kk
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%row=r    
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%col=c    
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%index = &
                                       LDT_domain(n)%gindex(c,r)
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%vegt=vegt
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%soilt=soilt
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%sand=sand
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%clay=clay
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%silt=silt

                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%elev=elev
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%slope=slope
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%aspect=aspect
                                  LDT_surface(n,sf_index)%tile(&
                                       kk_sf(sf_index))%fgrd=&
                                       LDT_LSMparam_struc(n)%landcover%value(c,r,vegt)
                                  if(soilt_selected) then 
                                     LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LDT_LSMparam_struc(n)%texture%value(c,r,soilt)
                                  elseif(soilf_selected) then 
                                     LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LDT_LSMparam_struc(n)%soilsfgrd%value(c,r,soilf_index)
                                  endif

                                  if(elev_selected) then 
                                     LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LDT_LSMparam_struc(n)%elevfgrd%value(c,r,elev_index)
                                  endif
                                  if(slope_selected) then 
                                     LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LDT_LSMparam_struc(n)%slopefgrd%value(c,r,slope_index)
                                  endif
                                  if(aspect_selected) then 
                                     LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LDT_surface(n,sf_index)%tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LDT_LSMparam_struc(n)%aspectfgrd%value(c,r,aspect_index)
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
    allocate(sumv(LDT_rc%ngrid(n)))
    sumv = 0 
    do t=1,LDT_rc%ntiles(n)
       c = LDT_domain(n)%tile(t)%col
       r = LDT_domain(n)%tile(t)%row
       gid = LDT_domain(n)%gindex(c,r)
       
       sumv(gid) = sumv(gid) + LDT_domain(n)%tile(t)%fgrd
    enddo
    
    do t=1,LDT_rc%ngrid(n)
       print*, t, sumv(gid)
    enddo
    stop
#endif
    call LDT_domain_setup(n)

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
! This function determines if the input surface type is among 
! the surface model types that are currently chosen in the 
! LDT run
!
!EOP
    logical       :: isSurfaceTypeSelected
    integer       :: i 

    isSurfaceTypeSelected = .false. 
    
    do i=1,LDT_rc%nsf_model_types
       if(LDT_rc%sf_model_type_select(i).eq.surface_type) then 
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
  function ifContainsSelectedSurfaceType(n,surface_types)
! !ARGUMENTS:     
    integer       :: n 
    integer       :: surface_types(LDT_LSMparam_struc(n)%sfctype%num_bins)
! 
! !DESCRIPTION: 
!  This function determines if any of the input list of surface types
!  is among the surface model types that are currently chosen in the 
!  LDT run
!EOP
    logical       :: ifContainsSelectedSurfaceType
    integer       :: i,k

    ifContainsSelectedSurfaceType = .false. 
    
    do k=1,LDT_LSMparam_struc(n)%sfctype%num_bins       
       do i=1,LDT_rc%nsf_model_types
          if(LDT_rc%sf_model_type_select(i).eq.surface_types(k)) then 
             ifContainsSelectedSurfaceType = .true. 
             exit
          endif
       enddo
    enddo

  end function ifContainsSelectedSurfaceType

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
    do t=1,LDT_LSMparam_struc(n)%sfctype%num_bins        
       if(LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0.and.&
            LDT_LSMparam_struc(n)%landcover%value(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(&
            nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,t)))) then
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
    if(LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0.and.&
         LDT_LSMparam_struc(n)%landcover%value(c,r,t).gt.0.0.and.&
         isSurfaceTypeSelected(&
         nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,t)))) then
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
       do t=1,LDT_LSMparam_struc(n)%texture%num_bins       
          if(LDT_LSMparam_struc(n)%texture%value(c,r,t).gt.0.0) then 
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
       do t=1,LDT_LSMparam_struc(n)%soilsfgrd%num_bins
          if(LDT_LSMparam_struc(n)%soilsfgrd%value(c,r,t).gt.0.0) then 
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
       do t=1,LDT_LSMparam_struc(n)%elevfgrd%num_bins       
          if(LDT_LSMparam_struc(n)%elevfgrd%value(c,r,t).gt.0.0) then 
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
       do t=1,LDT_LSMparam_struc(n)%slopefgrd%num_bins        
          if(LDT_LSMparam_struc(n)%slopefgrd%value(c,r,t).gt.0.0) then 
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
       do t=1,LDT_LSMparam_struc(n)%aspectfgrd%num_bins
          if(LDT_LSMparam_struc(n)%aspectfgrd%value(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_aspect = max(1,kk)

  end function compute_ntiles_aspect

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
    do t=1,LDT_LSMparam_struc(n)%sfctype%num_bins
       if(LDT_LSMparam_struc(n)%landcover%value(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(&
            nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,t)))) then 
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
    do t=1,LDT_LSMparam_struc(n)%sfctype%num_bins
       if(LDT_LSMparam_struc(n)%landcover%value(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(&
            nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,t)))) then 
          kk = kk + 1
          if(kk.eq.i) then 
             sf_index = nint(LDT_LSMparam_struc(n)%sfctype%value(c,r,t))
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
       do t=1, LDT_LSMparam_struc(n)%texture%num_bins
          if(LDT_LSMparam_struc(n)%texture%value(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                soilt = t
                return
             endif
          end if
       enddo
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
       do t=1,LDT_LSMparam_struc(n)%soilsfgrd%num_bins
          if(LDT_LSMparam_struc(n)%soilsfgrd%value(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                sand = LDT_LSMparam_struc(n)%sand%value(c,r,t)
                clay = LDT_LSMparam_struc(n)%clay%value(c,r,t)
                silt = LDT_LSMparam_struc(n)%silt%value(c,r,t)
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
       do t=1,LDT_LSMparam_struc(n)%elevfgrd%num_bins
          if(LDT_LSMparam_struc(n)%elevfgrd%value(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                elev = LDT_LSMparam_struc(n)%elevation%value(c,r,t)
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
       do t=1,LDT_LSMparam_struc(n)%slopefgrd%num_bins
          if(LDT_LSMparam_struc(n)%slopefgrd%value(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                slope = LDT_LSMparam_struc(n)%slope%value(c,r,t)
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
       do t=1,LDT_LSMparam_struc(n)%aspectfgrd%num_bins
          if(LDT_LSMparam_struc(n)%aspectfgrd%value(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                aspect = LDT_LSMparam_struc(n)%aspect%value(c,r,t)
                aspect_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_aspect_value

end module LDT_domainMod

