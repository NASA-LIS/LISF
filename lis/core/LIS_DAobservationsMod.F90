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
module LIS_DAobservationsMod

!BOP
!
! !MODULE: LIS_DAobservationsMod
!
! !DESCRIPTION:
!  The code in this file controls the handling of 
!  observations to be used for data assimilation. 
!
!   \item[LIS\_odeltas]
!    size of observation space buffers
!   \item[LIS\_ooffsets]
!    offsets of observation space buffers
!   \item[LIS\_ensOnGrid]
!    grid object for the ensemble space on grid
!   
! !REVISION HISTORY:
!
!  21 Jun 2006: Sujay Kumar; Initial implementation
!
  use ESMF
  use LIS_coreMod
  use LIS_domainMod, only : decompose_nx_ny, decompose_npes
  use LIS_logMod
  use LIS_mpiMod
  use map_utils
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public ::  LIS_initDAobservations
  public ::  LIS_readDAobservations
  public ::  LIS_perturb_DAobservations
  public ::  LIS_convertPatchSpaceToObsSpace
  public ::  LIS_convertPatchSpaceToObsEnsSpace
  public ::  LIS_gather_1dgrid_to_2dgrid_obs
  public ::  LIS_scatter_global_to_local_grid_obs
  public ::  LIS_getObsDomainResolutions
  public ::  LIS_writevar_gridded_obs
  public ::  LIS_convertObsVarToLocalSpace
  public ::  LIS_writevar_innov     ! writes innovations to a gridded file
  public ::  LIS_checkForValidObs
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_OBS_State
  public :: LIS_OBS_Pert_State

  public :: LIS_ews_obs_ind
  public :: LIS_ewe_obs_ind
  public :: LIS_nss_obs_ind
  public :: LIS_nse_obs_ind

  public :: LIS_ews_obs_halo_ind
  public :: LIS_ewe_obs_halo_ind
  public :: LIS_nss_obs_halo_ind
  public :: LIS_nse_obs_halo_ind

  public :: LIS_odeltas
  public :: LIS_ooffsets
  public :: LIS_obsVecGrid
  public :: LIS_obsEnsOnGrid
  public :: LIS_obs_domain
  public :: LIS_obsens_gdeltas
  public :: LIS_obsens_goffsets
  public :: LIS_obs_gdeltas
  public :: LIS_obs_goffsets

  type(ESMF_State), allocatable  :: LIS_OBS_State(:,:) 
  type(ESMF_State), allocatable  :: LIS_OBS_Pert_State(:,:)
  type(ESMF_Grid),  allocatable  :: LIS_obsEnsOnGrid(:,:)
  type(ESMF_Grid),  allocatable  :: LIS_obsVecGrid(:,:)
  integer, allocatable           :: LIS_obs_ngrids(:,:)

  integer, allocatable           :: LIS_obs_deltas(:,:)
  integer, allocatable           :: LIS_obs_offsets(:,:)
  integer, allocatable           :: LIS_obs_gdeltas(:,:)
  integer, allocatable           :: LIS_obs_goffsets(:,:)
  integer, allocatable           :: LIS_obsens_gdeltas(:,:)
  integer, allocatable           :: LIS_obsens_goffsets(:,:)

  integer, allocatable           :: LIS_ews_obs_ind(:,:)
  integer, allocatable           :: LIS_ewe_obs_ind(:,:)
  integer, allocatable           :: LIS_nss_obs_ind(:,:)
  integer, allocatable           :: LIS_nse_obs_ind(:,:)
  integer, allocatable           :: LIS_ews_obs_halo_ind(:,:)
  integer, allocatable           :: LIS_ewe_obs_halo_ind(:,:)
  integer, allocatable           :: LIS_nss_obs_halo_ind(:,:)
  integer, allocatable           :: LIS_nse_obs_halo_ind(:,:)
  integer, allocatable           :: LIS_odeltas(:,:),LIS_ooffsets(:,:)

  type, public :: lis_obs_domain_type
     real                       :: datares
     real,          allocatable :: lat(:)
     real,          allocatable :: lon(:)
     integer,       allocatable :: col(:)
     integer,       allocatable :: row(:)
     integer,       allocatable :: glb_col(:,:)
     integer,       allocatable :: glb_row(:,:)
     type(proj_info)            :: lisproj
     integer, allocatable       :: gindex(:,:)
     real,    allocatable       :: landmask(:,:)
     real                       :: minLat, maxLat, minLon, maxLon
     real                       :: stlat, stlon, truelat1
     real                       :: truelat2, truelon, orient
     real                       :: dx,dy,nlatcircles
     character*20               :: gridtype
     integer, allocatable       :: nbr_index(:)
     real,    allocatable       :: weight(:)
     real,    allocatable       :: rlat(:)
     real,    allocatable       :: rlon(:)
     integer                    :: max_obsngrid
  end type lis_obs_domain_type

  type(lis_obs_domain_type), allocatable :: LIS_obs_domain(:,:)
!EOP


contains

!BOP
! 
! !ROUTINE: LIS_initDAobservations
! \label{LIS_initDAobservations}
! 
! !INTERFACE: 
  subroutine LIS_initDAobservations
! 
! !DESCRIPTION: 
!  This routine allocates the structures required for managing 
!  observational data used for data assimilation. The routine
!  creates ESMF State objects used to incorporate 
!  observational data, and creates the domain definition and 
!  parallel decomposition on the observation grid. 
! 
!  The methods invoked are: 
!  \begin{description}
!  \item[create\_obsdomain\_objects](\ref{create_obsdomain_objects})
!   creates the objects and data structures required for the
!   definition of the observation grid and parallel domain 
!   decomposition. 
!  \item[readObsDataConfig](\ref{readobsdataconfig}) \newline
!    invokes the generic method in the registry to 
!    set up the structures to read the observation data, for 
!    assimilation.
!  \item[perturbinit](\ref{perturbinit}) \newline
!    invokes the generic method in the registry to initialize
!    the observation perturbation algorithm.
!  \item[perturbsetup](\ref{perturbsetup}) \newline
!    invokes the generic method in the registry to populate
!    the observation perturbation objects. 
!  \end{description}
!
! !USES:     
!EOP
    integer                  :: n 
    integer                  :: status
    character*1              :: nestid(2)
    character*1              :: caseid(3)
    character*100            :: temp
    integer                  :: max_index
    integer                  :: i, k
    logical                  :: name_found
    character*20             :: alglist(10)

    TRACE_ENTER("daobs_init")
    if(LIS_rc%ndas.gt.0) then 

       allocate(LIS_OBS_State(LIS_rc%nnest, LIS_rc%ndas))
       allocate(LIS_OBS_Pert_State(LIS_rc%nnest, LIS_rc%ndas))
       allocate(LIS_obsEnsOnGrid(LIS_rc%nnest,LIS_rc%ndas))
       allocate(LIS_obsVecGrid(LIS_rc%nnest,LIS_rc%ndas))

       call create_obsdomain_objects()

       do n=1,LIS_rc%nnest
          do k=1,LIS_rc%ndas
!-----------------------------------------------------------------------------
!    setup interpolation/upscaling weights
!    if LIS is at a finer resolution than the observations, then 
!    LIS fields will be upscaled to the observation space. 
!-----------------------------------------------------------------------------

             if(LIS_isatAfinerResolution(n,LIS_obs_domain(n,k)%datares)) then 
                
                allocate(LIS_obs_domain(n,k)%nbr_index(&
                     LIS_rc%lnc(n)*LIS_rc%lnr(n)))
                allocate(LIS_obs_domain(n,k)%weight(&
                     LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             
                call upscaleByWeightedAveraging_input(&
                     LIS_rc%gridDesc(n,:),&
                     LIS_rc%obs_gridDesc(k,:),&
                     LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                     LIS_obs_domain(n,k)%nbr_index,&
                     LIS_obs_domain(n,k)%weight)
             else
             
                allocate(LIS_obs_domain(n,k)%nbr_index(&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(LIS_obs_domain(n,k)%rlat(&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(LIS_obs_domain(n,k)%rlon(&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
             
                call neighbor_interp_input_withgrid(&
                     LIS_rc%gridDesc(n,:), &
                     LIS_rc%obs_gridDesc(k,:),&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                     LIS_obs_domain(n,k)%rlat, &
                     LIS_obs_domain(n,k)%rlon, &
                     LIS_obs_domain(n,k)%nbr_index)
                
             endif

             write(unit=temp,fmt='(i2.2)') n
             read(unit=temp,fmt='(2a1)') nestid
             
             write(unit=temp,fmt='(i3.3)') k
             read(unit=temp,fmt='(3a1)') caseid
             
             LIS_OBS_State(n,k) = ESMF_StateCreate(name=&
                  "Observations State"//nestid(1)//nestid(2) &
                  //"_"//caseid(1)//caseid(2)//caseid(3), rc=status)
             call LIS_verify(status, 'obs state create')
             
             
             LIS_OBS_Pert_State(n,k) = ESMF_StateCreate(&
                  name="Observations Perturbations"//nestid(1)//nestid(2)//&
                  "_"//caseid(1)//caseid(2)//caseid(3), rc=status)
             call LIS_verify(status, 'obs pert state create')
             
          enddo
       enddo
       write(LIS_logunit,*) '[INFO] Successfully created observations state'

       do k=1, LIS_rc%ndas
          call readObsDataConfig(trim(LIS_rc%daset(k))//char(0),k,&
               LIS_OBS_State(:,k), LIS_OBS_Pert_State(:,k))
       enddo

       max_index = -1
       do i=1,LIS_rc%nperts
          if(max_index.eq.-1.and.LIS_rc%perturb_obs(i).ne."none") then 
             max_index = 1
             alglist(max_index) = LIS_rc%perturb_obs(i)
          else
             name_found = .false. 
             do k=1,max_index
                if(LIS_rc%perturb_obs(i).ne."none".and.&
                     LIS_rc%perturb_obs(i).eq.alglist(k)) then
                   name_found = .true. 
                endif
             enddo
             if(.not.name_found.and.max_index.ne.-1) then 
                max_index = max_index + 1
                alglist(max_index) = LIS_rc%perturb_obs(i)
             endif
          endif
       enddo
       if(max_index.gt.0) then 
          do i=1,max_index
             call perturbinit(trim(alglist(i))//char(0), 3)
          enddo
       endif
       do i=1,LIS_rc%nperts
          if(LIS_rc%perturb_obs(i).ne."none") then 
             call perturbsetup(trim(LIS_rc%perturb_obs(i))//char(0), 3, i,&
                  LIS_OBS_State(:,i), LIS_OBS_Pert_State(:,i))
          endif
       enddo
    endif
    TRACE_EXIT("daobs_init")
  end subroutine LIS_initDAobservations

!BOP
! !ROUTINE: create_obsdomain_objects
! \label{create_obsdomain_objects}
!
! !INTERFACE: 
  subroutine create_obsdomain_objects()
! 
! !DESCRIPTION: 
!   This routine creates the objects required for the definition 
!   of the DA observation space grid. 
!  
!  The methods invoked are: 
!  \begin{description}
!  \item[readObsDomainInput](\ref{readObsDomainInput})
!    This routine reads the run time configuration options
!    related to the observation space.     
!  \item[make\_obs\_domain](\ref{make_obs_domain) \newline
!    This routine creates the required mappings between 
!    the observation space and the model space. 
!  \end{description}
!BOP
    integer              :: odeltas(LIS_rc%ndas)
    integer, allocatable :: deblklist(:,:,:)
    type(ESMF_DistGrid)  :: obsgridDG(LIS_rc%ndas), gridEnsDG(LIS_rc%ndas)
    integer              :: stid, enid
    integer              :: n, i, k    
    integer              :: status, ierr

    call readObsDomainInput()
    call make_obs_domain()
!              print*, 'reached here2 ',LIS_localPet
!          call mpi_barrier(LIS_mpi_comm, ierr)
!          call LIS_endrun


    do n=1,LIS_rc%nnest

       do k=1,LIS_rc%ndas
!          odeltas(k) = LIS_rc%obs_ngrid(k)*LIS_rc%nobtypes(k)
          odeltas(k) = LIS_rc%obs_ngrid(k)
       enddo
!-----------------------------------------------------------------------------
!  The grid, tile space sizes of the decomposed domains are gathered 
!  from each processor to compute the total sizes of the entire domain. 
!-----------------------------------------------------------------------------
#if (defined SPMD)
       do k=1,LIS_rc%ndas
          call MPI_ALLGATHER(odeltas(k),1,MPI_INTEGER,&
               LIS_odeltas(k,:),1,MPI_INTEGER,&
               LIS_mpi_comm,ierr)
       enddo
#else
       do k=1,LIS_rc%ndas
          LIS_odeltas(k,:) = odeltas(k)
       enddo
#endif
       if(LIS_masterproc) then 
          LIS_ooffsets(:,0) = 0 
          do i=1,LIS_npes-1
             do k=1,LIS_rc%ndas
                LIS_ooffsets(k,i) = LIS_ooffsets(k,i-1)+LIS_odeltas(k,i-1)
             enddo
          enddo
       endif

#if (defined SPMD)
       do k=1,LIS_rc%ndas
          call MPI_BCAST(LIS_ooffsets(k,:), LIS_npes, MPI_INTEGER,0, &
               LIS_mpi_comm, ierr)
       enddo
#endif
       do k=1,LIS_rc%ndas

          allocate(deblklist(1,2,LIS_npes))
          do i=0,LIS_npes-1
             stid = LIS_ooffsets(k,i)+1
             enid = stid + LIS_obs_ngrids(k,i)-1
             
             deblklist(:,1,i+1) = (/stid/)
             deblklist(:,2,i+1) = (/enid/)
          enddo
          obsgridDG(k) = ESMF_DistGridCreate(minIndex=(/1/), &
               maxIndex=(/LIS_rc%obs_glbngrid(k)/),&
               deBlockList=deblklist,rc=status)
          call LIS_verify(status, 'ESMF_DistGridCreate failed: obsgridDG')
          deallocate(deblklist)
       
          allocate(deblklist(2,2,LIS_npes))
          do i=0,LIS_npes-1
             stid = LIS_ooffsets(k,i)+1
             enid = stid+LIS_obs_ngrids(k,i)-1
             deblklist(:,1,i+1) = (/stid,1/)
             deblklist(:,2,i+1) = (/enid,LIS_rc%nensem(n)/)
          enddo

          gridEnsDG(k) = ESMF_DistGridCreate(minIndex=(/1,1/),&
               maxIndex=(/LIS_rc%obs_glbngrid(k),LIS_rc%nensem(n)/),&
               deBlockList=deblklist,rc=status)
          call LIS_verify(status, 'ESMF_DistGridCreate failed: gridEnsDG')
          
          LIS_obsVecGrid(n,k) = ESMF_GridCreate(name = "LIS Obs Grid Space",&
               coordTypeKind=ESMF_TYPEKIND_R4, distGrid = obsgridDG(k),&
               gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
          call LIS_verify(status,'ESMF_GridCreate failed: LIS_obsVecGrid')


          LIS_obsEnsOnGrid(n,k) = ESMF_GridCreate(&
               name = "LIS Obs Grid Ensemble Space",&
               coordTypeKind = ESMF_TYPEKIND_R4, distGrid = gridEnsDG(k),&
               gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,0/),rc=status)
          call LIS_verify(status,'ESMF_GridCreate failed: LIS_obsEnsOnGrid')
          
          deallocate(deblklist)
       enddo

    enddo

  end subroutine create_obsdomain_objects

!BOP
! 
! !ROUTINE: readObsDomainInput
! \label{readObsDomainInput}
! 
! !INTERFACE: 
  subroutine readObsDomainInput
! !USES: 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
! 
! !DESCRIPTION: 
!   This routine reads the configurable options for the definition 
!   of the DA observation space. The routine also computes
!   the parallel domain decomposition on the observation 
!   grid. 
! 
!BOP

    integer                    :: n,k
    integer                    :: ios
    integer                    :: ncId, nrId
    integer                    :: ncbId, nrbId
    integer                    :: latId,lonId
    integer                    :: latbId,lonbId
    integer                    :: gnc,gnr
    integer                    :: ftn
    integer                    :: gindex
    logical                    :: file_exists
    character*50               :: map_proj
    integer                    :: c,r,kk
    real,      allocatable     :: lat(:,:)
    real,      allocatable     :: lon(:,:)
    real,      allocatable     :: locallat(:,:)
    real,      allocatable     :: locallon(:,:)
    integer                    :: status

    allocate(LIS_obs_domain(LIS_rc%nnest,LIS_rc%ndas))
    allocate(LIS_rc%obs_gridDesc(LIS_rc%ndas, 50))
    
    allocate(LIS_odeltas(LIS_rc%ndas, 0:LIS_npes-1))
    allocate(LIS_ooffsets(LIS_rc%ndas, 0:LIS_npes-1))
    
    LIS_odeltas = 0
    LIS_ooffsets = 0

    LIS_rc%obs_gridDesc = 0
!    
!  All nests will run the same configuration of DA
!
    allocate(LIS_rc%lis_obs_map_proj(LIS_rc%ndas))
    allocate(LIS_rc%obs_gnc(LIS_rc%ndas))
    allocate(LIS_rc%obs_gnr(LIS_rc%ndas))
    allocate(LIS_rc%obs_lnc(LIS_rc%ndas))
    allocate(LIS_rc%obs_lnr(LIS_rc%ndas))
    allocate(LIS_rc%obs_lnc_red(LIS_rc%ndas))
    allocate(LIS_rc%obs_lnr_red(LIS_rc%ndas))    
    allocate(LIS_ews_obs_ind(LIS_rc%ndas,LIS_npes))
    allocate(LIS_ewe_obs_ind(LIS_rc%ndas,LIS_npes))
    allocate(LIS_nss_obs_ind(LIS_rc%ndas,LIS_npes))
    allocate(LIS_nse_obs_ind(LIS_rc%ndas,LIS_npes))

    allocate(LIS_ews_obs_halo_ind(LIS_rc%ndas,LIS_npes))
    allocate(LIS_ewe_obs_halo_ind(LIS_rc%ndas,LIS_npes))
    allocate(LIS_nss_obs_halo_ind(LIS_rc%ndas,LIS_npes))
    allocate(LIS_nse_obs_halo_ind(LIS_rc%ndas,LIS_npes))

    allocate(LIS_rc%obs_ngrid(LIS_rc%ndas))
    allocate(LIS_rc%obs_glbngrid(LIS_rc%ndas))
    allocate(LIS_rc%obs_glbngrid_red(LIS_rc%ndas))
    allocate(LIS_rc%obsdomainFile(LIS_rc%ndas))

    allocate(LIS_obs_deltas(LIS_rc%ndas,0:LIS_npes-1))
    allocate(LIS_obs_offsets(LIS_rc%ndas,0:LIS_npes-1))

    allocate(LIS_rc%obs_halox(LIS_rc%ndas))
    allocate(LIS_rc%obs_haloy(LIS_rc%ndas))

!--------------------------------------------------------------------------
! When observation is coarser than the model domain, a single halo
! width should be sufficient (except when spatial localization is 
! employed. For now, the observation halo widths are set to 1 obs
! grid cell. 
!--------------------------------------------------------------------------

    LIS_rc%obs_halox = 0
    LIS_rc%obs_haloy = 0 

    write(LIS_logunit,*)'[INFO] DA obs domain details:'
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "Data assimilation observation domain file:",&
         rc=status)

    do k=1,LIS_rc%ndas
       call ESMF_ConfigGetAttribute(LIS_config,&
            LIS_rc%obsdomainFile(k),rc=status)
       call LIS_verify(status,&
            "Data assimilation observation domain file: not defined")
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%ndas
          inquire(file=LIS_rc%obsdomainFile(k), exist=file_exists)
          if(file_exists) then 
             
             ios = nf90_open(path=LIS_rc%obsdomainFile(k),&
                  mode=NF90_NOWRITE,ncid=ftn)
             call LIS_verify(ios,'Error in nf90_open in readObsDomainInput')
             
             ios = nf90_inq_dimid(ftn,"east_west",ncId)
             call LIS_verify(ios,&
                  'Error in nf90_inq_dimid in readObsDomainInput:east_west')
             
             ios = nf90_inq_dimid(ftn,"north_south",nrId)
             call LIS_verify(ios,&
                  'Error in nf90_inq_dimid in readObsDomainInput:north_south')
             
             ios = nf90_inquire_dimension(ftn,ncId, len=gnc)
             call LIS_verify(ios,&
                  'Error in nf90_inquire_dimension in readObsDomainInput:ncId')
             
             ios = nf90_inquire_dimension(ftn,nrId, len=gnr)
             call LIS_verify(ios,&
                  'Error in nf90_inquire_dimension in readObsDomainInput:nrId')
             LIS_rc%obs_gnc(k) = gnc
             LIS_rc%obs_gnr(k) = gnr
             
             allocate(lat(gnc,gnr))
             allocate(lon(gnc,gnr))
             
             ios = nf90_inq_varid(ftn,'lat',latid)
             call LIS_verify(ios,&
                  'lat field not found in the DA obs domain file')
             
             ios = nf90_inq_varid(ftn,'lon',lonid)
             call LIS_verify(ios,&
                  'lon field not found in the DA obs domain file')
             
             ios = nf90_get_var(ftn,latid,lat)
             call LIS_verify(ios,&
                  'Error in nf90_get_var for latid in readObsDomainInput')
             
             ios = nf90_get_var(ftn,lonid,lon)
             call LIS_verify(ios,&
                  'Error in nf90_get_var for lonid in readObsDomainInput')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',map_proj)
             call LIS_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
                  LIS_obs_domain(n,k)%stlat)
             call LIS_verify(ios, &
                  'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
                  LIS_obs_domain(n,k)%stlon)
             call LIS_verify(ios, &
                  'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',LIS_obs_domain(n,k)%dx)
             call LIS_warning(ios, 'Error in nf90_get_att: DX')
             if(ios.ne.0) then 
                LIS_obs_domain(n,k)%dx = 0.0
             endif

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',LIS_obs_domain(n,k)%dy)
             call LIS_warning(ios, 'Error in nf90_get_att: DY')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'GRIDTYPE',&
                  LIS_obs_domain(n,k)%gridtype)
             call LIS_warning(ios, 'Error in nf90_get_att: GRIDTYPE')
             
             if(ios.ne.0) then 
                LIS_obs_domain(n,k)%gridtype = "none"
             endif

             call LIS_quilt_obs_domain(n, k, gnc, gnr)

             allocate(locallat(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
             allocate(locallon(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
             
             locallat = lat(&
                  LIS_ews_obs_halo_ind(k,LIS_localPet+1):&         
                  LIS_ewe_obs_halo_ind(k,LIS_localPet+1), &
                  LIS_nss_obs_halo_ind(k,LIS_localPet+1): &
                  LIS_nse_obs_halo_ind(k,LIS_localPet+1))
          
             locallon = lon(&
                  LIS_ews_obs_halo_ind(k,LIS_localPet+1):&         
                  LIS_ewe_obs_halo_ind(k,LIS_localPet+1), &
                  LIS_nss_obs_halo_ind(k,LIS_localPet+1): &
                  LIS_nse_obs_halo_ind(k,LIS_localPet+1))
             
             allocate(LIS_obs_domain(n,k)%lat(&
                  LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
             allocate(LIS_obs_domain(n,k)%lon(&
                  LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
             
             do r=1,LIS_rc%obs_lnr(k)
                do c=1,LIS_rc%obs_lnc(k)
                   gindex = c+(r-1)*LIS_rc%obs_lnc(k)
                   LIS_obs_domain(n,k)%lat(gindex) = locallat(c,r)
                   LIS_obs_domain(n,k)%lon(gindex) = locallon(c,r)
                enddo
             enddo                                     
             
             LIS_rc%obs_gridDesc(k,32) = gnc
             LIS_rc%obs_gridDesc(k,33) = gnr
             LIS_rc%obs_gridDesc(k,34) = lat(1,1)
             LIS_rc%obs_gridDesc(k,35) = lon(1,1)
             LIS_rc%obs_gridDesc(k,37) = lat(gnc,gnr)
             LIS_rc%obs_gridDesc(k,38) = lon(gnc,gnr)
             LIS_rc%obs_gridDesc(k,39) = LIS_obs_domain(n,k)%dx
             LIS_rc%obs_gridDesc(k,40) = LIS_obs_domain(n,k)%dy
             
             if(map_proj.eq."EQUIDISTANT CYLINDRICAL") then              
                LIS_rc%lis_obs_map_proj(k) = "latlon"
                
                LIS_rc%obs_gridDesc(k,4) = LIS_obs_domain(n,k)%stlat + &
                     (LIS_nss_obs_halo_ind(k,LIS_localPet+1)-1)*&
                     LIS_obs_domain(n,k)%dy
                LIS_rc%obs_gridDesc(k,5) = LIS_obs_domain(n,k)%stlon + &
                     (LIS_ews_obs_halo_ind(k,LIS_localPet+1)-1)*&
                     LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,7) = LIS_obs_domain(n,k)%stlat + & 
                     (LIS_nse_obs_halo_ind(k,LIS_localPet+1)-1)*&
                     LIS_obs_domain(n,k)%dy
                LIS_rc%obs_gridDesc(k,8) = LIS_obs_domain(n,k)%stlon + &
                     (LIS_ewe_obs_halo_ind(k,LIS_localPet+1)-1)*&
                     LIS_obs_domain(n,k)%dx
                
                write(LIS_logunit,*) '[INFO] local obs domain ',&
                     LIS_rc%obs_gridDesc(k,4),LIS_rc%obs_gridDesc(k,7),&
                  LIS_rc%obs_gridDesc(k,5),LIS_rc%obs_gridDesc(k,8), 'd: ',n
                
                LIS_rc%obs_gridDesc(k,1) = 0
                LIS_rc%obs_gridDesc(k,9) = LIS_obs_domain(n,k)%dx
                
                if(LIS_rc%obs_gridDesc(k,1).eq.0) then 
                   LIS_rc%obs_gridDesc(k,10) = LIS_obs_domain(n,k)%dy
                   LIS_rc%obs_gridDesc(k,6) = 128
                   LIS_rc%obs_gridDesc(k,11) = 64
                   LIS_rc%obs_gridDesc(k,20) = 64
                endif
                if(LIS_rc%obs_gridDesc(k,7).lt.LIS_rc%obs_gridDesc(k,4)) then
                   write(LIS_logunit,*) '[ERR] lat2 must be greater than lat1'
                   write(LIS_logunit,*) '[ERR] ',LIS_rc%obs_gridDesc(k,7),&
                        LIS_rc%obs_gridDesc(k,4)
                   write(LIS_logunit,*) '[ERR] Stopping run...'
                   call LIS_endrun
                endif
                if(LIS_rc%obs_gridDesc(k,8).lt.LIS_rc%obs_gridDesc(k,5)) then
                   write(LIS_logunit,*) '[ERR] lon2 must be greater than lon1'
                   write(LIS_logunit,*) '[ERR] ',LIS_rc%obs_gridDesc(k,8),&
                        LIS_rc%obs_gridDesc(k,5)
                   write(LIS_logunit,*) '[ERR] Stopping run...'
                   call LIS_endrun
                endif
                
                LIS_rc%obs_gridDesc(k,2) = nint((LIS_rc%obs_gridDesc(k,8)-&
                     LIS_rc%obs_gridDesc(k,5))&
                     /LIS_rc%obs_gridDesc(k,9))+ 1
                LIS_rc%obs_gridDesc(k,3) = nint((LIS_rc%obs_gridDesc(k,7)-&
                     LIS_rc%obs_gridDesc(k,4))&
                     /LIS_rc%obs_gridDesc(k,10)) + 1  
                
                
                !local domain of each processor
                call map_set(PROJ_LATLON, LIS_rc%obs_gridDesc(k,4),&
                     LIS_rc%obs_gridDesc(k,5),&
                     0.0, LIS_rc%obs_gridDesc(k,9),&
                     LIS_rc%obs_gridDesc(k,10), 0.0,&
                     LIS_rc%obs_lnc(k),&
                     LIS_rc%obs_lnr(k),LIS_obs_domain(n,k)%lisproj)
                
             elseif(map_proj.eq."LAMBERT CONFORMAL") then      

                LIS_rc%lis_obs_map_proj(k) = "lambert"
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',&
                     LIS_obs_domain(n,k)%truelat1)
                call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT1')
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT2',&
                     LIS_obs_domain(n,k)%truelat2)
                call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT2')
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',&
                     LIS_obs_domain(n,k)%truelon)
                call LIS_verify(ios, 'Error in nf90_get_att: STANDARD_LON')
                
                LIS_rc%obs_gridDesc(k,1) = 3
                LIS_rc%obs_gridDesc(k,2) = LIS_rc%obs_lnc(k)
                LIS_rc%obs_gridDesc(k,3) = LIS_rc%obs_lnr(k)
                LIS_rc%obs_gridDesc(k,4) = &
                     lat(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,5) = &
                     lon(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,6) = 8
                LIS_rc%obs_gridDesc(k,7) = LIS_obs_domain(n,k)%truelat2
                LIS_rc%obs_gridDesc(k,8) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,9) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,10) = LIS_obs_domain(n,k)%truelat1
                LIS_rc%obs_gridDesc(k,11) = LIS_obs_domain(n,k)%truelon
                
                call map_set(PROJ_LC,LIS_rc%obs_gridDesc(k,4),&
                     LIS_rc%obs_gridDesc(k,5),&
                     LIS_rc%obs_gridDesc(k,8)*1000.0,&
                     LIS_rc%obs_gridDesc(k,11),&
                     LIS_rc%obs_gridDesc(k,10),&
                     LIS_domain(n)%truelat2,&
                     LIS_rc%lnc(n),LIS_rc%lnr(n),&
                     LIS_obs_domain(n,k)%lisproj)
                
             elseif(map_proj.eq."MERCATOR") then              

                LIS_rc%lis_obs_map_proj(k) = "mercator"
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',&
                     LIS_obs_domain(n,k)%truelat1)
                call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT1')
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',&
                     LIS_obs_domain(n,k)%truelon)
                call LIS_verify(ios, 'Error in nf90_get_att: STANDARD_LON')
                
                LIS_rc%obs_gridDesc(k,1) = 1
                LIS_rc%obs_gridDesc(k,2) = LIS_rc%obs_lnc(k)
                LIS_rc%obs_gridDesc(k,3) = LIS_rc%obs_lnr(k)
                LIS_rc%obs_gridDesc(k,4) = &
                     lat(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,5) = &
                     lon(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,6) = 8
                LIS_rc%obs_gridDesc(k,7) = 0.0
                LIS_rc%obs_gridDesc(k,8) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,9) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,10) = LIS_obs_domain(n,k)%truelat1
                LIS_rc%obs_gridDesc(k,11) = LIS_obs_domain(n,k)%truelon
                
                call map_set(PROJ_MERC,LIS_rc%obs_gridDesc(k,4),&
                     LIS_rc%obs_gridDesc(k,5),&
                     LIS_rc%obs_gridDesc(k,8)*1000.0,&
                     LIS_rc%obs_gridDesc(k,11),&
                     LIS_rc%obs_gridDesc(k,10),&
                     0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),&
                     LIS_obs_domain(n,k)%lisproj)
             
             elseif(map_proj.eq."POLAR STEREOGRAPHIC") then              

                LIS_rc%lis_obs_map_proj(k) = "polar"
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',&
                     LIS_obs_domain(n,k)%truelat1)
                call LIS_verify(ios, 'Error in nf90_get_att: TRUELAT1')
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'ORIENT',&
                     LIS_obs_domain(n,k)%orient)
                call LIS_verify(ios, 'Error in nf90_get_att: ORIENT')
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',&
                     LIS_obs_domain(n,k)%truelon)
                call LIS_verify(ios, 'Error in nf90_get_att: STANDARD_LON')
                
                LIS_rc%obs_gridDesc(k,1) = 5
                LIS_rc%obs_gridDesc(k,2) = LIS_rc%obs_lnc(k)
                LIS_rc%obs_gridDesc(k,3) = LIS_rc%obs_lnr(k)
                LIS_rc%obs_gridDesc(k,4) = &
                     lat(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,5) = &
                     lon(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,6) = 8
                LIS_rc%obs_gridDesc(k,7) = LIS_obs_domain(n,k)%orient
                LIS_rc%obs_gridDesc(k,8) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,9) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,10) = LIS_obs_domain(n,k)%truelat1
                LIS_rc%obs_gridDesc(k,11) = LIS_obs_domain(n,k)%truelon
                
                call map_set(PROJ_PS,LIS_rc%obs_gridDesc(k,4),&
                     LIS_rc%obs_gridDesc(k,5),&
                     LIS_rc%obs_gridDesc(k,8)*1000.0,&
                     LIS_rc%obs_gridDesc(k,11),LIS_rc%obs_gridDesc(k,10),&
                     0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),&
                     LIS_obs_domain(n,k)%lisproj)
                
             elseif(map_proj.eq."EASE V2") then 

                LIS_rc%lis_obs_map_proj(k) = "ease_v2"

                LIS_rc%obs_gridDesc(k,1) = 9
                LIS_rc%obs_gridDesc(k,2) = LIS_rc%obs_lnc(k)
                LIS_rc%obs_gridDesc(k,3) = LIS_rc%obs_lnr(k)
                LIS_rc%obs_gridDesc(k,4) = &
                     lat(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,5) = &
                     lon(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,6) = 128
                LIS_rc%obs_gridDesc(k,7) = &
                     lat(LIS_ewe_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nse_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,8) = &
                     lon(LIS_ewe_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nse_obs_halo_ind(k,LIS_localPet+1))
                if(LIS_obs_domain(n,k)%gridtype.eq."M36") then 
                   LIS_rc%obs_gridDesc(k,9) = 4
                   LIS_rc%obs_gridDesc(k,10) = 0.36
                   LIS_obs_domain(n,k)%dx = 0.36
                   LIS_obs_domain(n,k)%dy = 0.36
                elseif(LIS_obs_domain(n,k)%gridtype.eq."M09") then 
                   LIS_rc%obs_gridDesc(k,9) = 5
                   LIS_rc%obs_gridDesc(k,10) = 0.09
                   LIS_obs_domain(n,k)%dx = 0.09
                   LIS_obs_domain(n,k)%dy = 0.09
                endif

                !local domain of each processor
                call map_set(PROJ_EASEV2,LIS_rc%obs_gridDesc(k,4),&
                     LIS_rc%obs_gridDesc(k,5),&
                     LIS_rc%obs_gridDesc(k,9),&
                     LIS_rc%obs_gridDesc(k,9),&
                     LIS_rc%obs_gridDesc(k,9),0.0,&
                     LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k),&
                     LIS_obs_domain(n,k)%lisproj)

             elseif(map_proj.eq."GAUSSIAN") then 

                LIS_rc%lis_obs_map_proj(k) = "gaussian"
                
                ios = nf90_get_att(ftn, NF90_GLOBAL, 'NUMBER OF LAT CIRCLES',&
                     LIS_obs_domain(n,k)%nlatcircles)
                call LIS_verify(ios, &
                     'Error in nf90_get_att: NUMBER_OF_LAT_CIRCLES')
                
                LIS_rc%obs_gridDesc(k,1) = 4
                LIS_rc%obs_gridDesc(k,2) = LIS_rc%obs_lnc(k)
                LIS_rc%obs_gridDesc(k,3) = LIS_rc%obs_lnr(k)
                LIS_rc%obs_gridDesc(k,4) = &
                     lat(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,5) = &
                     lon(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,6) = 128
                LIS_rc%obs_gridDesc(k,7) = &
                     lat(LIS_ewe_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nse_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,8) = &
                     lon(LIS_ewe_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nse_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,9) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,10) = LIS_obs_domain(n,k)%nlatcircles
                LIS_rc%obs_gridDesc(k,11) = 64
                LIS_rc%obs_gridDesc(k,20) = 64
                
             elseif(map_proj.eq."HRAP") then 
                
                LIS_rc%lis_obs_map_proj(k) = "hrap"

                LIS_rc%obs_gridDesc(k,1) = 8
                LIS_rc%obs_gridDesc(k,2) = LIS_rc%obs_lnc(k)
                LIS_rc%obs_gridDesc(k,3) = LIS_rc%obs_lnr(k)
                LIS_rc%obs_gridDesc(k,4) = &
                     lat(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,5) = &
                     lon(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,6) = 8
                LIS_rc%obs_gridDesc(k,7) = LIS_obs_domain(n,k)%orient
                LIS_rc%obs_gridDesc(k,8) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,9) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,10) = LIS_obs_domain(n,k)%truelat1
                LIS_rc%obs_gridDesc(k,11) = LIS_obs_domain(n,k)%truelon
             
                call map_set(PROJ_HRAP,LIS_rc%obs_gridDesc(k,4),&
                     LIS_rc%obs_gridDesc(k,5),&
                     LIS_rc%obs_gridDesc(k,8)*1000.0,&
                     LIS_rc%obs_gridDesc(k,11),&
                     LIS_rc%obs_gridDesc(k,10),&
                     0.0,LIS_rc%lnc(n),LIS_rc%lnr(n),&
                     LIS_obs_domain(n,k)%lisproj)
                
             elseif(map_proj.eq."UTM") then 
                
                LIS_rc%lis_obs_map_proj(k) = "utm"

                LIS_rc%obs_gridDesc(k,1) = 7
                LIS_rc%obs_gridDesc(k,2) = LIS_rc%obs_lnc(k)
                LIS_rc%obs_gridDesc(k,3) = LIS_rc%obs_lnr(k)
                LIS_rc%obs_gridDesc(k,4) = &
                     lat(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,5) = &
                     lon(LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nss_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,6) = 128
                LIS_rc%obs_gridDesc(k,7) = &
                     lat(LIS_ewe_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nse_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,8) = &
                     lon(LIS_ewe_obs_halo_ind(k,LIS_localPet+1),&
                     LIS_nse_obs_halo_ind(k,LIS_localPet+1))
                LIS_rc%obs_gridDesc(k,9) = LIS_obs_domain(n,k)%dx
                LIS_rc%obs_gridDesc(k,10) = LIS_obs_domain(n,k)%dy
                LIS_rc%obs_gridDesc(k,11) = 64
                LIS_rc%obs_gridDesc(k,20) = 64
             endif
             LIS_obs_domain(n,k)%datares = max(LIS_obs_domain(n,k)%dx, &
                  LIS_obs_domain(n,k)%dy)
             

             deallocate(lat)
             deallocate(lon)
             deallocate(locallat)
             deallocate(locallon)
             
          else
             write(LIS_logunit,*) '[ERR] ',&
                  LIS_rc%obsdomainFile(k), ' does not exist'
             write(LIS_logunit,*) '[ERR] program stopping ...'
             call LIS_endrun
          endif

          write(LIS_logunit,*)&
               '[INFO] --------------------Obs domain ',n,'----------------------'
          do kk=1,13
             write(LIS_logunit,*) &
                  '[INFO] (',kk,',',LIS_rc%obs_gridDesc(k,kk),')'
          enddo
          write(LIS_logunit,*)&
               '[INFO] --------------------------------------------------------------'
          
       enddo
    enddo
#endif

  end subroutine readObsDomainInput

!BOP
! !ROUTINE: make_obs_domain
! \label{make\_obs\_domain} 
!
! !INTERFACE: 
  subroutine make_obs_domain()
! !USES: 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
! 
! !DESCRIPTION: 
!  This routine relies on the LDT-generated observation domain 
!  file to create the objects for parallelization of the 
!  DA observation space. The objects required for mapping 
!  between the observation space and the model space are
!  also created in this routine. 
!  
!EOP
    integer              :: n,k,t,l
    logical              :: file_exists
    integer              :: nid, ncId, nrId, maskid
    integer              :: i,c,r,c1,r1,c2,r2,kk,nc,nr
    integer              :: ios, ierr
    integer              :: count,gid,gid_red
    integer              :: ntiles, stid
    integer              :: ngrids
    integer              :: gdeltas
    real, allocatable    :: gmask(:,:)
    integer, allocatable :: gtmp(:,:)
    integer, allocatable :: gtmp1(:)

    allocate(LIS_obs_ngrids(LIS_rc%ndas,0:LIS_npes-1))

    allocate(LIS_obs_gdeltas(LIS_rc%ndas,0:LIS_npes-1))
    allocate(LIS_obs_goffsets(LIS_rc%ndas,0:LIS_npes-1))
    allocate(LIS_obsens_gdeltas(LIS_rc%ndas,0:LIS_npes-1))
    allocate(LIS_obsens_goffsets(LIS_rc%ndas,0:LIS_npes-1))

    !read the obs landmask

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%ndas

          inquire(file=trim(LIS_rc%obsDomainFile(k)), exist=file_exists)
          if(file_exists) then 
             
             allocate(LIS_obs_domain(n,k)%landmask(&
                  LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
             allocate(LIS_obs_domain(n,k)%gindex(&
                  LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
             allocate(gmask(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))

             write(LIS_logunit,*)'[INFO] Reading obs landmask map from ',&
                  trim(LIS_rc%obsDomainFile(k))
             ios = nf90_open(path=trim(LIS_rc%obsDomainFile(k)),&
                  mode=NF90_NOWRITE,ncid=nid)
             call LIS_verify(ios,'Error in nf90_open in make_obs_domain')
             
             ios = nf90_inq_dimid(nid,"east_west",ncId)
             call LIS_verify(ios,'Error in nf90_inq_dimid in make_obs_domain')
             
             ios = nf90_inq_dimid(nid,"north_south",nrId)
             call LIS_verify(ios,'Error in nf90_inq_dimid in make_obs_domain')
             
             ios = nf90_inquire_dimension(nid,ncId, len=nc)
             call LIS_verify(ios,&
                  'Error in nf90_inquire_dimension in make_obs_domain')
             
             ios = nf90_inquire_dimension(nid,nrId, len=nr)
             call LIS_verify(ios,&
                  'Error in nf90_inquire_dimension in make_obs_domain')

             ios = nf90_inq_varid(nid,"LANDMASK",maskid)
             call LIS_verify(ios,'LANDMASK field not found in the obs domain file')

             ios = nf90_get_var(nid,maskid,LIS_obs_domain(n,k)%landmask)
             call LIS_verify(ios,'Error in nf90_get_var in make_obs_domain')
             
             ios = nf90_close(nid)
             call LIS_verify(ios,'Error in nf90_close in make_obs_domain')
             
             gmask(:,:) = &
                  LIS_obs_domain(n,k)%landmask(&
                  LIS_ews_obs_halo_ind(k,LIS_localPet+1):&         
                  LIS_ewe_obs_halo_ind(k,LIS_localPet+1), &
                  LIS_nss_obs_halo_ind(k,LIS_localPet+1): &
                  LIS_nse_obs_halo_ind(k,LIS_localPet+1))

             write(LIS_logunit,*)'[INFO] Successfully read obs landmask map '
             
          else
             write(LIS_logunit,*) '[ERR] obs domain landmask map: ',&
                  trim(LIS_rc%obsDomainFile(k)), ' does not exist'
             write(LIS_logunit,*) '[ERR] program stopping ...'
             call LIS_endrun
          endif
       
          LIS_rc%obs_ngrid(k)=0
          LIS_obs_domain(n,k)%gindex = -1
          kk = 1

          do r=1,LIS_rc%obs_lnr(k)
             do c=1,LIS_rc%obs_lnc(k)
                if(gmask(c,r).gt.0) then 
                   LIS_obs_domain(n,k)%gindex(c,r) = kk
                   LIS_rc%obs_ngrid(k) = LIS_rc%obs_ngrid(k)+1             
                   kk = kk + 1
                endif
             enddo
          enddo
          deallocate(gmask)
!-----------------------------------------------------------------------------
!  The grid space sizes of the decomposed domains are gathered 
!  from each processor to compute the total sizes of the entire domain. 
!-----------------------------------------------------------------------------
          ngrids = LIS_rc%obs_ngrid(k)
          gdeltas = LIS_rc%obs_ngrid(k)
          
#if (defined SPMD)
          call MPI_ALLGATHER(ngrids,1,MPI_INTEGER,&
               LIS_obs_ngrids(k,:),1,MPI_INTEGER,&
               LIS_mpi_comm,ierr)
          call MPI_ALLGATHER(gdeltas,1,MPI_INTEGER,&
               LIS_obs_gdeltas(k,:),1,MPI_INTEGER,&
               LIS_mpi_comm,ierr)
!          call MPI_ALLGATHER(LIS_rc%obs_lnc(k),1,MPI_INTEGER,&
!               LIS_obs_lnc(k,:),1,MPI_INTEGER,&
!               LIS_mpi_comm,ierr)
!          call MPI_ALLGATHER(LIS_rc%obs_lnr(k),1,MPI_INTEGER,&
!               LIS_obs_lnr(k,:),1,MPI_INTEGER,&
!               LIS_mpi_comm,ierr)
#else
          LIS_obs_ngrids(k,:) = ngrids
          LIS_obs_gdeltas(k,:) = gdeltas
!          LIS_obs_lnc(k,:) = LIS_rc%obs_lnc(k)
!          LIS_obs_lnr(k,:) = LIS_rc%obs_lnr(k)
          
#endif

!-----------------------------------------------------------------------------
! Create the global mapping
!-----------------------------------------------------------------------------
          LIS_rc%obs_glbngrid_red(k) = 0
          do r=1,LIS_rc%obs_gnr(k)
             do c=1,LIS_rc%obs_gnc(k)
                if(LIS_obs_domain(n,k)%landmask(c,r).gt.0) then 
                   LIS_rc%obs_glbngrid_red(k) =&
                        LIS_rc%obs_glbngrid_red(k) + 1
                   LIS_obs_domain(n,k)%landmask(c,r) = &
                        LIS_rc%obs_glbngrid_red(k)   
                endif
             enddo
          enddo

          if(LIS_masterproc) then 
             LIS_obs_domain(n,k)%max_obsngrid = 0 
             do l=1,LIS_npes
                LIS_obs_domain(n,k)%max_obsngrid = &
                     max(LIS_obs_domain(n,k)%max_obsngrid, &
                     LIS_obs_ngrids(k,l-1))
             enddo
             allocate(LIS_obs_domain(n,k)%glb_col(LIS_npes, &
                  LIS_obs_domain(n,k)%max_obsngrid))
             allocate(LIS_obs_domain(n,k)%glb_row(LIS_npes, &
                  LIS_obs_domain(n,k)%max_obsngrid))
             
             LIS_obs_domain(n,k)%glb_col = -1
             LIS_obs_domain(n,k)%glb_row = -1

!  Use landmask to store the global 2d to 1d mapping 
             do l=1,LIS_npes                
                count = 1
                do r=LIS_nss_obs_halo_ind(k,l),LIS_nse_obs_halo_ind(k,l)
                   do c=LIS_ews_obs_halo_ind(k,l),LIS_ewe_obs_halo_ind(k,l)
                      if(LIS_obs_domain(n,k)%landmask(c,r).gt.0) then 
                         LIS_obs_domain(n,k)%glb_col(l,count) = c
                         LIS_obs_domain(n,k)%glb_row(l,count) = r
                         count = count + 1
                      endif
                   enddo
                enddo
                      
             enddo

          endif

          if(LIS_masterproc) then 
             LIS_obs_goffsets(k,0) = 0 
             do i=1,LIS_npes-1
                LIS_obs_goffsets(k,i) = LIS_obs_goffsets(k,i-1)+&
                     LIS_obs_gdeltas(k,i-1)
             enddo
          end if

#if (defined SPMD)
          call MPI_BCAST(LIS_obs_goffsets(k,:), LIS_npes, MPI_INTEGER,0, &
               LIS_mpi_comm, ierr)
#endif
          gdeltas = LIS_rc%obs_ngrid(k)*LIS_rc%nensem(n)
          
#if (defined SPMD)
          call MPI_ALLGATHER(gdeltas,1,MPI_INTEGER,&
               LIS_obsens_gdeltas(k,:),1,MPI_INTEGER,&
               LIS_mpi_comm,ierr)

#else
          LIS_obs_gdeltas(k,:) = gdeltas
          
#endif

          if(LIS_masterproc) then 
             LIS_obsens_goffsets(k,0) = 0 
             do i=1,LIS_npes-1
                LIS_obsens_goffsets(k,i) = LIS_obsens_goffsets(k,i-1)+&
                     LIS_obsens_gdeltas(k,i-1)
             enddo
          end if

#if (defined SPMD)
          call MPI_BCAST(LIS_obsens_goffsets(k,:), LIS_npes, MPI_INTEGER,0, &
               LIS_mpi_comm, ierr)
#endif
          
          
          LIS_rc%obs_glbngrid(k) = 0 
          do i=0,LIS_npes-1
             LIS_rc%obs_glbngrid(k) = &
                  LIS_rc%obs_glbngrid(k) + LIS_obs_ngrids(k,i)
          enddo
          
          
          allocate(LIS_obs_domain(n,k)%col(LIS_rc%obs_ngrid(k)))
          allocate(LIS_obs_domain(n,k)%row(LIS_rc%obs_ngrid(k)))
          
          LIS_obs_domain(n,k)%col = -1
          LIS_obs_domain(n,k)%row = -1

          do r=1,LIS_rc%obs_lnr(k)
             do c=1,LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                   LIS_obs_domain(n,k)%col(LIS_obs_domain(n,k)%gindex(c,r)) = c
                   LIS_obs_domain(n,k)%row(LIS_obs_domain(n,k)%gindex(c,r)) = r
                endif
             enddo
          enddo

          
       enddo
    enddo
#endif

  end subroutine make_obs_domain
!BOP
! !ROUTINE: LIS_quilt_obs_domain
! \label{LIS_quilt_obs_domain}
! 
! !INTERFACE:
  subroutine LIS_quilt_obs_domain(n,k, nc, nr)
! !USES:     

! !ARGUMENTS: 
    integer,          intent(in)  :: n
    integer,          intent(in)  :: k
    integer,          intent(in)  :: nc
    integer,          intent(in)  :: nr
! 
! !DESCRIPTION: 
! This routine generates the quilted domain extents and sizes based on the 
! processor layout and the size of the specified halos for the 
! DA observation domain 
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

   LIS_rc%obs_lnc(k) = ewe_halo_ind - ews_halo_ind + 1
   LIS_rc%obs_lnr(k) = nse_halo_ind - nss_halo_ind + 1

   LIS_rc%obs_lnc_red(k)= ewe_ind - ews_ind + 1
   LIS_rc%obs_lnr_red(k)= nse_ind - nss_ind + 1

   write(unit=LIS_logunit,fmt=*)'[INFO] local obs domain',':(', &
                                LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k),')'
   write(unit=LIS_logunit,fmt=*)'[INFO] local obs domain without halo',':(', &
                                LIS_rc%obs_lnc_red(k),LIS_rc%obs_lnr_red(k),')'
   write(unit=LIS_logunit,fmt=*)'[INFO] running obs domain',':(',nc,nr,')'

   deltas = LIS_rc%obs_lnc_red(k)*LIS_rc%obs_lnr_red(k)

#if (defined SPMD)
   call MPI_ALLGATHER(deltas,1,MPI_INTEGER,LIS_obs_deltas(k,:),1,MPI_INTEGER,   &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(ews_ind,1,MPI_INTEGER,LIS_ews_obs_ind(k,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nss_ind,1,MPI_INTEGER,LIS_nss_obs_ind(k,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(ewe_ind,1,MPI_INTEGER,LIS_ewe_obs_ind(k,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nse_ind,1,MPI_INTEGER,LIS_nse_obs_ind(k,:),1,MPI_INTEGER, &
                      LIS_mpi_comm,ierr)

   call MPI_ALLGATHER(ews_halo_ind,1,MPI_INTEGER,LIS_ews_obs_halo_ind(k,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nss_halo_ind,1,MPI_INTEGER,LIS_nss_obs_halo_ind(k,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(ewe_halo_ind,1,MPI_INTEGER,LIS_ewe_obs_halo_ind(k,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
   call MPI_ALLGATHER(nse_halo_ind,1,MPI_INTEGER,LIS_nse_obs_halo_ind(k,:),1, &
                      MPI_INTEGER,LIS_mpi_comm,ierr)
#else
   LIS_obs_deltas(k,:)  = deltas
   LIS_ews_obs_ind(k,:) = ews_ind
   LIS_nss_obs_ind(k,:) = nss_ind
   LIS_ewe_obs_ind(k,:) = ewe_ind
   LIS_nse_obs_ind(k,:) = nse_ind

   LIS_ews_obs_halo_ind(k,:) = ews_halo_ind
   LIS_nss_obs_halo_ind(k,:) = nss_halo_ind
   LIS_ewe_obs_halo_ind(k,:) = ewe_halo_ind
   LIS_nse_obs_halo_ind(k,:) = nse_halo_ind
#endif   
   if(LIS_masterproc) then
      LIS_obs_offsets(k,0) = 0
      do i=1,LIS_npes-1
         LIS_obs_offsets(k,i) = LIS_obs_offsets(k,i-1)+LIS_obs_deltas(k,i-1)
      enddo
   end if
#if (defined SPMD)
   call MPI_BCAST(LIS_obs_offsets(k,:), LIS_npes, MPI_INTEGER,0, LIS_mpi_comm, ierr)
#endif

   end subroutine LIS_quilt_obs_domain

!BOP
! 
! !ROUTINE: LIS_readDAobservations
! \label{LIS_readDAobservations}
! 
! !INTERFACE: 
  subroutine LIS_readDAobservations(n)
! !USES: 


! !ARGUMENTS: 
    integer, intent(in) :: n

! 
! !DESCRIPTION: 
!  This routines calls the method to read the specific observation
!  source. 
! 
!  The methods invoked are: 
!  \begin{description}
!  \item[readDAObs](\ref{readdaobs}) \newline
!   invokes the method in the registry to read the observations and 
!   incorporate into the observation state.
!  \end{description}
!EOP
    integer          :: k 

    TRACE_ENTER("daobs_read")
    do k=1,LIS_rc%ndas
       call readDAObs(trim(LIS_rc%daset(k))//char(0), n, k, &
            LIS_OBS_State(n,k), LIS_OBS_Pert_State(n,k))
    enddo
    TRACE_EXIT("daobs_read")

  end subroutine LIS_readDAobservations

!BOP
! 
! !ROUTINE: LIS_perturb_DAobservations
! \label{LIS_perturb_DAobservations}
!
! !INTERFACE: 
  subroutine LIS_perturb_DAobservations(n)
! !USES: 

! !ARGUMENTS: 
    integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routines calls the method to perturb the observations based on the
!  specified choice. 
! 
!  The methods invoked are: 
!  \begin{description}
!  \item[perturbmethod](\ref{perturbmethod}) \newline
!   invokes the specified algorithm to perturb the observations. 
!  \item[writeDAObs](\ref{writedaobs}) \newline
!   invokes the registry method to write the perturbed and 
!   processed observations to disk. 
!  \end{description}
!EOP
    real                    :: curr_time
    logical                 :: data_status
    integer                 :: status
    integer                 :: k
    integer                 :: i,j,m
    integer                 :: Nobjs
    real                    :: delta
    character*100,    allocatable :: obs_state_objs(:)
    type(ESMF_Field), allocatable :: obs_field(:)
    real,             pointer     :: obs_temp(:,:)

    TRACE_ENTER("daobs_perturb")
    do k=1,LIS_rc%ndas
       if(LIS_rc%perturb_obs(k).ne."none") then 
          curr_time = float(LIS_rc%hr)*3600+60*float(LIS_rc%mn)+float(LIS_rc%ss)
!          if(mod(curr_time,real(LIS_rc%pertobsInterval(k))).eq.0)then
          call ESMF_AttributeGet(LIS_OBS_State(n,k),"Data Update Status",&
               data_status,rc=status)
          call LIS_verify(status, &
               "ESMF_AttributeGet: Data Update Status failed in LIS_perturb_DAobservations")
          if(data_status) then
             call perturbmethod(trim(LIS_rc%perturb_obs(k))//char(0),3,n,k, &
                  LIS_OBS_State(n,k), LIS_OBS_Pert_State(n,k))

             if(LIS_rc%pert_bias_corr.eq.1) then 
                
                call ESMF_StateGet(LIS_OBS_State(n,k),itemCount=Nobjs,rc=status)
                call LIS_verify(status, &
                     'ESMF_StateGet failed in LIS_DAobservationsMod')
                
                allocate(obs_state_objs(Nobjs))
                allocate(obs_field(Nobjs))
                call ESMF_StateGet(LIS_OBS_Pert_State(n,k),&
                     itemNameList=obs_state_objs,rc=status)
                call LIS_verify(status,&
                     'ESMF_StateGet failed in LIS_DAobservationsMod')        
                
                do i=1,Nobjs
                   call ESMF_StateGet(LIS_OBS_Pert_State(n,k),obs_state_objs(i),&
                        obs_field(i),rc=status)
                   call LIS_verify(status,&
                        'ESMF_StateGet failed in LIS_DAobservationsMod')
                   call ESMF_FieldGet(obs_field(i),localDE=0,farrayPtr=obs_temp,rc=status)
                   call LIS_verify(status,&
                        'ESMF_FieldGet failed in LIS_DAobservationsMod')
                   
                   do j=1, LIS_rc%obs_ngrid(k)
                      delta    = 0.0
                      obs_temp(j,LIS_rc%nensem(n)) = 0.0
                      delta = sum(obs_temp(j,1:LIS_rc%nensem(n)-1))/&
                           (LIS_rc%nensem(n)-1)
                      do m=1,LIS_rc%nensem(n)-1
                         obs_temp(j,m) = obs_temp(j,m) - delta
                      enddo
                   enddo
                enddo

                deallocate(obs_state_objs)
                deallocate(obs_field)
             endif

          endif
       endif

       if(LIS_rc%wobs(k).gt.0) then 
          call writeDAObs(trim(LIS_rc%daset(k))//char(0), n, k, &
               LIS_OBS_State(n,k))
       endif

    enddo
    TRACE_EXIT("daobs_perturb")
    
  end subroutine LIS_perturb_DAobservations

!BOP
! !ROUTINE: LIS_convertPatchSpaceToObsSpace
! \label{LIS_convertPatchSpaceToObsSpace}
! 
! !INTERFACE:
  subroutine LIS_convertPatchSpaceToObsSpace(&
       n,&
       k,&
       patch_index, &
       pvar, &
       ovar)
!
! !ARGUMENTS: 
    integer,          intent(in) :: n 
    integer,          intent(in) :: k
    integer,          intent(in) :: patch_index
    real                         :: pvar(LIS_rc%npatch(n,patch_index))
    real                         :: ovar(LIS_rc%obs_ngrid(k))
!
! !DESCRIPTION: 
!
!  This routine converts a variable in the patch space to the observation
!  grid space. If the observation space is at a coarser resolution than 
!  the patch space, then the variable is upscaled. On the other hand, the
!  variable is spatially interpolated to the observation space. 
! 
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [patch\_index]
!     index of the patch to which the variable belong to
!   \item [pvar]
!     variable in the patch space
!   \item [ovar]
!     variable in the observation space
!  \end{description}
!  
!EOP
    integer                      :: c,r,t,g,gid
    real                         :: lis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    integer                      :: nlis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: li(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    real                         :: obs_gvar(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                      :: iret

    TRACE_ENTER("daobs_patch2obs")
    lis_gvar  = 0.0
    nlis_gvar = 0

    do t=1,LIS_rc%npatch(n,patch_index)
       c = LIS_surface(n,patch_index)%tile(t)%col
       r = LIS_surface(n,patch_index)%tile(t)%row
       gid = c+(r-1)*LIS_rc%lnc(n)
       lis_gvar(gid)  = lis_gvar(gid) + pvar(t)
       nlis_gvar(gid) = nlis_gvar(gid) + 1
    enddo

    li = .false. 
    do g=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
       if(nlis_gvar(g).ne.0) then 
          lis_gvar(g)  = lis_gvar(g)/ &
               nlis_gvar(g) 
          li(g) = .true. 
       else
          lis_gvar(g) = LIS_rc%udef
       endif
    enddo

    if(LIS_isatAfinerResolution(n,LIS_obs_domain(n,k)%datares)) then     
       call upscaleByWeightedAveraging(&
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
            LIS_rc%udef, &
            LIS_obs_domain(n,k)%nbr_index, &
            LIS_obs_domain(n,k)%weight, &
            li, &
            lis_gvar, &
            lo, &
            obs_gvar)
    else
       call neighbor_interp(LIS_rc%obs_gridDesc(k,:), &
            li, lis_gvar, lo, obs_gvar,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
            LIS_obs_domain(n,k)%rlat, &
            LIS_obs_domain(n,k)%rlon, &
            LIS_obs_domain(n,k)%nbr_index, &
            LIS_rc%udef,iret)
    endif
    ovar = LIS_rc%udef
    do r=1,LIS_rc%obs_lnr(k)
       do c=1,LIS_rc%obs_lnc(k)          
          if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
             ovar(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                  obs_gvar(c+(r-1)*LIS_rc%obs_lnc(k))
          endif
       enddo
    enddo
    TRACE_EXIT("daobs_patch2obs")

  end subroutine LIS_convertPatchSpaceToObsSpace

!BOP
! !ROUTINE: LIS_convertPatchSpaceToObsEnsSpace
! \label{LIS_convertPatchSpaceToObsEnsSpace}
! 
! !INTERFACE:
  subroutine LIS_convertPatchSpaceToObsEnsSpace(&
       n,&
       k,&
       patch_index, &
       pvar, &
       ovar)
! !ARGUMENTS: 
    integer,          intent(in) :: n 
    integer,          intent(in) :: k
    integer,          intent(in) :: patch_index
    real                         :: pvar(LIS_rc%npatch(n,patch_index))
    real                         :: ovar(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION: 
!
!  This routine converts a variable in the patch space to the observation
!  ensemble grid space. If the observation space is at a coarser resolution than 
!  the patch space, then the variable is upscaled. On the other hand, the
!  variable is spatially interpolated to the observation space. These 
!  transformations are done separately for each ensemble member. 
! 
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [patch\_index]
!     index of the patch to which the variable belong to
!   \item [pvar]
!     variable in the patch space
!   \item [ovar]
!     variable in the observation ensemble space
!  \end{description}
!  
!EOP

    integer                      :: c,r,t,i,m,g,gid
    real                         :: lis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    integer                      :: nlis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: li(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    real                         :: obs_gvar(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                      :: iret

    ovar = LIS_rc%udef
    TRACE_ENTER("daobs_patch2obsEns")
    do m=1,LIS_rc%nensem(n)
       lis_gvar  = 0.0
       nlis_gvar = 0
       obs_gvar = LIS_rc%udef

       do i=1,LIS_rc%npatch(n,patch_index), LIS_rc%nensem(n)
          t = i+m-1
          c = LIS_surface(n,patch_index)%tile(t)%col
          r = LIS_surface(n,patch_index)%tile(t)%row
          gid = c+(r-1)*LIS_rc%lnc(n)
          lis_gvar(gid)  = lis_gvar(gid) + pvar(t)
          nlis_gvar(gid) = nlis_gvar(gid) + 1
       enddo
       
       li = .false. 
       do g=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
          if(nlis_gvar(g).ne.0) then 
             lis_gvar(g)  = lis_gvar(g)/ &
                  nlis_gvar(g) 
             li(g) = .true. 
          else
             lis_gvar(g) = LIS_rc%udef
          endif
       enddo

       if(LIS_isatAfinerResolution(n,LIS_obs_domain(n,k)%datares)) then     
          call upscaleByAveraging(&
               LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
               LIS_rc%udef, &
               LIS_obs_domain(n,k)%nbr_index, &
               li, &
               lis_gvar, &
               lo, &
               obs_gvar)
       else          
          call neighbor_interp(LIS_rc%obs_gridDesc(k,:), &
               li, lis_gvar, lo, obs_gvar,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
               LIS_obs_domain(n,k)%rlat, &
               LIS_obs_domain(n,k)%rlon, &
               LIS_obs_domain(n,k)%nbr_index, &
               LIS_rc%udef,iret)
       endif
       
       do r=1,LIS_rc%obs_lnr(k)
          do c=1,LIS_rc%obs_lnc(k)          
             if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                ovar(LIS_obs_domain(n,k)%gindex(c,r),m) = & 
                     obs_gvar(c+(r-1)*LIS_rc%obs_lnc(k))                
             endif
          enddo
       enddo
    enddo
 
    TRACE_EXIT("daobs_patch2obsEns")
  end subroutine LIS_convertPatchSpaceToObsEnsSpace

!BOP
!
! !ROUTINE: LIS_gather_1dgrid_to_2dgrid_obs
! \label{LIS_gather_1dgrid_to_2dgrid_obs}
!
! !INTERFACE:
subroutine LIS_gather_1dgrid_to_2dgrid_obs(n, k, gtmp, var)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer, intent(in)          :: n
   integer, intent(in)          :: k
   real, allocatable            :: gtmp(:,:)
   real, intent(in)             :: var(LIS_rc%obs_ngrid(k))

!
! !DESCRIPTION:
! This routine gathers a gridded 1d array in the DA observation space
! into a 2d global array (in the DA observation space). 
! The mapping process accounts for the observation halo.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   real, allocatable  :: gtmp1d(:)
   integer :: i,c,r,m,t,l
   integer :: ntiles, gid, count1
   integer :: ierr
   integer :: gdeltas

   TRACE_ENTER("daobs_oneD2twoD")
   if ( LIS_masterproc ) then 
      allocate(gtmp1d(LIS_rc%obs_glbngrid(k)))
      allocate(gtmp(LIS_rc%obs_gnc(k), LIS_rc%obs_gnr(k)))
      gtmp = 0.0     
   else
      allocate(gtmp1d(1))
      allocate(gtmp(1,1))
      gtmp = 0.0
   endif

#if (defined SPMD)     
   gdeltas = LIS_obs_gdeltas(k,LIS_localPet)
   call MPI_GATHERV(var,gdeltas,&
        MPI_REAL,gtmp1d,LIS_obs_gdeltas(k,:),LIS_obs_goffsets(k,:),&
        MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
   gtmp1d= var
#endif

   if ( LIS_masterproc ) then 
      gtmp = LIS_rc%udef
      count1=1
      
      do l=1,LIS_npes
         do t =1, LIS_obs_ngrids(k,l-1)
            c = LIS_obs_domain(n,k)%glb_col(l,t)
            r = LIS_obs_domain(n,k)%glb_row(l,t)
            if(LIS_obs_domain(n,k)%landmask(c,r).gt.0) then 
               gtmp(c,r) = gtmp1d(count1)
               count1 = count1 + 1
            endif
         enddo
      enddo
   endif
   deallocate(gtmp1d)
   TRACE_EXIT("daobs_oneD2twoD")
  
 end subroutine LIS_gather_1dgrid_to_2dgrid_obs


!BOP
!
! !ROUTINE: LIS_scatter_global_to_local_grid_obs
! \label{LIS_scatter_global_to_local_grid_obs}
!
! !INTERFACE:
 subroutine LIS_scatter_global_to_local_grid_obs(n, k, gtmp,ltmp)
! !USES: 

! !ARGUMENTS: 
   implicit none

   integer, intent(in)          :: n
   integer, intent(in)          :: k
   real, intent(in)             :: gtmp(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k))
   real                         :: ltmp(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))

!
! !DESCRIPTION:
! This routine distributes (scatters) a gridded global 2d array in the 
! DA observation space into a loal 2d array (in the DA observation 
! space). This process accounts for the observation halos.
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
!EOP

   real, allocatable :: gtmp1d(:)
   real, allocatable :: gtmp2d(:,:)
   integer :: i,c,r,m,t,l
   integer :: gid
   integer :: ierr

   TRACE_ENTER("daobs_gblTwo2locTwo")
   allocate(gtmp1d(LIS_rc%obs_gnc(k)*LIS_rc%obs_gnr(k)))
   gtmp1d = LIS_rc%udef

   if ( LIS_masterproc ) then 
      do r=1,LIS_rc%obs_gnr(k)
         do c=1,LIS_rc%obs_gnc(k)
            gid = c+(r-1)*LIS_rc%obs_gnc(k)
            gtmp1d(gid) = gtmp(c,r)
         enddo
      enddo
   endif
#if (defined SPMD)      
   call MPI_BCAST(gtmp1d, LIS_rc%obs_gnc(k)*LIS_rc%obs_gnr(k), MPI_REAL,0, &
        LIS_mpi_comm, ierr)
   call LIS_verify(ierr, 'MPI_BCAST failed in LIS_scatter_global_to_local_grid')
#endif
   allocate(gtmp2d(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
   do r=1,LIS_rc%obs_gnr(k)
      do c=1,LIS_rc%obs_gnc(k)
         gid = c+(r-1)*LIS_rc%obs_gnc(k)
         gtmp2d(c,r) = gtmp1d(gid)
      enddo
   enddo
   deallocate(gtmp1d)
   ltmp(:,:) = gtmp2d(LIS_ews_obs_halo_ind(k,LIS_localPet+1):&         
          LIS_ewe_obs_halo_ind(k,LIS_localPet+1), &
          LIS_nss_obs_halo_ind(k,LIS_localPet+1): &
          LIS_nse_obs_halo_ind(k,LIS_localPet+1))
   deallocate(gtmp2d)
   TRACE_EXIT("daobs_gblTwo2locTwo")

 end subroutine LIS_scatter_global_to_local_grid_obs


!BOP
!
! !ROUTINE: LIS_getObsDomainResolutions
! \label{LIS_getObsDomainResolutions}
!
! !INTERFACE:
  subroutine LIS_getObsDomainResolutions(k,dx,dy)

  implicit none
! !ARGUMENTS:
    integer, intent(in) :: k
    real                :: dx, dy
!
! !DESCRIPTION:
!
! This function determines spatial resolution of the DA 
! observation grid (in degrees)
!
! The arguments are:
!  \begin{description}
!  \item[k]
!    index of the DA observation instance
!  \item[dx]
!    spatial resolution along the east-west dimension
!  \item[dy]
!    spatial resolution along the north-south dimension
!  \end{description}
!EOP

    TRACE_ENTER("daobs_getDomRes")
    if ( LIS_rc%obs_gridDesc(k,1) .eq. 0 ) then
       dx = LIS_rc%obs_gridDesc(k,9)
       dy = LIS_rc%obs_gridDesc(k,10)
    elseif ( LIS_rc%obs_gridDesc(k,1) .eq. 1 ) then
       dx = LIS_rc%obs_gridDesc(k,8)/100.0
       dy = LIS_rc%obs_gridDesc(k,9)/100.0
    elseif ( LIS_rc%obs_gridDesc(k,1) .eq. 3 ) then
       dx = LIS_rc%obs_gridDesc(k,8)/100.0
       dy = LIS_rc%obs_gridDesc(k,9)/100.0
    elseif ( LIS_rc%obs_gridDesc(k,1) .eq. 5 ) then
       dx = LIS_rc%obs_gridDesc(k,8)/100.0
       dy = LIS_rc%obs_gridDesc(k,9)/100.0
    elseif ( LIS_rc%obs_gridDesc(k,1) .eq. 9 ) then
       dx = LIS_rc%obs_gridDesc(k,10)
       dy = LIS_rc%obs_gridDesc(k,10)
    else
       
       write(LIS_logunit,*) '[ERR] LIS_getObsDomainResolutions routine ' // &
                            '[ERR] is NOT supported'
       write(LIS_logunit,*) '[ERR] for this map projection.'
       call LIS_endrun()

    endif
    TRACE_EXIT("daobs_getDomRes")
  end subroutine LIS_getObsDomainResolutions

!BOP
! !ROUTINE: LIS_writevar_gridded_obs
! \label{LIS_writevar_gridded_obs}
! 
! !INTERFACE:
  subroutine LIS_writevar_gridded_obs(ftn, n, k, var)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)          :: ftn
    integer, intent(in)          :: n
    integer, intent(in)          :: k
    real, intent(inout)          :: var(LIS_rc%obs_ngrid(k))

! !DESCRIPTION:
!  Writes a real variable in the DA observation space to a binary 
!  sequential access, gridded file as either 2-dimensional gridded field. 
!
!  The arguments are: 
!  \begin{description}
!   \item [ftn]
!     unit number of the binary output file
!   \item [n]
!     index of the domain or nest.
!  \item[k]
!    index of the DA observation instance
!   \item [var]
!     variables being written, dimensioned in the 1-d DA observation grid.
!  \end{description}
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    real, allocatable :: gtmp2(:)
    integer       :: l
    integer       :: c,r,t,count1, ierr
    integer       :: gdeltas
    
    TRACE_ENTER("daobs_writeGridObs")
    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
       allocate(gtmp1(LIS_rc%obs_glbngrid(k)))
       gtmp = 0.0
       gtmp1 = 0.0
    else
       allocate(gtmp1(1))
       gtmp1 = 0.0
    endif
#if (defined SPMD)
    gdeltas = LIS_obs_gdeltas(k,LIS_localPet)
    call MPI_GATHERV(var,gdeltas,MPI_REAL,gtmp1,&
         LIS_obs_gdeltas(k,:),LIS_obs_goffsets(k,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp1 = var
#endif

    if(LIS_masterproc) then 
       gtmp = LIS_rc%udef
       count1=1

       do l=1,LIS_npes
          do t =1, LIS_obs_ngrids(k,l-1)
             c = LIS_obs_domain(n,k)%glb_col(l,t)
             r = LIS_obs_domain(n,k)%glb_row(l,t)
             if(LIS_obs_domain(n,k)%landmask(c,r).gt.0) then 
                if(gtmp(c,r).eq.LIS_rc%udef) then 
                   gtmp(c,r) = gtmp1(count1)
                endif
                count1 = count1 + 1
             endif
          enddo
       enddo
       write(ftn) gtmp
       deallocate(gtmp)
    endif
    deallocate(gtmp1)
    TRACE_EXIT("daobs_writeGridObs")
    
  end subroutine LIS_writevar_gridded_obs


!BOP
! !ROUTINE: LIS_writevar_innov
! \label{LIS_writevar_innov}
!
! !INTERFACE:
  subroutine LIS_writevar_innov(ftn, n, k, varid, var)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: varid
    real, intent(in)    :: var(LIS_rc%obs_ngrid(k)*LIS_rc%nobtypes(k))
!
! !DESCRIPTION:
!  Writes the innovations to a NetCDF file in a gridded format. 
!
!  The arguments are: 
!  \begin{description}
!   \item [ftn]
!     unit number of the binary output file
!   \item [n]
!     index of the domain or nest.
!  \item[k]
!    index of the DA observation instance
!  \item[varid]
!    id of the variable being written
!   \item [var]
!     variables being written, dimensioned in the 1-d observation space
!  \end{description}
!EOP

    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp2(:)
    real, allocatable :: gtmp1(:,:)
    integer :: ierr
    integer :: odeltas
    integer :: i,c,r,l,t,count1

    TRACE_ENTER("daobs_writeInnov")
    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%obs_glbngrid(k)*LIS_rc%nobtypes(k)))
       allocate(gtmp1(LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k)))
    else
       allocate(gtmp(1))
    endif
#if (defined SPMD)
    odeltas = LIS_odeltas(k,LIS_localPet)
    call MPI_GATHERV(var,odeltas,&
         MPI_REAL,gtmp,LIS_odeltas(k,:),LIS_ooffsets(k,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
    gtmp = var
#endif
    if(LIS_masterproc) then 
       do i=1,LIS_rc%nobtypes(k)
          gtmp1 = LIS_rc%udef
          count1=1

          do l=1,LIS_npes
             do t =1, LIS_obs_ngrids(k,l-1)
                c = LIS_obs_domain(n,k)%glb_col(l,t)
                r = LIS_obs_domain(n,k)%glb_row(l,t)
                if(LIS_obs_domain(n,k)%landmask(c,r).gt.0) then 
                   if(gtmp1(c,r).eq.LIS_rc%udef) then 
                      gtmp1(c,r) = gtmp(count1)
                   endif
                   count1 = count1 + 1
                endif
             enddo
          enddo
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
          ierr = nf90_put_var(ftn,varid,gtmp1,(/1,1,1/),&
               (/LIS_rc%obs_gnc(k),LIS_rc%obs_gnr(k),i/))
          call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')
#endif          
          
       enddo
       deallocate(gtmp1)
    endif
    deallocate(gtmp)
    TRACE_EXIT("daobs_writeInnov")
  end subroutine LIS_writevar_innov

!BOP
! !ROUTINE: LIS_convertObsVarToLocalSpace
! \label{LIS_convertObsVarToLocalSpace}
! 
! !INTERFACE:
  subroutine LIS_convertObsVarToLocalSpace(n,k, gvar,lvar)
! !USES:

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: n
    integer, intent(in)   :: k
    real                  :: gvar(LIS_rc%obs_glbngrid_red(k))
    real                  :: lvar(LIS_rc%obs_ngrid(k))
! !DESCRIPTION:
!  Converts a variable in the global 1-d DA observation space
!  to local 1-d DA observation space. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \item[k]
!    index of the DA observation instance
!   \item [gvar]
!    input variable dimensioned in the 1-d global DA observation space
!   \item [lvar]
!    output variable dimensioned in the 1-d local DA observation space
!  \end{description}
!EOP
    integer             :: count1
    integer             :: c,r,t
    integer             :: gid, tid, ntiles, stid
    integer             :: lmask

    TRACE_ENTER("daobs_gblOne2lclOne")
    do t=1, LIS_rc%obs_ngrid(k)
       c = LIS_obs_domain(n,k)%col(t) + LIS_ews_obs_halo_ind(k,LIS_localPet+1)-1
       r = LIS_obs_domain(n,k)%row(t) + LIS_nss_obs_halo_ind(k,LIS_localPet+1)-1

       lmask = int(LIS_obs_domain(n,k)%landmask(c,r))
       if(lmask.gt.0) then 
          lvar(t) = gvar(lmask)
       endif
    enddo

#if 0 
    count1 = 1
    stid = LIS_obs_goffsets(k,LIS_localPet) + 1

    do r=LIS_nss_obs_halo_ind(k,LIS_localPet+1),&
         LIS_nse_obs_halo_ind(k,LIS_localPet+1)
       do c=LIS_ews_obs_halo_ind(k,LIS_localPet+1),&
            LIS_ewe_obs_halo_ind(k,LIS_localPet+1)
          if(LIS_obs_domain(n,k)%landmask(c,r).gt.0) then 
             lvar(count1) = gvar(stid+count1) 
             count1 = count1 + 1
          endif
       enddo
    enddo
#endif   
    TRACE_EXIT("daobs_gblOne2lclOne")

  end subroutine LIS_convertObsVarToLocalSpace

!BOP
!
! !ROUTINE: LIS_checkForValidObs
! \label{LIS_checkForValidObs}
!
! !INTERFACE:
  subroutine  LIS_checkForValidObs(n, k,obs_1d,validObsFlag,obs_2d)
! !ARGUMENTS:     
    integer                 :: n
    integer                 :: k
    real                    :: obs_1d(LIS_rc%obs_ngrid(k))
    integer                 :: validObsFlag
    real                    :: obs_2d(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
! 
! !DESCRIPTION: 
! 
!  This routine checks to see if there are valid observations in the
!  local observation space. The routine also converts the observations
!  in the 1-d space to a 2-d grid. 
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [obs\_1d]
!     observations in the 1-d vector space
!   \item [validObsFlag]
!     flag to indicate if there are valid observations in the local 
!     processor's domain
!   \item [obs\_2d]
!     observations in the 2-d grid space
!  \end{description}
!EOP


    integer                 :: c,r

    TRACE_ENTER("daobs_checkObs")
    validObsFlag = 0 
    obs_2d       = LIS_rc%udef

    do r =1,LIS_rc%obs_lnr(k)
       do c =1,LIS_rc%obs_lnc(k)
          if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
             obs_2d(c,r) = obs_1d(LIS_obs_domain(n,k)%gindex(c,r))
             if(obs_2d(c,r).ne.LIS_rc%udef) then 
                validObsFlag = 1
             endif
          end if
       end do
    end do
    TRACE_EXIT("daobs_checkObs")

  end subroutine LIS_checkForValidObs

end module LIS_DAobservationsMod
