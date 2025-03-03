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
module LDT_coreMod
!BOP
!
! !MODULE: LDT_coreMod
! 
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for the operation of LDT 
!
!  \subsubsection{Overview}
!  This module contains the defintion and specification of the basic
!  datastructures that define a LDT instance. It consists of: 
!  \begin{description}
!   \item[LDT\_rc]
!    datastructure containg generic LDT run control (rc) variables
!   \item[LDT\_domain]
!    datastructure containing grid, tile spaces, 
!    grid projection, domain decomposition information.
!   \item[LDT\_histData]
!     datastructure defining the metadata for
!     model output.
!   \item[LDT\_config] 
!     instance of the configuration class. 
!   \item[LDT\_vm] 
!     instance of the virtual machine controlling
!     the LDT resources. 
!   \item[LDT\_masterproc]
!    logical variable to specify the master processor (true = for master, 
!    false for non-zero processors)
!   \item[LDT\_npes]
!    number of processors used in LDT
!   \item[LDT\_ews\_ind] 
!    starting East West index of each nest (excluding halo regions)
!   \item[LDT\_ewe\_ind] 
!    ending East West index of each nest (excluding halo regions)
!   \item[LDT\_nss\_ind] 
!    starting North-South index of each nest (excluding halo regions)
!   \item[LDT\_nse\_ind] 
!    ending North-South index of each nest (excluding halo regions)
!   \item[LDT\_ews\_ind\_halo] 
!    starting East West index of each nest (including halo regions)
!   \item[LDT\_ewe\_ind\_halo] 
!    ending East West index of each nest (including halo regions)
!   \item[LDT\_nss\_ind\_halo] 
!    starting North-South index of each nest (including halo regions)
!   \item[LDT\_nse\_ind\_halo] 
!    ending North-South index of each nest (including halo regions)
!   \item[LDT\_deltas] 
!    size of unmasked grid space buffers
!   \item[LDT\_offsets] 
!   offsets for unmasked grid space buffers
!   \item[LDT\_tdeltas]
!   size of tile space buffers
!   \item[LDT\_toffsets]  
!    offsets of tile space buffers
!   \item[LDT\_gdeltas]    
!    size of grid space buffers
!   \item[LDT\_goffsets]
!    offsets of grid space buffers
!   \item[LDT\_odeltas]
!    size of observation space buffers
!   \item[LDT\_ooffsets]
!    offsets of observation space buffers
!   \item[LDT\_ntiless]
!    size of tile spaces
!    \item[LDT\_ngrids]
!    size of grid spaces
!   \item[LDT\_vecTile]
!    grid object for the tile space-vector of tiles
!   \item[LDT\_vecGrid]
!    grid object for the grid space-vector of grid pts
!   \item[LDT\_ensOnGrid]
!    grid object for the ensemble space on grid
!  \end{description}
!
! !REVISION HISTORY: 
!  02 Oct 2008   Sujay Kumar   Initial Specification
! 
! !USES:       
  use ESMF
  use LDT_timeMgrMod
  use LDT_PRIV_rcMod
  use LDT_PRIV_gridMod
  use LDT_PRIV_tileMod
  use LDT_logMod
  use LDT_mpiMod
  use map_utils

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_configinit
  public :: LDT_coreInit      
  public :: LDT_ticktime
  public :: LDT_endofrun
  public :: LDT_isLDTatAfinerResolution
  public :: LDT_endofTimeWindow 
!  public :: LDT_TimeToRunNest   ! in LIS and not LDT
!  public :: LDT_resetTimeMgr    ! in LIS and not LDT
  public :: LDT_finalize        
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LDT_rc
  public :: LDT_domain 
  public :: LDT_surface
  public :: LDT_config
  public :: LDT_vm
  public :: LDT_masterproc
  public :: LDT_localPet
  public :: LDT_npes
  public :: LDT_ews_ind   
  public :: LDT_ewe_ind   
  public :: LDT_nss_ind   
  public :: LDT_nse_ind   
  public :: LDT_ews_halo_ind
  public :: LDT_ewe_halo_ind
  public :: LDT_nss_halo_ind
  public :: LDT_nse_halo_ind
  public :: LDT_ews_b_ind
  public :: LDT_ewe_b_ind
  public :: LDT_nss_b_ind
  public :: LDT_nse_b_ind

  public :: LDT_deltas     
  public :: LDT_offsets   
  public :: LDT_tdeltas    
  public :: LDT_toffsets   
  public :: LDT_gdeltas    
  public :: LDT_goffsets   
  public :: LDT_odeltas    
  public :: LDT_ooffsets   
  public :: LDT_patch_deltas
  public :: LDT_patch_offsets
  public :: LDT_ntiless       
  public :: LDT_ngrids     
  public :: LDT_npatches
! Three different types of grid objects are defined because ESMF3.1.0 does
! not support defining tiles on a single grid. 
  public :: LDT_vecTile 
  public :: LDT_vecPatch
  public :: LDT_vecGrid   
  public :: LDT_ensOnGrid 
!EOP

  type, public :: ldt_domain_type 
     type(proj_info)            :: ldtproj
     type(proj_info)            :: ldtglbproj
     integer,       allocatable :: gindex(:,:)
     real,          allocatable :: lat(:)
     real,          allocatable :: lon(:)
     real,          allocatable :: lat_b(:)
     real,          allocatable :: lon_b(:)
     type(griddec), allocatable :: grid(:)
     type(tiledec), allocatable :: tile(:)
     integer, allocatable       :: ntiles_pergrid(:)
     integer, allocatable       :: str_tind(:) !starting tile id for each grid
     integer, allocatable       :: datamask(:,:)
     real                       :: minLat, maxLat, minLon, maxLon
     real                       :: stlat, stlon, truelat1
     real                       :: truelat2, truelon, orient
     real                       :: dx,dy,nlatcircles
  end type ldt_domain_type

  type, public :: ldt_domain_sf_type 
     type(tiledec), allocatable :: tile(:)
     integer,       allocatable :: npatch_pergrid(:)
     integer,       allocatable :: str_patch_ind(:)
  end type ldt_domain_sf_type

  type(ldtrcdec), save   :: LDT_rc
  type(ldt_domain_type), allocatable :: LDT_domain(:)
  type(ldt_domain_sf_type), allocatable :: LDT_surface(:,:)
  type(ESMF_Config),save :: LDT_config
  type(ESMF_VM), save    :: LDT_vm
  type(ESMF_DELayout), save :: LDT_DElayout
  logical                :: LDT_masterproc
  integer                :: LDT_localPet
  integer                :: LDT_npes

  integer, allocatable   :: LDT_ews_ind(:,:), LDT_ewe_ind(:,:)
  integer, allocatable   :: LDT_nss_ind(:,:), LDT_nse_ind(:,:)
  integer, allocatable   :: LDT_ews_halo_ind(:,:), LDT_ewe_halo_ind(:,:)
  integer, allocatable   :: LDT_nss_halo_ind(:,:), LDT_nse_halo_ind(:,:)
  integer, allocatable   :: LDT_ews_b_ind(:,:), LDT_ewe_b_ind(:,:)
  integer, allocatable   :: LDT_nss_b_ind(:,:), LDT_nse_b_ind(:,:)
  integer, allocatable   :: LDT_deltas(:,:), LDT_offsets(:,:)
  integer, allocatable   :: LDT_tdeltas(:,:), LDT_toffsets(:,:)
  integer, allocatable   :: LDT_gdeltas(:,:), LDT_goffsets(:,:)
  integer, allocatable   :: LDT_odeltas(:,:), LDT_ooffsets(:,:)  ! original LDT
!  integer, allocatable   :: LDT_odeltas(:,:,:), LDT_ooffsets(:,:,:)  ! based on LIS
  integer, allocatable   :: LDT_ntiless(:,:),LDT_ngrids(:,:)
  integer, allocatable   :: LDT_npatches(:,:,:)
  integer, allocatable   :: LDT_patch_deltas(:,:,:), LDT_patch_offsets(:,:,:)

  type(ESMF_Grid), allocatable  :: LDT_vecTile(:)
  type(ESMF_Grid), allocatable  :: LDT_vecPatch(:,:)
  type(ESMF_Grid), allocatable  :: LDT_vecGrid(:)
  type(ESMF_Grid), allocatable  :: LDT_ensOnGrid(:)
  
contains

!BOP
! !ROUTINE: LDT_configinit
! \label{LDT_configinit}
!
! !INTERFACE: 
  subroutine LDT_configinit(configfile)
!
! !DESCRIPTION: 
!  Performs the initialization of LDT runtime configuration for an 
!  offline (uncoupled simulation).   
!
!  The Calling sequence is : 
!  \begin{description}
!    \item[spmd\_init\_offline] (\ref{spmd_init_offline}) \newline
!     performs SPMD initializations
!    \item[LDT\_log\_init] (\ref{LDT_log_init}) \newline
!     initializes the LDT log handler
!   \item[LDT\_readConfig] (\ref{LDT_readConfig}) \newline
!     reads the model independent options
!   \item[spmd\_setup] (\ref{spmd_setup}) \newline
!     allocates memory for SPMD variables 
!   \item[LDT\_timemgr\_init] (\ref{LDT_timemgr_init}) \newline
!    initializes the time manager
!   \end{description}
!EOP
    implicit none
    character(len=*) :: configfile

    integer :: status

    call spmd_init_offline(LDT_vm)
    call LDT_log_init(LDT_getNextUnitNumber())
!    call ldt_process_cmd_args
    call LDT_readConfig(configfile)

! - Original LDT code:
    call spmd_setup(LDT_rc%nnest,LDT_rc%max_model_types)
! - LIS code:
!    call spmd_setup(LDT_rc%nnest, LDT_rc%ndas, &
!         LDT_rc%max_model_types)
!    call LDT_timemgr_init(LDT_rc)
! - LIS

    LDT_rc%endtime = 0

    allocate(LDT_domain(LDT_rc%nnest))
    allocate(LDT_surface(LDT_rc%nnest,LDT_rc%max_model_types))
    allocate(LDT_rc%ntiles(LDT_rc%nnest))
    allocate(LDT_rc%ngrid(LDT_rc%nnest))

    allocate(LDT_rc%lnc(LDT_rc%nnest))
    allocate(LDT_rc%lnr(LDT_rc%nnest))
    allocate(LDT_rc%gnc_buf(LDT_rc%nnest))
    allocate(LDT_rc%gnr_buf(LDT_rc%nnest))
    allocate(LDT_rc%lnc_b(LDT_rc%nnest))
    allocate(LDT_rc%lnr_b(LDT_rc%nnest))
    allocate(LDT_rc%lnc_red(LDT_rc%nnest))
    allocate(LDT_rc%lnr_red(LDT_rc%nnest))
    allocate(LDT_rc%gnc(LDT_rc%nnest))
    allocate(LDT_rc%gnr(LDT_rc%nnest))

    allocate(LDT_rc%glbntiles(LDT_rc%nnest))
    allocate(LDT_rc%glbntiles_red(LDT_rc%nnest))
    allocate(LDT_rc%glbngrid(LDT_rc%nnest))
    allocate(LDT_rc%glbngrid_red(LDT_rc%nnest))
    allocate(LDT_rc%ncatg(LDT_rc%nnest))

    allocate(LDT_rc%glbnpatch(LDT_rc%nnest,LDT_rc%max_model_types))
    allocate(LDT_rc%glbnpatch_red(LDT_rc%nnest,LDT_rc%max_model_types))

    LDT_rc%gridDesc = 0 

    LDT_DElayout = ESMF_DELayoutCreate(LDT_vm,&
         deCountList=(/LDT_rc%npesx,LDT_rc%npesy/),&
         rc=status)

  end subroutine LDT_configinit

!BOP
! !ROUTINE: LDT_coreInit
! \label{LDT_coreInit}
!
! !INTERFACE: 
  subroutine LDT_coreInit

! !DESCRIPTION: 
!
! This "init" routine completes the initialization of the 
!  LDT's clock and alarms.
!
!EOP

! Checks to ensure that the LDT timestep and starting date
!  do not conflict:

!    if(LDT_rc%ts.gt.3600.and. LDT_rc%shr.ne.0) then 
!       write(LDT_logunit,*) 'The starting hour has a permanent offset to the '
!       write(LDT_logunit,*) 'LDT timestep ',LDT_rc%ts
!       write(LDT_logunit,*) 'Please adjust the starting hour in ldt.config'
!       call LDT_endrun()
!    elseif(LDT_rc%ts.gt.60.and.LDT_rc%smn.ne.0) then 
!       write(LDT_logunit,*) 'The starting minute has a permanent offset to the '
!       write(LDT_logunit,*) 'LDT timestep ',LDT_rc%ts
!       write(LDT_logunit,*) 'Please adjust the starting minute in lis.config'
!       call LDT_endrun()
!    endif

    call LDT_update_clock(LDT_rc%ts)

!    call LDT_finishDekadalAlarms(LDT_rc)

  end subroutine LDT_coreInit

!BOP
! !ROUTINE: LDT_ticktime
! \label{LDT_ticktime}
! 
! !INTERFACE: 
  subroutine LDT_ticktime()
!
! !DESCRIPTION: 
!
!  This routine calls the time manager to increment the runtime 
!  clock by the model timestep. 
!
!  The Calling sequence is : 
!  \begin{description}
!   \item[LDT\_advance\_timestep] (\ref{LDT_advance_timestep}) \newline
!    advances the clock 
!   \end{description}
!EOP

    call LDT_advance_timestep(LDT_rc)

!    if(mod(real(LDT_rc%hr)*3600+60*real(LDT_rc%mn)+float(LDT_rc%ss),&
!            real(LDT_rc%tavgInterval)).eq.0) then        
!       LDT_rc%computeFlag = .true.
!    else
!       LDT_rc%computeFlag = .false. 
!    endif

  end subroutine LDT_ticktime


!BOP
! !ROUTINE: LDT_endofrun
! \label{LDT_endofrun}
!      
! !INTERFACE:            
  function LDT_endofrun() result(finish)

! !ARGUMENTS:
    logical :: finish
! !DESCRIPTION:
!  This function checks to see if the runtime clock has reached the
!  specified stop time of the simulation. 
! 
!   The arguments are:
!   \begin{description}
!   \item [finish]
!     boolean value indicating if the end of simulation is reached. 
!   \end{description}
!    
!  The calling sequence is: 
!  \begin{description}
!   \item[LDT\_is\_last\_step] (\ref{LDT_is_last_step}) 
!    check if the clock has reached the stop time
!   \end{description}
!EOP
    integer :: ierr

    if(LDT_masterproc) then
       finish = LDT_is_last_step(LDT_rc)
    endif
#if (defined SPMD)      
    call MPI_BCAST(finish, 1, MPI_LOGICAL, 0, & 
         MPI_COMM_WORLD, ierr)
#endif
  end function LDT_endofrun  

!BOP
! !ROUTINE: LDT_endofTimeWindow
! \label{LDT_endofTimeWindow}
!      
! !INTERFACE:            
  function LDT_endofTimeWindow() result(finish)
! !USES:
    use LDT_logMod, only : LDT_logunit

! !ARGUMENTS:
    logical :: finish
! !DESCRIPTION:
!  This function checks to see if the runtime clock has reached the
!  specified stop time of the simulation. 
! 
!   The arguments are:
!   \begin{description}
!   \item [finish]
!     boolean value indicating if the end of simulation is reached. 
!   \end{description}
!    
!  The calling sequence is: 
!  \begin{description}
!   \item[LDT\_is\_last\_step] (\ref{LDT_is_last_step}) \newline
!    check if the clock has reached the stop time
!   \end{description}
!EOP
    integer         :: rc, ierr
    type(ESMF_Time) :: currTime

    finish = .false.
    if(LDT_masterproc) then
       call ESMF_ClockGet(LDT_clock,currTime=currTime, rc=rc)
       if(currTime.ge.LDT_twStopTime) then  
          finish = .true.
       endif
    endif
#if (defined SPMD)      
    call MPI_BCAST(finish, 1, MPI_LOGICAL, 0, &
         MPI_COMM_WORLD, ierr)
#endif
  end function LDT_endofTimeWindow


!BOP
! !ROUTINE: LDT_finalize
! \label{LDT_finalize}
!      
! !INTERFACE:            
  subroutine LDT_finalize

! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures. 
! 
! The calling sequence is: 
! \begin{description} 
!  \item[spmd\_finalize] (\ref{spmd_finalize}) \newline
!   cleanup SPMD structures
! \end{description}
!EOP
    call spmd_finalize
  end subroutine LDT_finalize



!BOP
! !ROUTINE: spmd_init_offline
! \label{spmd_init_offline} 
!
! !INTERFACE:
  subroutine spmd_init_offline(ldtvm)
!
! !DESCRIPTION:
!  Initializes MPI and retrieves the number of CPUs and the 
!  processors IDs. The MPI initialization is done in an
!  offline mode. In a coupled mode, the initialized
!  environment is passed to LDT from a parent component. 
! 
!EOP
    implicit none
    type(ESMF_VM)   :: ldtvm
    integer         :: ier        ! return error status      '

    call ESMF_Initialize(vm=ldtvm,&
         defaultCalKind=ESMF_CALKIND_GREGORIAN,&
         logkindflag=ESMF_LOGKIND_NONE,rc=ier)

    call ESMF_VMGet(ldtvm,localPet=LDT_localPet,petCount=LDT_npes,&
         rc=ier)

    if (LDT_localPet==0) then 
       LDT_masterproc = .true.
    else
       LDT_masterproc = .false.
    end if

  end subroutine spmd_init_offline

!BOP
! 
! !ROUTINE: spmd_setup
! \label{spmd_setup} 
! 
! !INTERFACE: 
  subroutine spmd_setup(nnest,nmodels)       ! LDT
!  subroutine spmd_setup(nnest,ndas,nmodels)   ! LIS

    implicit none
! !ARGUMENTS: 
    integer, intent(in):: nnest
!    integer, intent(in):: ndas   ! Found in LIS
    integer, intent(in):: nmodels
! 
! !DESCRIPTION: 
!  Allocates memory for the variables describing 
!  domain decomposition. 
!  
! The arguments are: 
! \begin{description}
! \item[nnest]
!  Number of nests or domains
! \end{description}
!EOP
    allocate(LDT_ews_ind(nnest,LDT_npes))
    allocate(LDT_ewe_ind(nnest,LDT_npes))
    allocate(LDT_nss_ind(nnest,LDT_npes))
    allocate(LDT_nse_ind(nnest,LDT_npes))

    allocate(LDT_ews_halo_ind(nnest,LDT_npes))
    allocate(LDT_ewe_halo_ind(nnest,LDT_npes))
    allocate(LDT_nss_halo_ind(nnest,LDT_npes))
    allocate(LDT_nse_halo_ind(nnest,LDT_npes))

    allocate(LDT_ews_b_ind(nnest,LDT_npes))
    allocate(LDT_ewe_b_ind(nnest,LDT_npes))
    allocate(LDT_nss_b_ind(nnest,LDT_npes))
    allocate(LDT_nse_b_ind(nnest,LDT_npes))

    allocate(LDT_deltas(nnest,0:LDT_npes-1))
    allocate(LDT_offsets(nnest,0:LDT_npes-1))
    allocate(LDT_tdeltas(nnest,0:LDT_npes-1))
    allocate(LDT_toffsets(nnest,0:LDT_npes-1))
    allocate(LDT_gdeltas(nnest,0:LDT_npes-1))
    allocate(LDT_goffsets(nnest,0:LDT_npes-1))
    allocate(LDT_odeltas(nnest,0:LDT_npes-1))  ! Original LDT
    allocate(LDT_ooffsets(nnest,0:LDT_npes-1))
!    allocate(LDT_odeltas(nnest,ndas, 0:LDT_npes-1))  ! Reflects latest LIS-7.1
!    allocate(LDT_ooffsets(nnest,ndas, 0:LDT_npes-1))

    allocate(LDT_ntiless(nnest,0:LDT_npes-1))
    allocate(LDT_ngrids(nnest,0:LDT_npes-1))
    
    allocate(LDT_npatches(nnest,nmodels,0:LDT_npes-1))
    allocate(LDT_patch_offsets(nnest,nmodels,0:LDT_npes-1))
    allocate(LDT_patch_deltas(nnest,nmodels,0:LDT_npes-1))

    allocate(LDT_vecTile(nnest))
    allocate(LDT_vecPatch(nnest,nmodels))
    allocate(LDT_vecGrid(nnest))
    allocate(LDT_ensOnGrid(nnest))

    LDT_ews_ind = 0
    LDT_ewe_ind = 0
    LDT_nss_ind = 0 
    LDT_nse_ind = 0 

    LDT_ews_halo_ind = 0
    LDT_ewe_halo_ind = 0
    LDT_nss_halo_ind = 0 
    LDT_nse_halo_ind = 0 

    LDT_ews_b_ind = 0
    LDT_ewe_b_ind = 0
    LDT_nss_b_ind = 0
    LDT_nse_b_ind = 0

    LDT_deltas = 0
    LDT_offsets = 0
    LDT_tdeltas = 0
    LDT_toffsets = 0
    LDT_gdeltas = 0
    LDT_goffsets = 0
    LDT_odeltas = 0
    LDT_ooffsets = 0
    LDT_ntiless = 0
    LDT_ngrids = 0

    LDT_npatches = 0 
    LDT_patch_offsets = 0 
    LDT_patch_deltas = 0 

  end subroutine spmd_setup

!BOP
! !ROUTINE: spmd_finalize
! \label{spmd_finalize} 
!
! !INTERFACE: 
  subroutine spmd_finalize
!
! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures. 
! 
!EOP
    implicit none
    integer :: ierr

    deallocate(LDT_ews_ind)
    deallocate(LDT_ewe_ind)
    deallocate(LDT_nss_ind)
    deallocate(LDT_nse_ind)

    deallocate(LDT_ews_halo_ind)
    deallocate(LDT_ewe_halo_ind)
    deallocate(LDT_nss_halo_ind)
    deallocate(LDT_nse_halo_ind)

    deallocate(LDT_ews_b_ind)
    deallocate(LDT_ewe_b_ind)
    deallocate(LDT_nss_b_ind)
    deallocate(LDT_nse_b_ind)

    deallocate(LDT_deltas)
    deallocate(LDT_offsets)
    deallocate(LDT_tdeltas)
    deallocate(LDT_toffsets)
    deallocate(LDT_gdeltas)
    deallocate(LDT_goffsets)
    deallocate(LDT_odeltas)
    deallocate(LDT_ooffsets)
    deallocate(LDT_ntiless)
    deallocate(LDT_ngrids)
    deallocate(LDT_vecTile)
    deallocate(LDT_vecPatch)
    deallocate(LDT_vecGrid)
    deallocate(LDT_ensOnGrid)

    deallocate(LDT_npatches)
    deallocate(LDT_patch_offsets)
    deallocate(LDT_patch_deltas)
#if ( defined COUPLED) 
#else
#if ( defined SPMD )
    call MPI_FINALIZE(ierr)
#endif
#endif
  end subroutine spmd_finalize

!BOP
! 
! !ROUTINE: LDT_isLDTatAfinerResolution
! 
! !INTERFACE:
  function LDT_isLDTatAfinerResolution(n,datares) result(finish)
!
! !DESCRIPTION:
! 
!  This function checks if the given resolution is finer than 
!  the LIS/target resolution. 
! 
!EOP
    integer :: n 
    real    :: datares
! !ARGUMENTS:
    logical :: finish
    
    if(LDT_rc%gridDesc(n,1).eq.0) then     ! lat-lon grid
       if(LDT_rc%gridDesc(n,10).lt.datares) then 
          finish = .true.
       else
          finish = .false. 
       endif
    elseif(LDT_rc%gridDesc(n,1).eq.9) then  ! EASE-grid
       if(LDT_rc%gridDesc(n,10).lt.datares) then 
          finish = .true.
       else
          finish = .false. 
       endif
    elseif(LDT_rc%gridDesc(n,1).eq.3) then  ! lambert grid
       if(LDT_rc%gridDesc(n,8)/100.0.lt.datares) then 
          finish = .true.
       else
          finish = .false. 
       endif
    else
       write(LDT_logunit,*) '[ERR] LDT_isLDTatAfinerResolution check not supported '
       write(LDT_logunit,*) '[ERR]  for this map projection.'
       call LDT_endrun()
    endif

  end function LDT_isLDTatAfinerResolution

end module LDT_coreMod
