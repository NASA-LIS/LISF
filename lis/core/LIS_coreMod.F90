!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
module LIS_coreMod
!BOP
!
! !MODULE: LIS_coreMod
!
! !DESCRIPTION:
!  The code in this file contains the basic datastructures and
!  control routines for the operation of LIS
!
!  \subsubsection{Overview}
!  This module contains the defintion and specification of the basic
!  datastructures that define a LIS instance. It consists of:
!  \begin{description}
!   \item[LIS\_rc]
!    datastructure containg generic LIS run control (rc) variables
!   \item[LIS\_domain]
!    datastructure containing grid, tile spaces,
!    grid projection, domain decomposition information.
!   \item[LIS\_histData]
!     datastructure defining the metadata for
!     model output.
!   \item[LIS\_config]
!     instance of the configuration class.
!   \item[LIS\_vm]
!     instance of the virtual machine controlling
!     the LIS resources.
!   \item[LIS\_masterproc]
!    logical variable to specify the master processor (true = for master,
!    false for non-zero processors)
!   \item[LIS\_npes]
!    number of processors used in LIS
!   \item[LIS\_ews\_ind]
!    starting East West index of each nest (excluding halo regions)
!   \item[LIS\_ewe\_ind]
!    ending East West index of each nest (excluding halo regions)
!   \item[LIS\_nss\_ind]
!    starting North-South index of each nest (excluding halo regions)
!   \item[LIS\_nse\_ind]
!    ending North-South index of each nest (excluding halo regions)
!   \item[LIS\_ews\_ind\_halo]
!    starting East West index of each nest (including halo regions)
!   \item[LIS\_ewe\_ind\_halo]
!    ending East West index of each nest (including halo regions)
!   \item[LIS\_nss\_ind\_halo]
!    starting North-South index of each nest (including halo regions)
!   \item[LIS\_nse\_ind\_halo]
!    ending North-South index of each nest (including halo regions)
!   \item[LIS\_deltas]
!    size of unmasked grid space buffers
!   \item[LIS\_offsets]
!   offsets for unmasked grid space buffers
!   \item[LIS\_tdeltas]
!   size of tile space buffers
!   \item[LIS\_toffsets]
!    offsets of tile space buffers
!   \item[LIS\_gdeltas]
!    size of grid space buffers
!   \item[LIS\_goffsets]
!    offsets of grid space buffers
!   \item[LIS\_ntiless]
!    size of tile spaces
!    \item[LIS\_ngrids]
!    size of grid spaces
!   \item[LIS\_vecTile]
!    grid object for the tile space-vector of tiles
!   \item[LIS\_vecGrid]
!    grid object for the grid space-vector of grid pts
!  \end{description}
!
! !REVISION HISTORY:
!  14 Nov 2002    Sujay Kumar  Initial Specification
!
! !USES:
  use ESMF
  use LIS_PRIV_rcMod
  use LIS_PRIV_gridMod
  use LIS_PRIV_tileMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_mpiMod
  use map_utils

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_config_init
  public :: LIS_core_init
  public :: LIS_ticktime
  public :: LIS_endofrun
  public :: LIS_endOfTimeWindow
  public :: LIS_TimeToRunNest
  public :: LIS_resetTimeMgr
  public :: LIS_finalize
  public :: LIS_isatAfinerResolution
  public :: LIS_howtoTransform
  public :: LIS_getDomainResolutions
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LIS_rc
  public :: LIS_domain
  public :: LIS_surface
  public :: LIS_config
  public :: LIS_initialized
  public :: LIS_vm
  public :: LIS_masterproc
  public :: LIS_localPet
  public :: LIS_npes
  public :: LIS_ews_ind
  public :: LIS_ewe_ind
  public :: LIS_nss_ind
  public :: LIS_nse_ind
  public :: LIS_ews_halo_ind
  public :: LIS_ewe_halo_ind
  public :: LIS_nss_halo_ind
  public :: LIS_nse_halo_ind
  public :: LIS_ews_b_ind
  public :: LIS_ewe_b_ind
  public :: LIS_nss_b_ind
  public :: LIS_nse_b_ind
  public :: LIS_deltas
  public :: LIS_offsets
  public :: LIS_tdeltas
  public :: LIS_toffsets
  public :: LIS_gdeltas
  public :: LIS_goffsets
  public :: LIS_patch_deltas
  public :: LIS_patch_offsets
  public :: LIS_ntiless
  public :: LIS_ngrids
  public :: LIS_npatches
! Three different types of grid objects are defined because ESMF3.1.0 does
! not support defining tiles on a single grid.
  public :: LIS_vecTile
  public :: LIS_vecPatch
  public :: LIS_vecGrid



!EOP
  type, public :: lis_domain_type
     real,          allocatable :: glat(:)
     real,          allocatable :: glon(:)
     real,          allocatable :: lat(:)
     real,          allocatable :: lon(:)
     real,          allocatable :: lat_b(:)
     real,          allocatable :: lon_b(:)
     type(griddec), allocatable :: grid(:)
     type(tiledec), allocatable :: tile(:)
     type(proj_info)            :: lisproj
     integer, allocatable       :: gindex(:,:)
     integer, allocatable       :: ntiles_pergrid(:)
     integer, allocatable       :: str_tind(:) !starting tile id for each grid
     real                       :: minLat, maxLat, minLon, maxLon
     real                       :: stlat, stlon, truelat1
     real                       :: truelat2, truelon, orient
     real                       :: dx,dy,nlatcircles
  end type lis_domain_type

  type, public :: lis_domain_sf_type
     type(tiledec), allocatable :: tile(:)
     integer,       allocatable :: npatch_pergrid(:)
     integer,       allocatable :: str_patch_ind(:)
  end type lis_domain_sf_type

  type(lisrcdec), save           :: LIS_rc
  type(lis_domain_type), allocatable :: LIS_domain(:)
  type(lis_domain_sf_type), allocatable :: LIS_surface(:,:)
  logical                        :: LIS_initialized = .false.
  type(ESMF_Config), save        :: LIS_config
  type(ESMF_VM), save            :: LIS_vm
  logical                        :: LIS_masterproc
  integer                        :: LIS_localPet
  integer                        :: LIS_npes
  logical                        :: LIS_init_esmf

  integer, allocatable           :: LIS_ews_ind(:,:),LIS_ewe_ind(:,:)
  integer, allocatable           :: LIS_nss_ind(:,:),LIS_nse_ind(:,:)
  integer, allocatable           :: LIS_ews_halo_ind(:,:),LIS_ewe_halo_ind(:,:)
  integer, allocatable           :: LIS_nss_halo_ind(:,:),LIS_nse_halo_ind(:,:)
  integer, allocatable           :: LIS_ews_b_ind(:,:),LIS_ewe_b_ind(:,:)
  integer, allocatable           :: LIS_nss_b_ind(:,:),LIS_nse_b_ind(:,:)
  integer, allocatable           :: LIS_deltas(:,:),LIS_offsets(:,:)
  integer, allocatable           :: LIS_tdeltas(:,:),LIS_toffsets(:,:)
  integer, allocatable           :: LIS_gdeltas(:,:),LIS_goffsets(:,:)
  integer, allocatable           :: LIS_ntiless(:,:),LIS_ngrids(:,:)
  integer, allocatable           :: LIS_npatches(:,:,:)
  integer, allocatable           :: LIS_patch_deltas(:,:,:),LIS_patch_offsets(:,:,:)

  type(ESMF_Grid),     allocatable   :: LIS_vecTile(:)
  type(ESMF_Grid),     allocatable   :: LIS_vecPatch(:,:)
  type(ESMF_Grid),     allocatable   :: LIS_vecGrid(:)

!BOPI
! !ROUTINE: LIS_config_init
! \label{LIS_config_init}
!
! !INTERFACE:
  interface LIS_config_init
! !PRIVATE MEMBER FUNCTIONS:
     module procedure lisconfig_generic
! !DESCRIPTION:
!   This interface provides routines for the setup of LIS
!   configuration management. It initializes the LIS Config utility,
!   reads the model independent specifications, and intializes the SPMD
!   parallel processing mode. This routine also initializes the
!   LIS time manager, and finally allocates memory for the LIS
!   generic datastructures.
!EOPI
  end interface

contains
!BOP
! !ROUTINE: lisconfig_generic
! \label{lisconfig_generic}
!
! !INTERFACE:
! Private name: call using LIS_config_init()
  subroutine lisconfig_generic(cmd_args, vm, clock, nx, ny, comm)

    implicit none
! !ARGUMENTS:
    logical, intent(in), optional         :: cmd_args
    type(ESMF_VM), intent(in), optional   :: vm
    type(ESMF_Clock), intent(in),optional :: clock
    integer, intent(in), optional         :: nx, ny
    integer, intent(in), optional         :: comm
!
! !DESCRIPTION:
!  Performs the initialization of LIS runtime configuration.
!
! The arguments are:
! \begin{description}
! \item[cmd\_args]
!   optional logical flag to indictate whether to process command line
!   arguments.  If present and .true., LIS will process command line
!   arguments.
! \item[vm]
!   optional virtual machine.  If present, LIS will use this virtual
!   machine instead of initializing a virtual machine via a call to
!   ESMF\_Initialize.  (Do not use with \var{nx} and \var{ny} or \var{comm}.)
! \item[clock]
!   optional ESMF clock.  Currently not used.
! \item[nx]
!   optional value representing number of processes in the x-direction.
!   If present, LIS will use this value for \var{LIS\_rc\%npesx} instead of
!   the value specified by \var{Number of processors along x:}
!   in the \file{lis.config} file.  Must be paired with \var{ny}.
!   (Do not use with \var{vm}.)
! \item[ny]
!   optional value representing number of processes in the y-direction.
!   If present, LIS will use this value for \var{LIS\_rc\%npesy} instead of
!   the value specified by \var{Number of processors along y:}
!   in the \file{lis.config} file.  Must be paired with \var{nx}.
!   (Do not use with \var{vm}.)
! \item[comm]
!   optional value representing the MPI communicator to use.
!   If present, LIS will use this MPI communicator
!   when initializing a virtual machine via a call to
!   ESMF\_Initialize.
!   (Do not use with \var{vm}.)
! \end{description}
!
!  The Calling sequence is :
!  \begin{description}
!    \item[spmd\_init\_generic] (\ref{spmd_init_generic}) \newline
!     performs SPMD initializations
!    \item[LIS\_log\_init] (\ref{LIS_log_init}) \newline
!     initializes the LIS log handler
!    \item[lis\_process\_cmd\_args] (\ref{LIS_process_cmd_args}) \newline
!     processes any command line arguments to LIS
!   \item[LIS\_readConfig] (\ref{LIS_readConfig}) \newline
!     reads the model independent options
!   \item[spmd\_setup] (\ref{spmd_setup}) \newline
!     allocates memory for SPMD variables
!   \item[LIS\_timemgr\_init] (\ref{LIS_timemgr_init}) \newline
!    initializes the time manager
!   \end{description}
!EOP

    integer :: status

    if ( present(vm) .and. &
         ( present(comm) .or. present(nx) .or. present(ny) ) ) then
       write(*,*) "[ERR] Do not mix vm with comm or nx, ny."
       write(*,*) "[ERR] Aborting."
       call LIS_endrun()
    endif

    if ( present(vm) ) then
       ! use incoming VM and its MPI communicator for LIS
       LIS_vm = vm
       LIS_init_esmf = .false.
    else
       LIS_init_esmf = .true.
    endif

    if ( present(comm) ) then
       call spmd_init_generic(LIS_vm, LIS_init_esmf, comm)
    else
       call spmd_init_generic(LIS_vm, LIS_init_esmf)
    endif

    call lis_log_init(LIS_getNextUnitNumber())

    if ( present(cmd_args) ) then
       if ( cmd_args .eqv. .true. ) then
          call lis_process_cmd_args
       endif
    endif

    call LIS_readConfig

    call spmd_setup(LIS_rc%nnest, &
         LIS_rc%max_model_types)

    call LIS_timemgr_init(LIS_rc)
    !if ( present(clock) ) then
    !overwrite the default clock:
    !LIS_clock = clock <-- this is not the proper way to copy a clock
    !LIS_clock = ESMF_ClockCreate(clock) <-- this copies the clock
    !I am not overwriting LIS' clock.
    !endif

    LIS_rc%endtime = 0

    if ( present(nx) .and. present(ny) ) then
       LIS_rc%npesx = nx
       LIS_rc%npesy = ny
    endif

    allocate(LIS_domain(LIS_rc%nnest))
    allocate(LIS_surface(LIS_rc%nnest, LIS_rc%max_model_types))
    allocate(LIS_rc%ntiles(LIS_rc%nnest))
    allocate(LIS_rc%ngrid(LIS_rc%nnest))
    allocate(LIS_rc%lnc(LIS_rc%nnest))
    allocate(LIS_rc%lnr(LIS_rc%nnest))
    allocate(LIS_rc%lnc_b(LIS_rc%nnest))
    allocate(LIS_rc%lnr_b(LIS_rc%nnest))
    allocate(LIS_rc%lnc_red(LIS_rc%nnest))
    allocate(LIS_rc%lnr_red(LIS_rc%nnest))
    allocate(LIS_rc%gnc(LIS_rc%nnest))
    allocate(LIS_rc%gnr(LIS_rc%nnest))
    allocate(LIS_rc%gnc_b(LIS_rc%nnest))
    allocate(LIS_rc%gnr_b(LIS_rc%nnest))
    allocate(LIS_rc%gridDesc(LIS_rc%nnest,50))
    allocate(LIS_rc%glbntiles(LIS_rc%nnest))
    allocate(LIS_rc%glbntiles_red(LIS_rc%nnest))
    allocate(LIS_rc%glbngrid(LIS_rc%nnest))
    allocate(LIS_rc%glbngrid_red(LIS_rc%nnest))
    allocate(LIS_rc%ncatg(LIS_rc%nnest))

    allocate(LIS_rc%glbnpatch(LIS_rc%nnest,LIS_rc%max_model_types))
    allocate(LIS_rc%glbnpatch_red(LIS_rc%nnest,LIS_rc%max_model_types))

    allocate(LIS_rc%minLat(LIS_rc%nnest))
    allocate(LIS_rc%minLon(LIS_rc%nnest))
    allocate(LIS_rc%maxLat(LIS_rc%nnest))
    allocate(LIS_rc%maxLon(LIS_rc%nnest))
    
    LIS_rc%gridDesc = 0
    LIS_rc%use_twelve = .false.
    LIS_rc%reset_flag = .false.
    LIS_rc%run_model  = .true.

  end subroutine lisconfig_generic

!BOP
! !ROUTINE: LIS_core_init
! \label{LIS_core_init}
!
! !INTERFACE:
  subroutine LIS_core_init
! !DESCRIPTION:
!
! Completes the initialization of the LIS' clock and alarms.
!EOP

    TRACE_ENTER("core_init")
!checks to ensure that the LIS timestep and starting date does not
!conflict

!    if(LIS_rc%ts.gt.3600.and. LIS_rc%shr.ne.0) then
!       write(LIS_logunit,*) 'The starting hour has a permanent offset to the '
!       write(LIS_logunit,*) 'LIS timestep ',LIS_rc%ts
!       write(LIS_logunit,*) 'Please adjust the starting hour in lis.config'
!       call LIS_endrun()
!    elseif(LIS_rc%ts.gt.60.and.LIS_rc%smn.ne.0) then
!       write(LIS_logunit,*) 'The starting minute has a permanent offset to the '
!       write(LIS_logunit,*) 'LIS timestep ',LIS_rc%ts
!       write(LIS_logunit,*) 'Please adjust the starting minute in lis.config'
!       call LIS_endrun()
!    endif
    call LIS_update_clock(LIS_rc%ts)
    call LIS_finishDekadalAlarms(LIS_rc)
    
    call LIS_updateAlarmSetups()

    write(LIS_logunit,*) '[INFO]: LIS timestep ',LIS_rc%ts
    TRACE_EXIT("core_init")
  end subroutine LIS_core_init

!BOP
! !ROUTINE: LIS_ticktime
! \label{LIS_ticktime}
!
! !INTERFACE:
  subroutine LIS_ticktime()
!
! !DESCRIPTION:
!  This routine calls the time manager to increment the runtime
!  clock by the model timestep.
!
!  The Calling sequence is :
!  \begin{description}
!   \item[LIS\_advance\_timestep] (\ref{LIS_advance_timestep}) \newline
!    advances the clock
!   \end{description}
!EOP
    call LIS_advance_timestep(LIS_rc)
  end subroutine LIS_ticktime


!BOP
! !ROUTINE: LIS_endofrun
! \label{LIS_endofrun}
!
! !INTERFACE:
  function LIS_endofrun() result(finish)

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
!   \item[LIS\_is\_last\_step] (\ref{LIS_is_last_step}) \newline
!    check if the clock has reached the stop time
!   \end{description}
!EOP
    integer :: ierr

    if(LIS_masterproc) then
       finish = LIS_is_last_step(LIS_rc)
    endif
#if (defined SPMD)
    call MPI_BCAST(finish, 1, MPI_LOGICAL, 0, &
         LIS_mpi_comm, ierr)
#endif
  end function LIS_endofrun

!BOP
! !ROUTINE: LIS_endofTimeWindow
! \label{LIS_endofTimeWindow}
!
! !INTERFACE:
  function LIS_endofTimeWindow() result(finish)
! !USES:
    use LIS_logMod, only : LIS_logunit
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
!   \item[LIS\_is\_last\_step] (\ref{LIS_is_last_step}) \newline
!    check if the clock has reached the stop time
!   \end{description}
!EOP
    integer         :: rc, ierr
    type(ESMF_Time) :: currTime


    finish = .false.

    if(LIS_masterproc) then
       call ESMF_ClockGet(LIS_clock,currTime=currTime, rc=rc)
       if(currTime.ge.LIS_twStopTime) then  !BZ changed gt to ge
          finish = .true.
       endif
! Assuming that all nests are in sync
!       if(LIS_rc%endtime.eq.1.and.LIS_rc%DAincrMode(1).eq.1) then 
!          finish = .true. 
!       endif

    endif
#if (defined SPMD)
    call MPI_BCAST(finish, 1, MPI_LOGICAL, 0, &
         LIS_mpi_comm, ierr)
#endif
  end function LIS_endofTimeWindow

!BOP
! !ROUTINE: LIS_timeToRunNest
! \label{LIS_timeToRunNest}
!
! !INTERFACE:
  function LIS_timeToRunNest(n) result(check)
! !ARGUMENTS:
    integer, intent(in) :: n
    logical :: check
!
! !DESCRIPTION:
!  Check to see if it is time to run the nest based on the
!  model timestep for that particular domain
!
!   \begin{description}
!   \item [n]
!      index of the nest or domain
!   \item [check]
!     boolean value indicating if the time to run the nest has
!     reached or not
!   \end{description}
!
!EOP
    real    :: curr_time

!multiple timesteps for different nests are not supported
!if time step is > 1 day
    if(LIS_rc%nts(n).ge.86400) then
       check = .true.
    else
       curr_time = float(LIS_rc%hr)*3600+60*float(LIS_rc%mn)+float(LIS_rc%ss)
       if(mod(curr_time,real(LIS_rc%nts(n))).eq.0) then
          check = .true.
       else
          check = .false.
       endif
    endif

  end function LIS_timeToRunNest

!BOP
! !ROUTINE: LIS_resetTimeMgr
! \label{LIS_resetTimeMgr}
!
! !INTERFACE:
  subroutine LIS_resetTimeMgr
!
! !DESCRIPTION:
!   This routine resets the time manager clock back to the starting
!    time.
!
! The calling sequence is:
! \begin{description}
!  \item[LIS\_resetclock] (\ref{LIS_resetclock})
! \end{description}
!EOP
    call LIS_resetclock(LIS_rc)
  end subroutine LIS_resetTimeMgr

!BOP
! !ROUTINE: LIS_finalize
! \label{LIS_finalize}
!
! !INTERFACE:
  subroutine LIS_finalize

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
    call spmd_finalize(LIS_init_esmf)
  end subroutine LIS_finalize


!BOP
! !ROUTINE: spmd_init_generic
! \label{spmd_init_generic}
!
! !INTERFACE:
  subroutine spmd_init_generic(lisvm, init_esmf, liscomm)
!
! !DESCRIPTION:
!  Determines whether to initialize an ESMF virtual machine or to use a
!  given virtual machine.  Likewise, determines whether to initialize
!  an MPI communicator or to use a given communicator.
!  Then retrieves the number of CPUs and the
!  processors IDs.
!
!   The arguments are:
!   \begin{description}
!   \item [lisvm]
!     ESMF virtual machine.  If given, LIS will use this virtual machine.
!     Else, LIS will create this virtual machine.
!   \item [init\_esmf]
!     boolean.  Specifies whether to call ESMF\_Initialize.  If \var{lisvm}
!     is given then \var{init\_esmf} must be set to \var{.false.} in the
!     calling routine.  Else, it must be set to \var{.true.}.
!   \item [liscomm]
!     optional MPI communicator.  If given, LIS will use this communicator.
!     Else, LIS will use the communicator from the lisvm virtual machine.
!   \end{description}
!
!   If \var{liscomm} is given, then this routine assumes that a virtual
!   machine must be created.  Do not use both an initialized ESMF virtual
!   machine and an initialized MPI communicator with this routine.
!
!EOP
    implicit none
    type(ESMF_VM)     :: lisvm
    logical           :: init_esmf
    integer, optional :: liscomm

    integer           :: ier        ! return error status


    ! If liscomm is present, then assume that a virtual machine
    ! must be initialized.  Use the given liscomm as the MPI communicator.
    !
    ! If liscomm is not present, then either initialize or use the lisvm
    ! virtual machine based on the init_esmf flag.

    if ( present(liscomm) ) then
       LIS_mpi_comm = liscomm

       call ESMF_Initialize(vm=lisvm,&
            mpiCommunicator=liscomm,&
            defaultCalKind=ESMF_CALKIND_GREGORIAN,&
            logkindflag=ESMF_LOGKIND_NONE,rc=ier)
    else
       if ( init_esmf ) then
          call ESMF_Initialize(vm=lisvm,&
               defaultCalKind=ESMF_CALKIND_GREGORIAN,&
               logkindflag=ESMF_LOGKIND_NONE,rc=ier)
       endif
    endif

    call ESMF_VMGet(lisvm,localPet=LIS_localPet,petCount=LIS_npes,&
                    mpiCommunicator=LIS_mpi_comm,rc=ier)
    if (ier /= ESMF_SUCCESS) then
       print *, "LIS: spmd_init_generic: error getting VM information."
    endif

    if (LIS_localPet==0) then
       LIS_masterproc = .true.
    else
       LIS_masterproc = .false.
    endif

  end subroutine spmd_init_generic


!BOP
! !ROUTINE: spmd_init_coupled
! \label{spmd_init_coupled}
!
! !INTERFACE:
  subroutine spmd_init_coupled()
!
! !DESCRIPTION:
!  Obtains the number of processors and the processors IDs from
!  an already initialized MPI state. The initialization is assumed
!  to be performed in a parent component.
!
!EOP
    implicit none
    integer ier        ! return error status

#if (defined SPMD)
    call mpi_comm_rank(LIS_mpi_comm, LIS_localPet, ier)
    call mpi_comm_size(LIS_mpi_comm, LIS_npes, ier)
#else
    LIS_localPet = 0
    LIS_npes = 1
#endif
    if (LIS_localPet==0) then
       LIS_masterproc = .true.
    else
       LIS_masterproc = .false.
    end if
    return

  end subroutine spmd_init_coupled

!BOP
!
! !ROUTINE: spmd_setup
! \label{spmd_setup}
!
! !INTERFACE:
  subroutine spmd_setup(nnest, nmodels)

    implicit none
! !ARGUMENTS:
    integer, intent(in):: nnest
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
    allocate(LIS_ews_ind(nnest,LIS_npes))
    allocate(LIS_ewe_ind(nnest,LIS_npes))
    allocate(LIS_nss_ind(nnest,LIS_npes))
    allocate(LIS_nse_ind(nnest,LIS_npes))

    allocate(LIS_ews_halo_ind(nnest,LIS_npes))
    allocate(LIS_ewe_halo_ind(nnest,LIS_npes))
    allocate(LIS_nss_halo_ind(nnest,LIS_npes))
    allocate(LIS_nse_halo_ind(nnest,LIS_npes))

    allocate(LIS_ews_b_ind(nnest,LIS_npes))
    allocate(LIS_ewe_b_ind(nnest,LIS_npes))
    allocate(LIS_nss_b_ind(nnest,LIS_npes))
    allocate(LIS_nse_b_ind(nnest,LIS_npes))

    allocate(LIS_deltas(nnest,0:LIS_npes-1))
    allocate(LIS_offsets(nnest,0:LIS_npes-1))
    allocate(LIS_tdeltas(nnest,0:LIS_npes-1))
    allocate(LIS_toffsets(nnest,0:LIS_npes-1))
    allocate(LIS_gdeltas(nnest,0:LIS_npes-1))
    allocate(LIS_goffsets(nnest,0:LIS_npes-1))
    allocate(LIS_ntiless(nnest,0:LIS_npes-1))
    allocate(LIS_ngrids(nnest,0:LIS_npes-1))

    allocate(LIS_npatches(nnest,nmodels,0:LIS_npes-1))
    allocate(LIS_patch_offsets(nnest,nmodels,0:LIS_npes-1))
    allocate(LIS_patch_deltas(nnest,nmodels,0:LIS_npes-1))

    allocate(LIS_vecTile(nnest))
    allocate(LIS_vecPatch(nnest,nmodels))
    allocate(LIS_vecGrid(nnest))

    LIS_ews_ind = 0
    LIS_ewe_ind = 0
    LIS_nss_ind = 0
    LIS_nse_ind = 0

    LIS_ews_halo_ind = 0
    LIS_ewe_halo_ind = 0
    LIS_nss_halo_ind = 0
    LIS_nse_halo_ind = 0

    LIS_ews_b_ind = 0
    LIS_ewe_b_ind = 0
    LIS_nss_b_ind = 0
    LIS_nse_b_ind = 0

    LIS_deltas = 0
    LIS_offsets = 0
    LIS_tdeltas = 0
    LIS_toffsets = 0
    LIS_gdeltas = 0
    LIS_goffsets = 0
    LIS_ntiless = 0
    LIS_ngrids = 0

    LIS_npatches = 0
    LIS_patch_offsets = 0
    LIS_patch_deltas = 0
  end subroutine spmd_setup

!BOP
! !ROUTINE: spmd_finalize
! \label{spmd_finalize}
!
! !INTERFACE:
  subroutine spmd_finalize(fin_esmf)
!
! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures.
!
!   The arguments are:
!   \begin{description}
!   \item [fin\_esmf]
!     boolean.  Specifies whether to call ESMF\_Finalize.
!   \end{description}
!
!EOP
    implicit none
    logical :: fin_esmf
    integer :: ierr

    deallocate(LIS_ews_ind)
    deallocate(LIS_ewe_ind)
    deallocate(LIS_nss_ind)
    deallocate(LIS_nse_ind)

    deallocate(LIS_ews_halo_ind)
    deallocate(LIS_ewe_halo_ind)
    deallocate(LIS_nss_halo_ind)
    deallocate(LIS_nse_halo_ind)

    deallocate(LIS_ews_b_ind)
    deallocate(LIS_ewe_b_ind)
    deallocate(LIS_nss_b_ind)
    deallocate(LIS_nse_b_ind)

    deallocate(LIS_deltas)
    deallocate(LIS_offsets)
    deallocate(LIS_tdeltas)
    deallocate(LIS_toffsets)
    deallocate(LIS_gdeltas)
    deallocate(LIS_goffsets)
    deallocate(LIS_ntiless)
    deallocate(LIS_ngrids)
    deallocate(LIS_vecTile)
    deallocate(LIS_vecPatch)
    deallocate(LIS_vecGrid)

    deallocate(LIS_npatches)
    deallocate(LIS_patch_offsets)
    deallocate(LIS_patch_deltas)
    if ( fin_esmf ) then
       call ESMF_Finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)
    endif
#if ( defined COUPLED)
#else
#if ( defined SPMD )
    call MPI_FINALIZE(ierr)
#endif
#endif
  end subroutine spmd_finalize

!BOP
!
! !ROUTINE: LIS_isatAfinerResolution
! \label{LIS_isatAfinerResolution}
!
! !INTERFACE:
  function LIS_isatAfinerResolution(n,datares) result(finish)

  implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
    real, intent(in)    :: datares
!
! !DESCRIPTION:
!
! This function determines whether LIS' running domain resolution
! is finer than the given datares resolution.
!
! This function returns .true. when LIS' running domain resolution
! is finer than datares; otherwise, it returns .false..
!
! Note that this function supports checking only running domains with
! a cylindrical lat/lon projection.  For other projections, this
! function returns ".true.".  Although this is not correct, it does mimic
! the behaviour of LIS before the use of this function.
!
! The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[datares]
!    the resolution to compare against
!  \end{description}
!EOP
    logical :: finish

    real    :: targetres

    if ( LIS_rc%gridDesc(n,1) .eq. 0 ) then
       ! Compare the larger of dlon and dlat against datares.
       if ( LIS_rc%gridDesc(n,9) > LIS_rc%gridDesc(n,10) ) then
          targetres = LIS_rc%gridDesc(n,9)
       else
          targetres = LIS_rc%gridDesc(n,10)
       endif

       if ( targetres .lt. datares ) then
          finish = .true.
       else
          finish = .false.
       endif
    elseif ( LIS_rc%gridDesc(n,1) .eq. 3 ) then !lambert
       targetres = LIS_rc%gridDesc(n,8)/100.0

       if ( targetres .lt. datares ) then
          finish = .true.
       else
          finish = .false.
       endif
       
    else
       write(LIS_logunit,*) '[WARN] : LIS_isatAfinerResolution check ' // &
                            '[WARN] is NOT supported'
       write(LIS_logunit,*) '[WARN] for this map projection.'
       write(LIS_logunit,*) '[WARN] Returning ".true.".'
       !call LIS_endrun()
       finish = .true.
    endif
  end function LIS_isatAfinerResolution

!BOP
!
! !ROUTINE: LIS_howtoTransform
! \label{LIS_howtoTransform}
!
! !INTERFACE:
  function LIS_howtoTransform(n,datares) result(finish)

  implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
    real, intent(in)    :: datares
!
! !DESCRIPTION:
!
! This function determines whether LIS' running domain resolution
! is finer than, comparable to, or coarser than the given datares resolution.
! This is used to determine how to transform data at the given datares
! resolution to LIS' running domain.
!
!
! datares should be the maximum of the dlon and dlat for the input grid.
!
! This function returns "interpolate" when LIS' running domain resolution
! is finer than datares, "upscale" when LIS' running domain resolution
! is coarser, "neighbor" if comparable.
!
! Note that this function supports checking only running domains with
! a cylindrical lat/lon projection.  For other projections, this
! function returns "interpolate".  Although this is not correct, it does mimic
! the behaviour of LIS before the use of this function.
!
! The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[datares]
!    the resolution to compare against
!  \end{description}
!EOP
    character(len=16) :: finish
    real    :: targetres, dom
    
    if ( LIS_rc%gridDesc(n,1) .eq. 0 ) then
       ! Compare the larger of dlon and dlat against datares.
       if ( LIS_rc%gridDesc(n,9) > LIS_rc%gridDesc(n,10) ) then
          targetres = LIS_rc%gridDesc(n,9)
       else
          targetres = LIS_rc%gridDesc(n,10)
       endif

       dom = max(targetres, datares)

       if ( ( abs(targetres - datares) / dom ) <= 0.2 ) then 
          finish = 'neighbor'
       else
          if ( targetres <= datares ) then
             finish = 'interpolate'
          else
             finish = 'upscale'
          endif
       endif
    else
       write(LIS_logunit,*) '[WARN] LIS_howtoTransform check ' // &
                            '[WARN] is NOT supported'
       write(LIS_logunit,*) '[WARN] for this map projection.'
       write(LIS_logunit,*) '[WARN] Returning "interpolate".'
       !call LIS_endrun()
       finish = 'interpolate'
    endif
  end function LIS_howtoTransform


!BOP
!
! !ROUTINE: LIS_getDomainResolutions
!
! !INTERFACE:
  subroutine LIS_getDomainResolutions(n,dx,dy)

  implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
    real                :: dx, dy
!
! !DESCRIPTION:
!
! This function determines spatial resolution of the LIS running grid
! (in degrees)
!
!
! The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[dx]
!    spatial resolution along the east-west dimension
!  \item[dy]
!    spatial resolution along the north-south dimension
!  \end{description}
!EOP

    if ( LIS_rc%gridDesc(n,1) .eq. 0 ) then
       dx = LIS_rc%gridDesc(n,9)
       dy = LIS_rc%gridDesc(n,10)
    elseif ( LIS_rc%gridDesc(n,1) .eq. 1 ) then
       dx = LIS_rc%gridDesc(n,8)/100.0
       dy = LIS_rc%gridDesc(n,9)/100.0
    elseif ( LIS_rc%gridDesc(n,1) .eq. 3 ) then
       dx = LIS_rc%gridDesc(n,8)/100.0
       dy = LIS_rc%gridDesc(n,9)/100.0
    elseif ( LIS_rc%gridDesc(n,1) .eq. 5 ) then
       dx = LIS_rc%gridDesc(n,8)/100.0
       dy = LIS_rc%gridDesc(n,9)/100.0
    else
       
       write(LIS_logunit,*) '[ERR] LIS_getDomainResolutions routine ' // &
                            '[ERR] is NOT supported'
       write(LIS_logunit,*) '[ERR] for this map projection.'
       call LIS_endrun()

    endif
  end subroutine LIS_getDomainResolutions


!BOP
! 
! !ROUTINE: LIS_updateAlarmSetups
! \label{LIS_updateAlarmSetups}
! 
! !INTERFACE: 
  subroutine LIS_updateAlarmSetups()
!
! !DESCRIPTION: 
!   This subroutine registers alarms after 
!   the updated LIS clock timestep has been determined
! 
!EOP

    integer                 :: i
    character*3             :: fda

    do i=1,LIS_rc%ndas

       write(fda,'(i3.3)') i

       call LIS_registerAlarm("LIS DA output "//trim(fda),&
            LIS_rc%ts,&
            LIS_rc%daOutInterval(i))
       
    enddo

    if(LIS_rc%nperts.gt.0) then 
       call LIS_registerAlarm("LIS pert restart alarm",&
            LIS_rc%ts ,&  
            LIS_rc%pertrestartInterval)
    endif
    
  end subroutine LIS_updateAlarmSetups
end module LIS_coreMod
