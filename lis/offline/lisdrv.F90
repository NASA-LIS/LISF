!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#ifdef USE_PFIO
#include "MAPL_ErrLog.h"
#include "unused_dummy.H"
#endif

#include "LIS_misc.h"
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
!BOP
!
! !ROUTINE: lisdrv
!  \label{lisdrv}
!   Main program for LIS in an offline simulation
!  
! !DESCRIPTION: 
!  Main driver program for LIS. It performs four main functions
!
!  \begin{description}
!  \item[LIS\_config\_init]  
!      calls the routines to read the runtime configurations
!  \item[LIS\_Init]  
!      calls the initializations based on the runmode
!  \item[LIS\_run]
!      calls the run routines based on the runmode
!  \item[LIS\_finalize] 
!     calls the cleanup routines 
!  \end{description}
!
!  The routines invoked are :
!  \begin{description}
!   \item[LIS\_config\_init](\ref{LIS_config_init}) \newline
!    call to initialize configuration tool and read model independent
!    options
!   \item[lisinit](\ref{lisinit}) \newline
!    call to initialize lis based on the runmode
!   \item[lisrun](\ref{lisrun}) \newline
!    call to run lis based on the runmode
!   \item[lisfinalize](\ref{lisfinalize}) \newline
!    call to cleanup lis structures based on the runmode
!  \end{description}
! !REVISION HISTORY: 
!  14Nov02    Sujay Kumar  Initial Specification
!  21Oct05    Sujay Kumar  Modified to include the runmodes. Switched
!                          to a init,run,finalize mode
!  13Aug23    Jules Kouatchou Introduce preprocessing directives for calls
!                             of subroutines with and without PFIO components.
!                             A new code that contains the PFIO statements spawn 
!                             all the MPI processes, creates the compute nodes 
!                             and IO nodes (based on the command line parameters) 
!                             and to drive the model.
!
program lisdrv
!
! !USES:       
      use LIS_coreMod,    only : LIS_config_init, LIS_rc, LIS_masterproc
      use LIS_ftimingMod, only : Ftiming_init, Ftiming_Output
      use LIS_ftimingMod, only : Ftiming_On, Ftiming_Off
      use LIS_logMod

#ifdef ESMF_TRACE
      use ESMF
#endif


#ifdef USE_PFIO
      use mpi
      use ESMF
      use MAPL
      use FLAP
      use pFIO_ClientManagerMod, only: o_Clients
      !use pFIO_ClientManagerMod, only: i_Clients

      implicit none

      !type(MAPL_FargparseCLI) :: cmd_line
      type(MAPL_FlapCLI)      :: cmd_line
      type(MAPL_CapOptions)   :: cap_options
      type(ServerManager)     :: ioserver_manager
      type(SplitCommunicator) :: split_comm

      integer :: client_comm,rank, npes, ierror, provided, required
      integer :: npes_world, npes_model
      integer :: pe_id
      integer :: status, rc
      integer :: subcommunicator

      integer :: ftim_lun, lun_array(1)
      integer :: root_cpu
!EOP
!------------------------------------------------------------------------------
!BOC  
      ! Read and parse the command line, and set parameters
      ! If you have extra options you make a procedure as seen below and add arguments
      ! there and pass in here
      cmd_line = MAPL_FlapCLI(description = 'LIS', &
                              authors     = 'LIS', &
                              extra       = extra_options)

      cap_options = MAPL_CapOptions(cmd_line)

      ! MPI is not initialized yet. Only Check if it is.
      ! MPI will be initialized under "model"
      call initialize_mpi(MPI_COMM_WORLD)

      call MPI_Comm_size(MPI_COMM_WORLD, npes_world, ierror)
      if (cap_options%npes_model == -1) then
          cap_options%npes_model = npes_world
      endif

      ! Initialize the IO Server Manager using parameters defined above
      subcommunicator = create_member_subcommunicator(MPI_COMM_WORLD, 1, npes_world)

      IF (subcommunicator /= MPI_COMM_NULL) THEN
         call initialize_ioserver(subcommunicator)

         call ioserver_manager%get_splitcomm(split_comm)

         SELECT CASE(split_comm%get_name())
         CASE('model')
            ! Get the model MPI communicator
            client_comm = split_comm%get_subcommunicator()

            ! Get the PE id
            call MPi_Comm_rank(client_comm, pe_id, ierror)
            if (pe_id == 0) PRINT*,"Start running model"

            ! Get the number of PEs used for the model
            call MPi_Comm_size(client_comm, npes_model, ierror)
            IF (npes_model /= cap_options%npes_model) STOP "sanity check failed"

            if (pe_id == 0) PRINT*,"Start running LIS_config_init"
            !call LIS_config_init(comm=client_comm)
            call LIS_config_init(cmd_args=.true., comm=client_comm)

            if (pe_id == 0) PRINT*,"Start running lisinit"
            TRACE_ENTER("LIS_init")

            if (LIS_rc%do_ftiming) then
               call Ftiming_Init ( )
               call Ftiming_On("Total Model")
               call Ftiming_On("LIS_init")
            endif

            call lisinit(trim(LIS_rc%runmode)//char(0))

            if (LIS_rc%do_ftiming) call Ftiming_Off("LIS_init")

            TRACE_EXIT("LIS_init")

            ! if there are multiple oserver, split it into large and small pool
            call o_Clients%split_server_pools()

            if (pe_id == 0) PRINT*,"Start running lisrun"
            TRACE_ENTER("LIS_run")

            if (LIS_rc%do_ftiming) call Ftiming_On("LIS_run")
            call lisrun(trim(LIS_rc%runmode)//char(0))
            if (LIS_rc%do_ftiming) then
               call Ftiming_Off("LIS_run")
               call Ftiming_Off("Total Model")
            endif

            TRACE_EXIT("LIS_run")
            if (pe_id == 0) PRINT*,"Done with lisrun"

            if (LIS_rc%do_ftiming) then
               ftim_lun = 499
               if ( LIS_masterproc ) then
                  ftim_lun = LIS_getNextUnitNumber()
               endif

               call Ftiming_Output(ftim_lun)

               if ( LIS_masterproc ) then
                  call LIS_releaseUnitNumber(ftim_lun)
               endif
            endif

            if (pe_id == 0) PRINT*,"Start running lisfinalize"
            call lisfinalize(trim(LIS_rc%runmode)//char(0))

            !call i_Clients%terminate()
            call o_Clients%terminate()

            call ESMF_Finalize(endflag=ESMF_END_KEEPMPI, rc=status)

         END SELECT
      END IF

      call ioserver_manager%finalize()

      call MPI_finalize(ierror)

!-------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
      subroutine extra_options(options, rc)

         type (command_line_interface), intent(inout) :: options
         integer, intent(out), optional :: rc

         call options%add( switch='--file', &
                           switch_ab='-f', &
                           help='LIS configuration file', &
                           def='lis.config', &
                           act='store', &
                           error=status)
         !_VERIFY(status)

         !_RETURN(_SUCCESS)
         if (present(rc)) rc = 0
      end subroutine extra_options
!      subroutine extra_options(parser, rc)
!         type (ArgParser), intent(inout) :: parser
!         integer, intent(out), optional :: rc
!
!         call parser%add_argument('-f', '--file', &
!              help='LIS configuration file', &
!              type='string', &
!              default='lis.config', &
!              action='store')
!
!         !_RETURN(_SUCCESS)
!         if (present(rc)) rc = 0
!
!      end subroutine extra_options
!------------------------------------------------------------------------------
   integer function create_member_subcommunicator(comm, n_members, npes_member, rc) result(subcommunicator)
      use MAPL_SimpleCommSplitterMod
      integer, intent(in) :: comm
      integer, intent(in) :: n_members, npes_member
      integer, optional, intent(out) :: rc

      type(SimpleCommSplitter) :: splitter
      type (SplitCommunicator) :: split_comm

      subcommunicator = MPI_COMM_NULL ! in case of failure
      splitter = SimpleCommSplitter(comm, n_members, npes_member)
      split_comm = splitter%split(rc=status)
      !_VERIFY(status)
      subcommunicator = split_comm%get_subcommunicator()

   end function create_member_subcommunicator
!
!------------------------------------------------------------------------------
!
     ! Check if MPI is already initialized. 
     ! If MPI is not initialized, then we initialize the MPI execution 
     ! environment with a single thread.
     subroutine initialize_mpi(comm)
         integer, intent(in) :: comm
         logical :: mpi_already_initialized
         integer :: ierror
         integer :: provided

         call MPI_Initialized(mpi_already_initialized, ierror)
         if (.not. mpi_already_initialized) then
            call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierror)
            _VERIFY(ierror)
            _ASSERT(provided == MPI_THREAD_SINGLE, "MPI_THREAD_SINGLE not supported by this MPI.")
         end if

         !call MPI_Comm_rank(comm, pe_rank, ierror); _VERIFY(ierror)
         !call MPI_Comm_size(comm, npes_world, ierror); _VERIFY(ierror)
      end subroutine initialize_mpi
!
!------------------------------------------------------------------------------
!
      ! Initialize the IO server and capture the command line configuration
      subroutine initialize_ioserver(comm)
         integer, intent(in) :: comm
         call ioserver_manager%initialize(comm, &
                       application_size     = cap_options%npes_model, &
                       nodes_input_server   = cap_options%nodes_input_server, &
                       nodes_output_server  = cap_options%nodes_output_server, &
                       npes_input_server    = cap_options%npes_input_server, &
                       npes_output_server   = cap_options%npes_output_server, &
                       oserver_type         = cap_options%oserver_type, &
                       npes_backend_pernode = cap_options%npes_backend_pernode, &
                       isolate_nodes        = cap_options%isolate_nodes, &
                       fast_oclient         = cap_options%fast_oclient, &
                       with_profiler        = cap_options%with_io_profiler, &
                    rc=status)
         !_VERIFY(status)
      end subroutine initialize_ioserver
!
!------------------------------------------------------------------------------
!
#else
#if ( defined SPMD )
      use LIS_mpiMod,     only : LIS_mpi_comm
      use mpi
#endif
!
!EOP
     implicit none

      integer :: ftim_lun, lun_array(1)
      integer :: root_cpu, ierr
!BOC
      call LIS_config_init(cmd_args=.true.)

      TRACE_ENTER("LIS_init")

      if (LIS_rc%do_ftiming) then
         call Ftiming_Init ( )
         call Ftiming_On("Total Model")
         call Ftiming_On("LIS_init")
      endif
      if ( LIS_masterproc ) PRINT*,"Start running lisinit"
      call lisinit(trim(LIS_rc%runmode)//char(0))
      if ( LIS_masterproc ) PRINT*,"Done with lisinit"
      if (LIS_rc%do_ftiming) call Ftiming_Off("LIS_init")

      TRACE_EXIT("LIS_init")

      TRACE_ENTER("LIS_run")

      if (LIS_rc%do_ftiming) call Ftiming_On("LIS_run")
      if ( LIS_masterproc ) PRINT*,"Start running lisrun"
      call lisrun(trim(LIS_rc%runmode)//char(0))
      if ( LIS_masterproc ) PRINT*,"Done with lisrun"
      if (LIS_rc%do_ftiming) then
         call Ftiming_Off("LIS_run")
         call Ftiming_Off("Total Model")
      endif

      TRACE_EXIT("LIS_run")

      if (LIS_rc%do_ftiming) then
         if ( LIS_masterproc ) PRINT*,"Collecting the timing statistics"
         ftim_lun = -999
         if ( LIS_masterproc ) then
            ftim_lun = LIS_getNextUnitNumber()
         endif

         root_cpu = 0

         lun_array(1) = ftim_lun
#if ( defined SPMD )
         call MPI_Bcast (lun_array, 1, MPI_INTEGER, root_cpu, LIS_mpi_comm, ierr)
#endif
         ftim_lun     = lun_array(1)

         call Ftiming_Output(ftim_lun)
         if ( LIS_masterproc ) PRINT*,"Done with the timing statistics"
      endif

      if ( LIS_masterproc ) PRINT*,"Start running lisfinalize"
      call lisfinalize(trim(LIS_rc%runmode)//char(0))

  !call LIS_config_init(cmd_args=.true.)
  !TRACE_ENTER("LIS_init")
  !call lisinit(trim(LIS_rc%runmode)//char(0))
  !TRACE_EXIT("LIS_init")
  !TRACE_ENTER("LIS_run")
  !call lisrun(trim(LIS_rc%runmode)//char(0))
  !TRACE_EXIT("LIS_run")
  !call lisfinalize(trim(LIS_rc%runmode)//char(0))
!EOC
#endif
end program lisdrv


