!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

MODULE clm2_shr_msg_mod

#if (! defined HIDE_SHR_MSG)
#if (! defined HIDE_MPI)

   use clm2_shr_kind_mod, only: SHR_KIND_IN, SHR_KIND_R8       ! defines real & integer kinds
   use clm2_shr_sys_mod, only: shr_sys_getenv, shr_sys_flush, shr_sys_abort, shr_sys_chdir, LIS_mpi_comm  ! system calls

   IMPLICIT none
   
   PRIVATE
!
! PUBLIC: Interfaces and data
!
   public clm2_shr_msg_stdio, shr_msg_init, shr_msg_groups, shr_msg_errstr, shr_msg_finalize, &
          clm2_shr_msg_send_i, shr_msg_send_r, shr_msg_send_c, shr_msg_recv_i, shr_msg_recv_r, &
          clm2_shr_msg_recv_c, shr_msg_wait, shr_msg_test
   !-----------------------------------------------------------------
   ! tid's for inter-model msg passing in LIS_mpi_comm
   !-----------------------------------------------------------------
   integer(SHR_KIND_IN), public :: SHR_MSG_TID_CPL   ! cpl model task id
   integer(SHR_KIND_IN), public :: SHR_MSG_TID_ATM   ! atm model task id 
   integer(SHR_KIND_IN), public :: SHR_MSG_TID_OCN   ! ocn model task id 
   integer(SHR_KIND_IN), public :: SHR_MSG_TID_ICE   ! ice model task id 
   integer(SHR_KIND_IN), public :: SHR_MSG_TID_LND   ! lnd model task id 

   !-----------------------------------------------------------------
   ! intra-model communicators
   !-----------------------------------------------------------------
   integer(SHR_KIND_IN), public :: SHR_MSG_COMM_ATM  ! atm model communicator
   integer(SHR_KIND_IN), public :: SHR_MSG_COMM_CPL  ! cpl model communicator
   integer(SHR_KIND_IN), public :: SHR_MSG_COMM_ICE  ! ice model communicator
   integer(SHR_KIND_IN), public :: SHR_MSG_COMM_LND  ! lnd model communicator
   integer(SHR_KIND_IN), public :: SHR_MSG_COMM_OCN  ! ocn model communicator

   !-----------------------------------------------------------------
   ! inter-model message tags
   !-----------------------------------------------------------------
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2AI = 11 ! c ->a, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2A  = 10 ! c ->a 
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_A2CI = 20 ! c<- a, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_A2C  = 21 ! c<- a

   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2II = 31 ! c ->i, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2I  = 30 ! c ->i
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_I2CI = 40 ! c<- i, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_I2C  = 41 ! c<- i

   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2OI = 51 ! c ->o, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2O  = 50 ! c ->o
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_O2CI = 60 ! c<- o, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_O2C  = 61 ! c<- o

   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2LI = 71 ! c ->l, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_C2L  = 70 ! c ->l
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_L2CI = 80 ! c<- l, initial
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_TAG_L2C  = 81 ! c<- l

   !-----------------------------------------------------------------
   ! interface compatiblility ID numbers for cpl/atm interface
   !-----------------------------------------------------------------
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MAJ_UNDEF =     0
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MAJ_V00   =  1000
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MAJ_V01   =  1001
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MAJ_V02   =  1002
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MAJ_V03   =  1003
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MAJ_V04   =  1004
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MAJ_V05   =  1005

   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_UNDEF =     0
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V00   = 10000
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V01   = 10001
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V02   = 10002
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V03   = 10003
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V04   = 10004
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V05   = 10005
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V06   = 10006
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V07   = 10007
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V08   = 10008
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_A_MIN_V09   = 10009

   !-----------------------------------------------------------------
   ! interface compatiblility ID numbers for cpl/ice interface
   !-----------------------------------------------------------------
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MAJ_UNDEF =     0
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MAJ_V00   =  2000
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MAJ_V01   =  2001
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MAJ_V02   =  2002
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MAJ_V03   =  2003
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MAJ_V04   =  2004
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MAJ_V05   =  2005

   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_UNDEF =     0
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V00   = 20000
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V01   = 20001
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V02   = 20002
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V03   = 20003
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V04   = 20004
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V05   = 20005
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V06   = 20006  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V07   = 20007  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V08   = 20008  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_I_MIN_V09   = 20009  

   !-----------------------------------------------------------------
   ! interface compatiblility ID numbers for cpl/lnd interface
   !-----------------------------------------------------------------
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MAJ_UNDEF =     0
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MAJ_V00   =  3002
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MAJ_V01   =  3001
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MAJ_V02   =  3002
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MAJ_V03   =  3003
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MAJ_V04   =  3004
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MAJ_V05   =  3005

   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_UNDEF =     0
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V00   = 30000
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V01   = 30001
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V02   = 30002
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V03   = 30003
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V04   = 30004
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V05   = 30005
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V06   = 30006
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V07   = 30007
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V08   = 30008
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_L_MIN_V09   = 30009

   !-----------------------------------------------------------------
   ! interface compatiblility ID numbers for cpl/ocn interface
   !-----------------------------------------------------------------
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MAJ_UNDEF =     0  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MAJ_V01   =  4001  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MAJ_V02   =  4002  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MAJ_V03   =  4003  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MAJ_V04   =  4004  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MAJ_V05   =  4005  

   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_UNDEF =     0  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V00   = 40000  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V01   = 40001  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V02   = 40002  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V03   = 40003  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V04   = 40004  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V05   = 40005  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V06   = 40006  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V07   = 40007  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V08   = 40008  
   integer(SHR_KIND_IN), public, parameter :: SHR_MSG_O_MIN_V09   = 40009  
! mpi library include file, in clm2_shr_sys_mod
#include <mpif.h>


CONTAINS

!===============================================================================

SUBROUTINE clm2_shr_msg_stdio(model)

   !--- arguments ---
   character(len=*),intent(in) :: model ! used to construct env varible name

   !--- local ---
   character(len=  8)   :: var_dir   ! env variable name wrt cwd
   character(len=  8)   :: var_in    ! env variable name wrt stdin
   character(len=  8)   :: var_out   ! env variable name wrt stdout
   character(len=256)   :: str_dir   ! env variable value for cwd
   character(len=256)   :: str_in    ! env variable value for stdin file
   character(len=256)   :: str_out   ! env variable value for stdout file
   integer(SHR_KIND_IN) :: rcode_dir ! 0 => no error on cwd    system call
   integer(SHR_KIND_IN) :: rcode_in  ! 0 => no error on stdin  system call
   integer(SHR_KIND_IN) :: rcode_out ! 0 => no error on stdout system call

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_stdio) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   1) change the cwd (current working directory) and 
!   2) redirect stdin & stdout (units 5 & 6) to named files,
!   where the desired cwd & files are specified by env variables.
!
!   Normally this is done to work around limitations in the execution syntax of
!   common MPI implementations.  For example, SGI's mpirun syntax is not 
!   flexible enough (?) to allow MPMD models to select different execution
!   directories or to redirect stdin & stdout on the command line.  Such 
!   functionality is highly desireable for the CSM purposes.  
!   ie. mpirun can't handle this:
!   unix> cd /usr/tmp/jdoe/csm/case01/atm ; atm < atm.parm > atm.log &
!   unix> cd /usr/tmp/jdoe/csm/case01/cpl ; cpl < cpl.parm > cpl.log &
!   etc.
!
! ASSUMPTIONS:
! o the env variables are:
!   <model>_dir - the desired cwd, eg. setenv atm_dir /usr/tmp/atm
!   <model>_in  - the stdin  file, eg. setenv atm_in  /usr/tmp/atm.parm
!   <model>_out - the stdout file, eg. setenv atm_out /usr/tmp/atm.log
!   where <model> is a character string, typically one of 
!   "cpl", "atm", "lnd", "ice", or "ocn".
!
!-------------------------------------------------------------------------------
 
   !--- Construct env variable names & (attempt to) get their values ---
   str_dir = '  '
   str_in  = '  '
   str_out = '  '

   var_dir = trim(model) // "_dir "
   var_in  = trim(model) // "_in  "
   var_out = trim(model) // "_out "

   call clm2_shr_sys_getenv(var_dir,str_dir,rcode_dir)
   if (rcode_dir == 0 .and. len_trim(str_dir) == 0) rcode_dir = 999

   call clm2_shr_sys_getenv(var_in ,str_in ,rcode_in )
   if (rcode_in  == 0 .and. len_trim(str_in ) == 0) rcode_in  = 999

   call clm2_shr_sys_getenv(var_out,str_out,rcode_out)
   if (rcode_out == 0 .and. len_trim(str_out) == 0) rcode_out = 999

   !--- 1st change cwd, then open units 5 & 6 (if getenv was succesful) ---
   if (rcode_dir == 0) call clm2_shr_sys_chdir(str_dir,rcode_dir)
   if (rcode_in  == 0) open(unit=5,file=str_in ,status='UNKNOWN')
   if (rcode_out == 0) open(unit=6,file=str_out,position='APPEND')
 
   !--- print informational messages (re-set unit 6 before doing this) ---
   write(6,F00) 'chdir & open units 5 & 6 for model = ',model

   if (rcode_dir == 0) then
     write(6,F00) 'changed cwd to      ',trim(str_dir)
   else
     write(6,F00) 'cwd *not* changed'
   endif
   if (rcode_in == 0) then
     write(6,F00) 'unit 5 connected to ',trim(str_in)
   else
     write(6,F00) 'unit 5 has *not* been redirected'
   endif
   if (rcode_out == 0) then
     write(6,F00) 'unit 6 connected to ',trim(str_out)
   else
     write(6,F00) 'unit 6 has *not* been redirected'
   endif

   call clm2_shr_sys_flush(6)
 
END SUBROUTINE clm2_shr_msg_stdio

!===============================================================================

SUBROUTINE clm2_shr_msg_init(model)

   !--- arguments ---
   character(len=3),intent(in) :: model ! used to construct env varible name

   !----- local -----
   integer(SHR_KIND_IN) :: rcode

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_init) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: initialize MPI
!-------------------------------------------------------------------------------

   call mpi_init(rcode)
   if (rcode /= MPI_SUCCESS) then
      write(6,F00) 'MPI ERROR:  error code: ',rcode
      write(6,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
   write(6,F00) 'connected to mpi fabric'

   call clm2_shr_sys_flush(6)

END SUBROUTINE clm2_shr_msg_init

!===============================================================================

SUBROUTINE clm2_shr_msg_groups (model)

   !--- arguments ---
   character(len=3), intent(in) :: model  ! model name: atm,lnd,ocn,ice,cpl

   !--- local ---
   character(len=3)             :: cmodel ! model name, temporary

   integer(SHR_KIND_IN), dimension(3) :: range ! array for creating groups

   integer(SHR_KIND_IN) :: group_world   ! group id for LIS_mpi_comm
   integer(SHR_KIND_IN) :: group_cpl     ! group of processors assigned to cpl
   integer(SHR_KIND_IN) :: group_atm     ! group of processors assigned to atm
   integer(SHR_KIND_IN) :: group_ocn     ! group of processors assigned to ocn
   integer(SHR_KIND_IN) :: group_ice     ! group of processors assigned to ice
   integer(SHR_KIND_IN) :: group_lnd     ! group of processors assigned to lnd

   integer(SHR_KIND_IN) :: n             ! dummy loop counter
   integer(SHR_KIND_IN) :: rcode         ! return code for MPI calls
   integer(SHR_KIND_IN) :: nprocs_world  ! total processor (task) count
   integer(SHR_KIND_IN) :: my_rank_world ! my rank (tid) wrt to LIS_mpi_comm
   integer(SHR_KIND_IN) :: atm_rank_min  ! processor range for each component
   integer(SHR_KIND_IN) :: atm_rank_max
   integer(SHR_KIND_IN) :: cpl_rank_min
   integer(SHR_KIND_IN) :: cpl_rank_max
   integer(SHR_KIND_IN) :: ice_rank_min
   integer(SHR_KIND_IN) :: ice_rank_max
   integer(SHR_KIND_IN) :: lnd_rank_min
   integer(SHR_KIND_IN) :: lnd_rank_max
   integer(SHR_KIND_IN) :: ocn_rank_min
   integer(SHR_KIND_IN) :: ocn_rank_max
   character(len=MPI_MAX_ERROR_STRING) :: str ! mpi error string

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_groups) ',4a)"
   character(len=*),parameter :: F01 = "('(clm2_shr_msg_groups) ',a,i4,a,i4)"


!-------------------------------------------------------------------------------
! PURPOSE:
!   this routine queries all the components of the full coupled system and sets
!   up proper communicators groups and task ids for each component model
!
!   this routine should be called after mpi_init, but before setting up any 
!   internal mpi setups (since these will require the internal communicators 
!   returned by this routine)
!
! ASSUMPTION:
!   the ranks of all tasks associated with any one model are contiguous (!)
!-------------------------------------------------------------------------------

   !-------------------------------------------------------------
   ! determine processor rank wrt LIS_mpi_comm
   !-------------------------------------------------------------
   call mpi_comm_rank  (LIS_mpi_comm, my_rank_world, rcode)
   write(6,F01) 'establish intra-model communicators, '// &
   'my LIS_mpi_comm rank is ',my_rank_world

   !-------------------------------------------------------------
   ! determine which group of processes assigned to each model
   !-------------------------------------------------------------
   call mpi_comm_size (LIS_mpi_comm, nprocs_world, rcode)

   atm_rank_min = nprocs_world
   atm_rank_max = 0
   ocn_rank_min = nprocs_world
   ocn_rank_max = 0
   ice_rank_min = nprocs_world
   ice_rank_max = 0
   lnd_rank_min = nprocs_world
   lnd_rank_max = 0
   cpl_rank_min = nprocs_world
   cpl_rank_max = 0

   !-------------------------------------------------------------
   ! each processor broadcasts its model to all other processors 
   !-------------------------------------------------------------
   do n=0,nprocs_world-1
     if (n == my_rank_world) then
       cmodel = model
     else
       cmodel = 'unk'
     endif

     call mpi_bcast(cmodel, 3, MPI_CHARACTER, n, LIS_mpi_comm, rcode)
     call clm2_shr_msg_errstr(rcode,'mpi_bcast(cmodel)')

     select case(cmodel)
     case ('atm')
       atm_rank_min = min(atm_rank_min, n)
       atm_rank_max = max(atm_rank_max, n)
     case ('cpl')
       cpl_rank_min = min(cpl_rank_min, n)
       cpl_rank_max = max(cpl_rank_max, n)
     case ('ice')
       ice_rank_min = min(ice_rank_min, n)
       ice_rank_max = max(ice_rank_max, n)
     case ('lnd')
       lnd_rank_min = min(lnd_rank_min, n)
       lnd_rank_max = max(lnd_rank_max, n)
     case ('ocn')
       ocn_rank_min = min(ocn_rank_min, n)
       ocn_rank_max = max(ocn_rank_max, n)
     case ('dud')
     case default
       write(6,F00) 'Unknown model ',cmodel,' in world communicator group'
       write(6,F00) 'Model must be atm, cpl, ice, lnd, ocn, dud'
       call clm2_shr_sys_abort()
     end select

   end do

   !----------------------------------------------------------------------
   ! for each component, assume the first task in the world comm group
   ! is the task that will communicate coupler/model messages
   !----------------------------------------------------------------------

   SHR_MSG_TID_ATM = atm_rank_min
   SHR_MSG_TID_CPL = cpl_rank_min
   SHR_MSG_TID_ICE = ice_rank_min
   SHR_MSG_TID_LND = lnd_rank_min
   SHR_MSG_TID_OCN = ocn_rank_min

   !----------------------------------------------------------------------
   ! create subroup and communicators for each models internal 
   ! communciations, note that MPI_COMM_CREATE must be called by all 
   ! processes in LIS_mpi_comm so this must be done by all models
   ! consistently and in the same order.
   !----------------------------------------------------------------------

   call mpi_comm_group(LIS_mpi_comm, group_world, rcode)
   call clm2_shr_msg_errstr(rcode,'mpi_comm_group (LIS_mpi_comm,group_world)')

   range(3) = 1

   range(1) = atm_rank_min
   range(2) = atm_rank_max
   call mpi_group_range_incl(group_world, 1, range, group_atm, rcode)

   range(1) = ocn_rank_min
   range(2) = ocn_rank_max
   call mpi_group_range_incl(group_world, 1, range, group_ocn, rcode)

   range(1) = ice_rank_min
   range(2) = ice_rank_max
   call mpi_group_range_incl(group_world, 1, range, group_ice, rcode)

   range(1) = lnd_rank_min
   range(2) = lnd_rank_max
   call mpi_group_range_incl(group_world, 1, range, group_lnd, rcode)

   range(1) = cpl_rank_min
   range(2) = cpl_rank_max
   call mpi_group_range_incl(group_world, 1, range, group_cpl, rcode)

   call mpi_comm_create(LIS_mpi_comm, group_atm, SHR_MSG_COMM_ATM, rcode)
   call mpi_comm_create(LIS_mpi_comm, group_ocn, SHR_MSG_COMM_OCN, rcode)
   call mpi_comm_create(LIS_mpi_comm, group_ice, SHR_MSG_COMM_ICE, rcode)
   call mpi_comm_create(LIS_mpi_comm, group_lnd, SHR_MSG_COMM_LND, rcode)
   call mpi_comm_create(LIS_mpi_comm, group_cpl, SHR_MSG_COMM_CPL, rcode)
 
   !----------------------------------------------------------------------
   ! record inter-model TID's and intra-model communicator groups
   !----------------------------------------------------------------------
   if (  my_rank_world == SHR_MSG_TID_ATM  &
   .or.  my_rank_world == SHR_MSG_TID_CPL  &
   .or.  my_rank_world == SHR_MSG_TID_ICE  &
   .or.  my_rank_world == SHR_MSG_TID_LND  &
   .or.  my_rank_world == SHR_MSG_TID_OCN  ) then

       write(6,F01) 'SHR_MSG_TID_ATM=',SHR_MSG_TID_ATM
       write(6,F01) 'SHR_MSG_TID_CPL=',SHR_MSG_TID_CPL
       write(6,F01) 'SHR_MSG_TID_ICE=',SHR_MSG_TID_ICE
       write(6,F01) 'SHR_MSG_TID_LND=',SHR_MSG_TID_LND
       write(6,F01) 'SHR_MSG_TID_OCN=',SHR_MSG_TID_OCN

       write(6,F01) 'group_atm=',atm_rank_min,'...',atm_rank_max
       write(6,F01) 'group_cpl=',cpl_rank_min,'...',cpl_rank_max
       write(6,F01) 'group_ice=',ice_rank_min,'...',ice_rank_max
       write(6,F01) 'group_lnd=',lnd_rank_min,'...',lnd_rank_max
       write(6,F01) 'group_ocn=',ocn_rank_min,'...',ocn_rank_max

   end if

   if ( model == 'atm') then
      call mpi_comm_rank(SHR_MSG_COMM_ATM, n, rcode)
      write(6,F01) &
      & 'LIS_mpi_comm rank',my_rank_world,'<=> SHR_MSG_COMM_ATM rank',n
   else if ( model == 'cpl')  then
      call mpi_comm_rank(SHR_MSG_COMM_CPL, n, rcode)
      write(6,F01) &
      & 'LIS_mpi_comm rank',my_rank_world,'<=> SHR_MSG_COMM_CPL rank',n
   else if ( model == 'ice')  then
      call mpi_comm_rank(SHR_MSG_COMM_ICE, n, rcode)
      write(6,F01) &
      & 'LIS_mpi_comm rank',my_rank_world,'<=> SHR_MSG_COMM_ICE rank',n
   else if ( model == 'lnd')  then
      call mpi_comm_rank(SHR_MSG_COMM_LND, n, rcode)
      write(6,F01) &
      & 'LIS_mpi_comm rank',my_rank_world,'<=> SHR_MSG_COMM_LND rank',n
   else if ( model == 'ocn') then
      call mpi_comm_rank(SHR_MSG_COMM_OCN, n, rcode)
      write(6,F01) &
      & 'LIS_mpi_comm rank',my_rank_world,'<=> SHR_MSG_COMM_OCN rank',n
   else if ( model == 'dud') then
      write(6,F01) &
      & 'LIS_mpi_comm rank',my_rank_world,'<=> not in any group rank',-1
   else
      write(6,F00) 'invalid model =',model
   end if

   call clm2_shr_sys_flush(6)

END SUBROUTINE clm2_shr_msg_groups

!===============================================================================

SUBROUTINE clm2_shr_msg_errstr(rcode,instr)

   !----- arguments -----
   integer(SHR_KIND_IN)                :: rcode   ! mpi error code
   character(len=*)                    :: instr   ! input string
   character(len=MPI_MAX_ERROR_STRING) :: errstr  ! mpi error string

   !----- local -----
   integer(SHR_KIND_IN)                :: n       ! length of mpi string
   integer(SHR_KIND_IN)                :: rank    ! mpi rank
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_errstr) ',a,i3,1x,4a)"

!-------------------------------------------------------------------------------
! PURPOSE: if there is an error, print the corresponding error message
!-------------------------------------------------------------------------------

   if (rcode /= MPI_SUCCESS) then
     call mpi_comm_rank(LIS_mpi_comm, rank, rcode)
     call mpi_error_string(rcode,errstr,n,rcode)
     write(*,F00) 'wrld rank=',rank,trim(instr),' => MPI error: ',trim(errstr)
   end if

END SUBROUTINE clm2_shr_msg_errstr

!===============================================================================

SUBROUTINE clm2_shr_msg_finalize

   !----- local -----
   integer(SHR_KIND_IN)                :: rcode   ! mpi error code

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_finalize) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: finalize MPI
!-------------------------------------------------------------------------------

   call mpi_finalize(rcode)
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR:  error code: ',rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
   write(*,F00) 'disconnecting from mpi'

END SUBROUTINE clm2_shr_msg_finalize

!===============================================================================

SUBROUTINE clm2_shr_msg_send_i(arr, n, tid, msgid, async_in, reqid)
 
   !----- arguments -----
   integer(SHR_KIND_IN) :: n        ! array size
   integer(SHR_KIND_IN) :: arr(n)   ! array
   integer(SHR_KIND_IN) :: tid      ! destination tid
   integer(SHR_KIND_IN) :: msgid    ! message id
   logical, intent(in),               optional :: async_in ! async flag
   integer(SHR_KIND_IN), intent(out), optional :: reqid    ! request id

   !----- local -----
   integer(SHR_KIND_IN)      :: rcode        ! return code
   logical                   :: async        ! flag for async mpi

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_send_i) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: send the integer :: array "arr(n)" to task "tid" with id "msgid".
!-------------------------------------------------------------------------------

   async   = .false. ; if (PRESENT(async_in)) async   = async_in

!  write(*,F00) 'sending int  msg ', msgid,' to model ',tid,' size=',n,' async=',async
   if (async) then
     call mpi_isend(arr,n,MPI_INTEGER,tid,msgid,LIS_mpi_comm,reqid,rcode)
   else
     call mpi_send (arr,n,MPI_INTEGER,tid,msgid,LIS_mpi_comm,      rcode)
   endif
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR: id=',msgid,'  error code: ',rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_send_i

!===============================================================================

SUBROUTINE clm2_shr_msg_send_r(arr, n, tid, msgid, async_in, reqid)
   
   !----- arguments -----
   integer(SHR_KIND_IN) :: n        ! array size
   real   (SHR_KIND_R8) :: arr(n)   ! array
   integer(SHR_KIND_IN) :: tid      ! destination tid
   integer(SHR_KIND_IN) :: msgid    ! message id
   logical, intent(in),               optional :: async_in ! async flag
   integer(SHR_KIND_IN), intent(out), optional :: reqid    ! request id

   !----- local -----
   integer(SHR_KIND_IN)      :: rcode        ! return code
   logical                   :: async        ! flag for async mpi

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_send_r) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: send the real array "arr(n)" to task "tid" with id "msgid".
!-------------------------------------------------------------------------------

   async   = .false. ; if (PRESENT(async_in)) async   = async_in

!  write(*,F00) 'sending real :: msg ', msgid,' to model ',tid,' size=',n,' async=',async
   if (async) then
     call mpi_isend(arr,n,MPI_REAL8  ,tid,msgid,LIS_mpi_comm,reqid,rcode)
   else
     call mpi_send (arr,n,MPI_REAL8  ,tid,msgid,LIS_mpi_comm,      rcode)
   endif
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR: id=',msgid,'  error code: ',rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_send_r

!===============================================================================

SUBROUTINE clm2_shr_msg_send_c(arr, n, tid, msgid, async_in, reqid)

   !----- arguments -----
   integer(SHR_KIND_IN) :: n        ! array size
   character            :: arr(n)   ! array
   integer(SHR_KIND_IN) :: tid      ! destination tid
   integer(SHR_KIND_IN) :: msgid    ! message id
   logical, intent(in),               optional :: async_in ! async flag
   integer(SHR_KIND_IN), intent(out), optional :: reqid    ! request id

   !----- local -----
   integer(SHR_KIND_IN)      :: rcode        ! return code
   logical                   :: async        ! flag for async mpi

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_send_c) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: Send the character array "arr(n)" to task "tid" with id "msgid".
!-------------------------------------------------------------------------------

   async   = .false. ; if (PRESENT(async_in)) async   = async_in

!  write(*,F00) 'sending char msg ', msgid,' to model ',tid,' size=',n,' async=',async
   if (async) then
     call mpi_isend(arr,n,MPI_CHARACTER,tid,msgid,LIS_mpi_comm,reqid,rcode)
   else
     call mpi_send (arr,n,MPI_CHARACTER,tid,msgid,LIS_mpi_comm,      rcode)
   endif
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR: id=',msgid,'  error code: ',rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_send_c

!===============================================================================

SUBROUTINE clm2_shr_msg_recv_i(arr, n, tid, msgid, async_in, reqid)
         
   !----- arguments -----
   integer(SHR_KIND_IN) :: n        ! array size
   integer(SHR_KIND_IN) :: arr(n)   ! array
   integer(SHR_KIND_IN) :: tid      ! destination tid
   integer(SHR_KIND_IN) :: msgid    ! message id
   logical, intent(in),               optional :: async_in ! async flag
   integer(SHR_KIND_IN), intent(out), optional :: reqid    ! request id

   !----- local -----
   integer(SHR_KIND_IN) :: rcode                    ! return code
   integer(SHR_KIND_IN) :: status(MPI_STATUS_SIZE)  ! request status
   logical                   :: async        ! flag for async mpi

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_recv_i) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: receive the integer array "arr(n)" from task "tid" with id "msgid".
!-------------------------------------------------------------------------------

   async   = .false. ; if (PRESENT(async_in)) async   = async_in

!  write(*,F00) 'recving int  msg ', msgid,' fm model ',tid,' size=',n,' async=',async
   if (async) then
     call mpi_irecv(arr, n, MPI_INTEGER, tid, msgid, LIS_mpi_comm, reqid,  rcode)
   else
     call mpi_recv (arr, n, MPI_INTEGER, tid, msgid, LIS_mpi_comm, status, rcode)
   endif
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR: id=',msgid,'  error code: ',rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_recv_i

!===============================================================================

SUBROUTINE clm2_shr_msg_recv_r(arr, n, tid, msgid, async_in, reqid)
          
   !----- arguments -----
   integer(SHR_KIND_IN) :: n        ! array size
   real   (SHR_KIND_R8) :: arr(n)   ! array
   integer(SHR_KIND_IN) :: tid      ! destination tid
   integer(SHR_KIND_IN) :: msgid    ! message id
   logical, intent(in),               optional :: async_in ! async flag
   integer(SHR_KIND_IN), intent(out), optional :: reqid    ! request id

   !----- local -----
   integer(SHR_KIND_IN) :: rcode                    ! return code
   integer(SHR_KIND_IN) :: status(MPI_STATUS_SIZE)  ! request status
   logical                   :: async        ! flag for async mpi

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_recv_r) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: receive the real array "arr(n)" from task "tid" with id "msgid".
!-------------------------------------------------------------------------------

   async   = .false. ; if (PRESENT(async_in)) async   = async_in

!  write(*,F00) 'recving real :: msg ', msgid,' fm model ',tid,' size=',n,' async=',async
   if (async) then
     call mpi_irecv(arr, n, MPI_REAL8, tid, msgid, LIS_mpi_comm, reqid,  rcode)
   else
     call mpi_recv (arr, n, MPI_REAL8, tid, msgid, LIS_mpi_comm, status, rcode)
   endif
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR: id=',msgid,'  error code: ',rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_recv_r

!===============================================================================
         
SUBROUTINE clm2_shr_msg_recv_c(arr, n, tid, msgid, async_in, reqid)

   !----- arguments -----
   integer(SHR_KIND_IN) :: n        ! array size
   character            :: arr(n)   ! array
   integer(SHR_KIND_IN) :: tid      ! destination tid
   integer(SHR_KIND_IN) :: msgid    ! message id
   logical, intent(in),               optional :: async_in ! async flag
   integer(SHR_KIND_IN), intent(out), optional :: reqid    ! request id
  
   !----- local -----
   integer(SHR_KIND_IN) :: rcode                    ! return code
   integer(SHR_KIND_IN) :: status(MPI_STATUS_SIZE)  ! request status
   logical                   :: async        ! flag for async mpi

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_recv_c) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: receive the character array "arr(n)" from task "tid" with id "msgid".
!-------------------------------------------------------------------------------

   async   = .false. ; if (PRESENT(async_in)) async   = async_in

!  write(*,F00) 'recving char msg ', msgid,' fm model ',tid,' size=',n
   if (async) then
     call mpi_irecv(arr, n, MPI_CHARACTER, tid, msgid, LIS_mpi_comm,reqid, rcode)
   else
     call mpi_recv (arr, n, MPI_CHARACTER, tid, msgid, LIS_mpi_comm,status,rcode)
   endif
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR: id=',msgid,'  error code: ' ,rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_recv_c

!===============================================================================
         
SUBROUTINE clm2_shr_msg_wait(reqid)

   !----- arguments -----
   integer(SHR_KIND_IN) :: reqid    ! request id
  
   !----- local -----
   integer(SHR_KIND_IN) :: rcode                    ! return code
   integer(SHR_KIND_IN) :: status(MPI_STATUS_SIZE)  ! request status

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_wait) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: receive the character array "arr(n)" from task "tid" with id "msgid".
!-------------------------------------------------------------------------------


!  write(*,F00) 'waiting on request ', reqid
   call mpi_wait(reqid,status,rcode)
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR:  error code: ' ,rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_wait

!===============================================================================

SUBROUTINE clm2_shr_msg_test(reqid,flag)

   !----- arguments -----
   integer(SHR_KIND_IN) :: reqid    ! request id
   logical              :: flag     ! flag for success
  
   !----- local -----
   integer(SHR_KIND_IN) :: rcode                    ! return code
   integer(SHR_KIND_IN) :: status(MPI_STATUS_SIZE)  ! request status

   !--- formats ---
   character(len=*),parameter :: F00 = "('(clm2_shr_msg_test) ',a,i4)"

!-------------------------------------------------------------------------------
! PURPOSE: receive the character array "arr(n)" from task "tid" with id "msgid".
!-------------------------------------------------------------------------------


!  write(*,F00) 'waiting on request ', reqid
   call mpi_test(reqid,flag,status,rcode)
   if (rcode /= MPI_SUCCESS) then
      write(*,F00) 'MPI ERROR:  error code: ' ,rcode
      write(*,F00) 'WARNING: calling clm2_shr_sys_abort()'
      call clm2_shr_sys_abort()
   end if
!  write(*,F00) 'done'

END SUBROUTINE clm2_shr_msg_test

!===============================================================================

#endif
#endif
END MODULE clm2_shr_msg_mod




