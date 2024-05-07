!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module hyssib_lsmMod
!BOP
!
! !MODULE: hyssib_lsmMod
!
! !DESCRIPTION:
!  This module provides the definition of derived data type used to 
!  control the operation of Hyssib LSM. It also provides the entry method
!  for the initialization of Hyssib-specific variables. The derived
!  data type {\tt hyssib\_struc} includes the variables that specify
!  the runtime options and other control variables as described below:
! \begin{description}
!  \item[rfile]
!    name of the hyssib restart file
!  \item[vfile]
!    name of the static vegetation parameter table
!  \item[topostdfile]
!    name of the file for std. dev. of topography
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[hyssibopen]
!    variable to keep track of opened files
!  \item[numout]
!    number of output times 
!  \item[nvegp]
!    number of static vegetation parameters in the table
!  \item[nvegip]
!    number of vegetation based albedo parameters in the table
!  \item[npr]
!    number of prognostic state variables used in data assimilation
!  \item[zh]
!   reference height of T and q forcing
!  \item[zm]
!   reference height of u and v forcing
!  \item[initsm]
!    initial soil moisture for a cold start run
!  \item[initTemp]
!    initial soil temperature for a cold start run
!  \item[rstInterval]
!   restart writing interval
!  \item[hyssib]
!   hyssib LSM tile specific variables
! \end{description} 
!
! !REVISION HISTORY:
!     Apr 2003: Sujay Kumar, Initial Code
!  21 Feb 2004: David Mocko, Conversion from NOAH to HY-SSiB
!  25 Aug 2007: Chuck Alonge, Updates for LIS 5.0 Compliance
!  27 Oct 2010: David Mocko, changes for HY-SSiB in LIS6.1
!
! !USES:
  use hyssib_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
   
   implicit none
   
   PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: hyssib_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: hyssib_struc
!EOP
  type, public :: hyssib_type_dec 

     character(len=LIS_CONST_PATH_LEN) :: rfile !Restart file
     character(len=LIS_CONST_PATH_LEN) :: vfile !Vegetation Params
     character(len=LIS_CONST_PATH_LEN) :: afile !Albedo and Radiation Params (Veg)
     character(len=LIS_CONST_PATH_LEN) :: topostdfile
     integer          :: count
     integer          :: hyssibopen    
     integer          :: numout      
     integer          :: nvegp   !Set to 20 in module 
     integer          :: nvegip  !Set to 11 in module
     integer          :: npr
     real             :: zh         
     real             :: zm
     real             :: initsm          
     real             :: initTemp        
     real             :: rstemp !Falling rain snow critical Temp (K)
     real             :: rstInterval  
     integer          :: forc_count
     real             :: ts
     type(hyssibdec), allocatable :: hyssib(:)
  end type hyssib_type_dec

  type(hyssib_type_dec), allocatable :: hyssib_struc(:)
  SAVE
contains
!BOP
!
! !ROUTINE: hyssib_lsm_ini
! \label{hyssib_lsm_ini}
!
! !INTERFACE:
    subroutine hyssib_lsm_ini()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_logMod,       only : LIS_verify
! !DESCRIPTION:
!  Reads in runtime HY-SSiB parameters, allocates memory for variables
!
!  The routines invoked are: 
!  \begin{description}
!   \item[hyssib\_readcrd](\ref{hyssib_readcrd}) \newline
!    reads the runtime options for Hyssib LSM
!  \end{description}
!EOP
    implicit none 
    integer :: n,i
    integer                 :: status
    character*3   :: fnest
  
    allocate(hyssib_struc(LIS_rc%nnest))
    call hyssib_readcrd()
    do n=1,LIS_rc%nnest
       allocate(hyssib_struc(n)%hyssib(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       hyssib_struc(n)%numout = 0 
       hyssib_struc(n)%rstemp = 273.15 

       hyssib_struc(n)%forc_count = 0 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          hyssib_struc(n)%hyssib(i)%tair = 0 
          hyssib_struc(n)%hyssib(i)%qair = 0 
          hyssib_struc(n)%hyssib(i)%swdown = 0 
          hyssib_struc(n)%hyssib(i)%lwdown = 0
          hyssib_struc(n)%hyssib(i)%uwind = 0
          hyssib_struc(n)%hyssib(i)%vwind = 0 
          hyssib_struc(n)%hyssib(i)%psurf = 0 
          hyssib_struc(n)%hyssib(i)%rainf_in = 0 
          hyssib_struc(n)%hyssib(i)%rainf_cp = 0 
       enddo

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, hyssib_struc(n)%ts)

       write(fnest,'(i3.3)') n

       call LIS_registerAlarm("Hyssib model alarm "//trim(fnest),&
            hyssib_struc(n)%ts,&
            hyssib_struc(n)%ts)

       call LIS_registerAlarm("Hyssib restart alarm "//trim(fnest),&
            hyssib_struc(n)%ts,&
            hyssib_struc(n)%rstInterval)

    enddo

  end subroutine hyssib_lsm_ini

end module hyssib_lsmMod

