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
module noah271_lsmMod
!BOP
!
! !MODULE: noah271_lsmMod
!
! !DESCRIPTION:
!  
! This module provides the definition of derived data type used to 
! control the operation of Noah2.7.1 LSM. It also provides the entry method
! for the initialization of Noah2.7.1-specific variables. The derived
!  data type {\tt noah2.7.1\_struc} includes the variables that specify
!  the runtime options and other control variables as described below:
!
! \begin{description}
!  \item[rfile]
!    name of the noah2.7.1 restart file
!  \item[vfile]
!    name of the static vegetation parameter table
!  \item[sfile]
!    name of the soil parameter table
!  \item[useptf]
!    whether PTFs should be used to map soil properties
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[noah271open]
!    variable to keep track of opened files
!  \item[soilscheme]
!    type of soil mapping scheme usd (1-zobler,2-statsgo)
!  \item[numout]
!    number of output times 
!  \item[nslay]
!    number of soil layers
!  \item[nvegp]
!    number of static vegetation parameters in the table
!  \item[nsoilp]
!    number of soil parameters in the table
!  \item[nstxts]
!    number of soil texture classes in the classification scheme
!  \item[param\_rst]
!    flag to indicate if parameters are to be read from a restart file
!    (e.g. calibrated parameters output by the optUE codes)
!  \item[prstfile]
!    parameter restart filename (e.g. GA output file)
!  \item[forcing\_ch]
!    flag indicating if the heat exchange coefficient is present in the 
!    forcing inputs
!  \item[inittemp]
!    initial soil temperatures for a cold start run
!  \item[initsm]
!    initial total soil moisture for a cold start run
!  \item[initsmliq]
!    initial liquid soil moisture for a cold start run
!  \item[initskintemp]
!    initial value of skin temperature for a cold start run
!  \item[initcanopywater]
!    initial value of canopy water for a cold start run
!  \item[initsnowdepth]
!    initial value of snow depth for a cold start run
!  \item[initsnoweqiv]
!    initial value of snow water equivalent for a cold start run
!  \item[outInterval]
!    output writing interval
!  \item[rstInterval]
!   restart writing interval
!  \item[rstAlarm]
!   restart alarm object 
!  \item[outAlarm]
!   output alarm object 
!  \item[zh]
!   reference height of T and q forcing
!  \item[zm]
!   reference height of u and v forcing
!  \item[forcing\_ch]
!   flag indicating whether to use aerodynamic conductance field from forcing
!  \item[lyrthk]
!   thickness of soil layers
!  \item[noah]
!   Noah2.7.1 LSM specific variables
! \end{description} 
!
! !REVISION HISTORY:
! Apr 2003; Sujay Kumar, Initial Code
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
!
! !USES:        
  use noah271_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noah271_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noah271_struc
!EOP
  type, public ::  noah271_type_dec 

     character(len=LIS_CONST_PATH_LEN) :: rfile
     character(len=LIS_CONST_PATH_LEN) :: vfile
     character(len=LIS_CONST_PATH_LEN) :: sfile
     integer                    :: useptf
     integer                    :: count
     integer                    :: noah271open
     integer                    :: soilscheme
     integer                    :: numout
     integer                    :: nslay
     integer                    :: nvegp
     integer                    :: nsoilp
     integer                    :: nstxts
     real                       :: ts
     real                       :: rstInterval
     integer                    :: param_rst
     character(len=LIS_CONST_PATH_LEN) :: prstfile
     real                       :: zh
     real                       :: zm
     integer                    :: forcing_ch
     integer                    :: forc_count
     real, allocatable              :: lyrthk(:)
     real, allocatable              :: inittemp(:)
     real, allocatable              :: initsm(:)
     real, allocatable              :: initsmliq(:)
     real                       :: initskintemp
     real                       :: initcanopywater
     real                       :: initsnowdepth
     real                       :: initsnowequiv
     type(noah271dec), allocatable :: noah(:)
  end type noah271_type_dec

  type(noah271_type_dec), allocatable :: noah271_struc(:)
  SAVE
contains
!BOP
! 
! !ROUTINE: noah271_lsm_ini
! \label{noah271_lsm_ini}
! 
! !INTERFACE:
  subroutine noah271_lsm_ini()
! !USES:
    use LIS_coreMod,      only : LIS_rc
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar,  &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_logMod,       only : LIS_verify
! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for Noah2.7.1-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for Noah2.7.1 from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[noah271\_readcrd](\ref{noah271_readcrd}) \newline
!    reads the runtime options for Noah2.7.1 LSM
!  \end{description}
!EOP
    implicit none
    integer                 :: i,n
    integer                 :: status
    character*3   :: fnest

    allocate(noah271_struc(LIS_rc%nnest))
    call noah271_readcrd()
    do n=1,LIS_rc%nnest
       allocate(noah271_struc(n)%noah(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          allocate(noah271_struc(n)%noah(i)%stc(noah271_struc(n)%nslay))
          allocate(noah271_struc(n)%noah(i)%smc(noah271_struc(n)%nslay))
          allocate(noah271_struc(n)%noah(i)%sh2o(noah271_struc(n)%nslay))
          allocate(noah271_struc(n)%noah(i)%relsmc(noah271_struc(n)%nslay))
       enddo
       noah271_struc(n)%numout = 0 

       noah271_struc(n)%forc_count = 0 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noah271_struc(n)%noah(i)%tair = 0 
          noah271_struc(n)%noah(i)%qair = 0 
          noah271_struc(n)%noah(i)%swdown = 0 
          noah271_struc(n)%noah(i)%lwdown = 0
          noah271_struc(n)%noah(i)%uwind = 0
          noah271_struc(n)%noah(i)%vwind = 0 
          noah271_struc(n)%noah(i)%psurf = 0 
          noah271_struc(n)%noah(i)%rainf = 0 
          noah271_struc(n)%noah(i)%rainf_c = 0 
          noah271_struc(n)%noah(i)%ch = 0 
       enddo

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, noah271_struc(n)%ts)

       write(fnest,'(i3.3)') n

       call LIS_registerAlarm("Noah271 model alarm "//trim(fnest),&
            noah271_struc(n)%ts,&
            noah271_struc(n)%ts)

       call LIS_registerAlarm("Noah271 restart alarm "//trim(fnest),&
            noah271_struc(n)%ts,&
            noah271_struc(n)%rstInterval)

       ! Initialize min/max values to implausible values.
       noah271_struc(n)%noah(:)%tair_agl_min = 999.0
       noah271_struc(n)%noah(:)%xice = LIS_rc%udef

       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noah271_struc(n)%noah(i)%cqs2 = LIS_rc%udef
          noah271_struc(n)%noah(i)%qsfc = LIS_rc%udef
          noah271_struc(n)%noah(i)%sca = 0.0
       enddo

       LIS_sfmodel_struc(n)%nsm_layers = noah271_struc(n)%nslay
       LIS_sfmodel_struc(n)%nst_layers = noah271_struc(n)%nslay
       allocate(LIS_sfmodel_struc(n)%lyrthk(noah271_struc(n)%nslay))
       do i = 1,noah271_struc(n)%nslay
          LIS_sfmodel_struc(n)%lyrthk(i) = noah271_struc(n)%lyrthk(i)*100.0
       enddo
       LIS_sfmodel_struc(n)%ts = noah271_struc(n)%ts

    enddo
  end subroutine noah271_lsm_ini

end module noah271_lsmMod
