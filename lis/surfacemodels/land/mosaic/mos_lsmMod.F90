!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module mos_lsmMod
!BOP
!
! !ROUTINE: mos_lsmMod
!
! !DESCRIPTION:
!  
! This module provides the definition of derived data type used to 
! control the operation of Mosaic LSM. It also provides the entry method
! for the initialization of Mosaic-specific variables. The derived
! data type {\tt mos\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
!  \item[mos\_rfile]
!    name of the mosaic restart file
!  \item[mos\_vfile]
!    name of the static vegetation parameter table
!  \item[mos\_sfile]
!    name of the soil parameter table
!  \item[mos\_mvfile]
!    name of monthly vegetation parameter file
!  \item[mos\_pfile]
!    name of general parameter file (levs, depths)
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[mosopen]
!    variable to keep track of opened files
!  \item[numout]
!    number of output times 
!  \item[mos\_nvegp]
!    number of static vegetation parameters in the table
!  \item[mos\_nmvegp]
!    number of monthly vegetation parameters in the table
!  \item[mos\_ism]
!    initial soil moisture for a cold start run
!  \item[mos\_it]
!    initial soil temperature for a cold start run
!  \item[outInterval]
!    output writing interval
!  \item[rstInterval]
!   restart writing interval
!  \item[forcing\_z]
!   flag indicating whether to use the observation height from the forcing
!  \item[forcing\_ch]
!   flag indicating whether to use aerodynamic conductance field from forcing
!  \item[mos]
!   Mosaic LSM specific variables
! \end{description} 

!
! !REVISION HISTORY:
! Jun 2003; Jon Gottschalck, Initial Code
! Sept 2007: Sujay Kumar, Upgraded the data structures for LIS 5.0
!
! !INTERFACE:

! !USES:        
  use mos_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: mos_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: mos_struc
!EOP       
  type, public ::  mos_type_dec
     integer               :: mosopen         
     integer               :: numout         
     integer               :: mos_nvegp  
     integer               :: mos_nmvegp      
     integer               :: mos_nsoilp      
     character(len=LIS_CONST_PATH_LEN) :: mos_rfile 
     character(len=LIS_CONST_PATH_LEN) :: mos_vfile 
     character(len=LIS_CONST_PATH_LEN) :: mos_pfile 
     character(len=LIS_CONST_PATH_LEN) :: mos_mvfile 
     character(len=LIS_CONST_PATH_LEN) :: mos_sfile
     integer               :: mos_nstxts
     integer               :: usedsoilmap
     real                  :: dpthlyr1
     real                  :: dpthlyr2
     real                  :: dpthlyr3
     real                  :: mos_ism            
     real                  :: mos_it
     real                  :: rstInterval     
     integer               :: count
     integer               :: forcing_z
     integer               :: forcing_ch
     integer               :: forc_count
     real                  :: ts
     type(mosdec), allocatable :: mos(:)
  end type mos_type_dec
  
  type(mos_type_dec), allocatable :: mos_struc(:)
  
  SAVE
contains
!BOP  
!
! !ROUTINE: mos_lsm_ini
! \label{mos_lsm_ini}
!
! !INTERFACE:
  subroutine mos_lsm_ini()
! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_logMod,       only : LIS_verify
! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for mosaic-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for mosaic from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readmoscrd](\ref{readmoscrd}) \newline
!    reads the runtime options for Mosaic LSM
!  \end{description}
!EOP  
    integer :: n,i
    integer                 :: status
    character*3         :: fnest

    allocate(mos_struc(LIS_rc%nnest))
    call readmoscrd()
    
    do n=1,LIS_rc%nnest
       allocate(mos_struc(n)%mos(LIS_rc%npatch(n,LIS_rc%lsm_index)))

       mos_struc(n)%forc_count = 0 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          mos_struc(n)%mos(i)%tair = 0 
          mos_struc(n)%mos(i)%qair = 0 
          mos_struc(n)%mos(i)%swdown = 0 
          mos_struc(n)%mos(i)%lwdown = 0
          mos_struc(n)%mos(i)%uwind = 0
          mos_struc(n)%mos(i)%vwind = 0 
          mos_struc(n)%mos(i)%psurf = 0 
          mos_struc(n)%mos(i)%rainf = 0 
          mos_struc(n)%mos(i)%rainf_c = 0 
          mos_struc(n)%mos(i)%obsz = 0 
          mos_struc(n)%mos(i)%water1 = 0 
          mos_struc(n)%mos(i)%water2 = 0 
          mos_struc(n)%mos(i)%water3 = 0 
       enddo

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, mos_struc(n)%ts)

       write(fnest,'(i3.3)') n    
       call LIS_registerAlarm("Mosaic model alarm "//trim(fnest),&
            mos_struc(n)%ts,&
            mos_struc(n)%ts)

       call LIS_registerAlarm("Mosaic restart alarm "//trim(fnest),&
            mos_struc(n)%ts,&
            mos_struc(n)%rstInterval)

       LIS_sfmodel_struc(n)%nsm_layers = 3
       LIS_sfmodel_struc(n)%nst_layers = 3
       allocate(LIS_sfmodel_struc(n)%lyrthk(3))
       LIS_sfmodel_struc(n)%lyrthk(1) = 2.0
       LIS_sfmodel_struc(n)%lyrthk(2) = 98.0
       LIS_sfmodel_struc(n)%lyrthk(3) = 200.0

       LIS_sfmodel_struc(n)%ts = mos_struc(n)%ts
    enddo

  end subroutine mos_lsm_ini
end module mos_lsmMod



