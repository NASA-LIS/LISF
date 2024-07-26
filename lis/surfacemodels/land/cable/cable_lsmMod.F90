!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      module cable_lsmMod
!BOP
!
! !MODULE: cable_lsmMod
! \label{cable_lsmMod}
!
! !DESCRIPTION:
! This module provides the definition of the derived data type used to
! control the operation of CABLE LSM.  It also provides the entry method
! for the initialization of CABLE-specific variables.  The derived data
! type {\tt cable\_struc} includes the variables that specify the
! runtime options and other control variables as described below:
!
!  \begin{description}
!   \item[rfile]
!    name of the CABLE restart file
!   \item[vfile]
!    name of the static vegetation parameter table
!   \item[sfile]
!    name of the soil parameter table
!   \item[rstInterval]
!    restart writing interval
!   \item[cable]
!    CABLE LSM specific variables
!  \end{description}
!
! !REVISION HISTORY:
!     Apr 2003: Sujay Kumar, Initial Code
!  23 Oct 2007: Kristi Arsenault, Updated for V5.0
!  25 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  10 Jun 2012: Sujay Kumar, updated the implementation for LIS7
!
! !USES:
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use cable_module
!
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------
  public :: cable_lsm_ini
!-----------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------
  public :: cable_struc
!EOP
  type, public :: cable_type_dec
     
     character(len=LIS_CONST_PATH_LEN)              :: rfile
     character(len=LIS_CONST_PATH_LEN)              :: vfile
     character(len=LIS_CONST_PATH_LEN)              :: sfile
     character*25               :: canopyflag
     character*25               :: photosynflag
     character*25               :: soilflag
     character*25               :: slilitterflag
     character*25               :: sliisotopeflag
     character*25               :: slicoupledflag

     integer :: forc_count
     real    :: ts
     real    :: rstInterval
     
     integer :: fixedvegtype
     integer :: fixedsoiltype
     real :: fixedalbsoil
     real :: fixedco2
     real :: refheight

     logical :: verbose
     integer :: tileprint
     
     integer :: ktau
     integer :: nvegt
     integer :: nsoilt
     
     type(ESMF_Alarm) :: rstAlarm
     type(ESMF_Alarm) :: outAlarm
     type(cabledec), allocatable :: cable(:)
  end type cable_type_dec

  type(cable_type_dec), allocatable :: cable_struc(:)
  
  SAVE
contains
!BOP
!
! !ROUTINE: cable_lsm_ini
! \label{cable_lsm_ini}
!
! !INTERFACE:
  subroutine cable_lsm_ini()
! !USES:
    use LIS_coreMod,        only : LIS_rc
    use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
    use LIS_timeMgrMod
    use LIS_logMod
! !DESCRIPTION:
!
!EOP
    implicit none
    
    integer :: n
    integer                 :: i,yr, mo, da, hr, mn, ss
    integer                 :: status
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: alarmTime, alarmTime1
    type(ESMF_TimeInterval) :: alarmInterval, alarmInterval1, deltaT
    character*3             :: fnest

    allocate(cable_struc(LIS_rc%nnest))
    
    call cable_readcrd()
    do n = 1,LIS_rc%nnest
       allocate(cable_struc(n)%cable(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       
!-----------------------------------------------------------------------
! Restart Alarm
! set the alarm time to be the first 0z instance from the current time
!-----------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, cable_struc(n)%ts)

       write(fnest,'(i3.3)') n

       call LIS_registerAlarm("CABLE model alarm "//trim(fnest), &
            cable_struc(n)%ts,cable_struc(n)%ts)

       call LIS_registerAlarm("CABLE restart alarm "//trim(fnest),&
            cable_struc(n)%ts,cable_struc(n)%rstInterval)


       LIS_sfmodel_struc(n)%nsm_layers = 6
       LIS_sfmodel_struc(n)%nst_layers = 6
       allocate(LIS_sfmodel_struc(n)%lyrthk(6))

       LIS_sfmodel_struc(n)%lyrthk(1) = 0.022
       LIS_sfmodel_struc(n)%lyrthk(2) = 0.058
       LIS_sfmodel_struc(n)%lyrthk(3) = 0.154
       LIS_sfmodel_struc(n)%lyrthk(4) = 0.409
       LIS_sfmodel_struc(n)%lyrthk(5) = 1.085
       LIS_sfmodel_struc(n)%lyrthk(6) = 2.872

       LIS_sfmodel_struc(n)%ts = cable_struc(n)%ts
    enddo
  end subroutine cable_lsm_ini
  
end module cable_lsmMod
