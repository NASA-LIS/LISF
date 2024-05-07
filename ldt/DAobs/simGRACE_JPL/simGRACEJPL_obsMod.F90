!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: simGRACEJPL_obsMod
! 
! !DESCRIPTION: 
!  This module provides the observation plugin for the 
!  simulated GRACE observations from JPL
!   
! !REVISION HISTORY: 
!  24 Feb 2015: Sujay Kumar, Initial Specification
!
module simGRACEJPL_obsMod 
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: simGRACEJPL_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: simGRACEJPLobs
!EOP
  type, public :: simgracetwsobsdec
     integer :: nest 
     character*50  :: config
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character(len=LDT_CONST_PATH_LEN) :: odir

     integer       :: reftime
     integer       :: tdims
     character(len=LDT_CONST_PATH_LEN) :: gracedir
     logical       :: startMode
     integer       :: gracenc, gracenr
     real, allocatable :: tvals(:)
     real, allocatable :: time_bounds(:,:)
     real, allocatable :: lwe_thickness(:,:,:)
     real, allocatable :: twsavg(:,:)
     real, allocatable :: lisavg(:,:)
     integer, allocatable :: nlisavg(:,:)

     integer            :: b_syr
     integer            :: b_eyr

     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real, allocatable      :: w11(:),w12(:)
     real, allocatable      :: w21(:),w22(:)
  end type simgracetwsobsdec

  type(simgracetwsobsdec) :: simGRACEJPLobs

contains
  
!BOP
! 
! !ROUTINE: simGRACEJPL_obsInit
! \label{simGRACEJPL_obsInit}
! 
! !INTERFACE: 
  subroutine simGRACEJPL_obsinit()
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading ECV soil moisture data. 
! 
!EOP
    real                    :: run_dd(8)
    integer                 :: rc
    integer                 :: n 
    real                    :: gridDesci(20)

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%gracedir, &
         label="Simulated GRACE data directory:",rc=rc)
    call LDT_verify(rc,'Simulated GRACE data directory: not defined')

! supported options - 'default', 'follow-on', 'GRACE-2'
    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%config, &
         label="Simulated GRACE configuration:",rc=rc)
    call LDT_verify(rc,'Simulated GRACE configuration: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%odir, &
         label="LIS TWS output directory:",rc=rc)
    call LDT_verify(rc,'LIS TWS output directory: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%format, &
         label="LIS TWS output format:",rc=rc)
    call LDT_verify(rc,'LIS TWS output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%wopt, &
         label="LIS TWS output methodology:",rc=rc)
    call LDT_verify(rc,'LIS TWS output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%wstyle, &
         label="LIS TWS output naming style:",rc=rc)
    call LDT_verify(rc,'LIS TWS output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%map_proj, &
         label="LIS TWS output map projection:",rc=rc)
    call LDT_verify(rc,'LIS TWS output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%nest, &
         label="LIS TWS output nest index:",rc=rc)
    call LDT_verify(rc,'LIS TWS output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%b_syr, &
         label="Simulated GRACE baseline starting year:",rc=rc)
    call LDT_verify(rc,'simulated GRACE baseline start year: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,simGRACEJPLobs%b_eyr, &
         label="Simulated GRACE baseline ending year:",rc=rc)
    call LDT_verify(rc,'simulated GRACE baseline ending year: not defined')

    if(simGRACEJPLobs%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(6),rc=rc)       
       
    elseif(simGRACEJPLobs%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS TWS output domain y-dimension size: not defined')

    elseif(simGRACEJPLobs%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS TWS output domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(8),rc=rc)

    endif

    simGRACEJPLobs%startMode = .true. 

    simGRACEJPLobs%gracenc = 360
    simGRACEJPLobs%gracenr = 180

!  which variable we want in the DA obs computations. 
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%TWS_obs, "mm",1,1)
    LDT_DAobsData(1)%tws_obs%selectStats = 1

    call LDT_get_julhr(2006,1,1,0,0,0,simGRACEJPLobs%reftime)
    
    n =1 

    allocate(simGRACEJPLobs%lisavg(LDT_rc%lnc(n),LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%nlisavg(LDT_rc%lnc(n),LDT_rc%lnr(n)))

    simGRACEJPLobs%lisavg = 0 
    simGRACEJPLobs%nlisavg = 0 

    gridDesci = 0 
    gridDesci(1) = 0
    gridDesci(2) = 360
    gridDesci(3) = 180
    gridDesci(4) = -89.5
    gridDesci(5) = -179.5
    gridDesci(6) = 128
    gridDesci(7) = 89.5
    gridDesci(8) = 179.5
    gridDesci(9) = 1.0
    gridDesci(10) = 1.0
    gridDesci(20) = 64

    allocate(simGRACEJPLobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(simGRACEJPLobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

    call bilinear_interp_input(n, gridDesci,&
         simGRACEJPLobs%n11,&
         simGRACEJPLobs%n12,&
         simGRACEJPLobs%n21,&
         simGRACEJPLobs%n22,&
         simGRACEJPLobs%w11,&
         simGRACEJPLobs%w12,&
         simGRACEJPLobs%w21,&
         simGRACEJPLobs%w22)

    call system('mkdir -p '//(LDT_rc%odir))

    if((simGRACEJPLobs%config.eq."default").or.&
         (simGRACEJPLobs%config.eq."follow-on")) then 
       simGRACEJPLobs%tdims = 12
       
       allocate(simGRACEJPLobs%tvals(simGRACEJPLobs%tdims))
!    simGRACEJPLobs%tvals = (/31,59,90,120,151,181,212,243,273,304,334,365/)
       simGRACEJPLobs%tvals = (/0,31,59,90,120,151,181,212,243,273,304,334/)
    elseif(simGRACEJPLobs%config.eq."GRACE-2") then 
       simGRACEJPLobs%tdims = 28       
       allocate(simGRACEJPLobs%tvals(simGRACEJPLobs%tdims))
!       simGRACEJPLobs%tvals = (/12,25,38,51,64,77,90,103,116,129,142,&
!            155,168,181,194,207,220,233,246,259,272,285,298,311,324,&
!            337,350,364/)

       simGRACEJPLobs%tvals = (/0,13,26,39,52,65,78,91,104,117,130,143,&
            156,169,182,195,208,221,234,247,260,273,286,299,312,325,&
            338,351/)
    endif

  end subroutine simGRACEJPL_obsinit
     
end module simGRACEJPL_obsMod
