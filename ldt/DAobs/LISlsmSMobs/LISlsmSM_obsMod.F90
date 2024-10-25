!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISlsmSM_obsMod
!BOP
! 
! !MODULE: LISlsmSM_obsMod
! 
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation 
!  output as "observations" for data assimilation. This
!  plugin is typically used to handle the computations of 
!  scaling factors such as cumulative distribution function
!  (CDF) for use in DA
! 
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
!
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISlsmSM_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lsmsmobs
!
!EOP
  
  type, public :: lsmsmobsdec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     real          :: datares
     real          :: run_dd(8)
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character(len=LDT_CONST_PATH_LEN) :: odir
     character*20  :: security_class
     character*20  :: distribution_class
     character*20  :: data_category
     character*20  :: area_of_data
     character*20  :: write_interval
!--------------------------------------------------------
!  interpolation/upscaling weights
!--------------------------------------------------------
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,  allocatable     :: w11(:)
     real,  allocatable     :: w12(:)
     real,  allocatable     :: w21(:)
     real,  allocatable     :: w22(:)

  end type lsmsmobsdec

  type(lsmsmobsdec)  :: lsmsmobs

contains

!BOP
! !ROUTINE: LISlsmSM_obsInit
! \label{LISlsmSM_obsInit}
! 
! !INTERFACE: 
  subroutine LISlsmSM_obsInit()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_logMod

    implicit none
! 
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model output (from a LIS simulation) as observations.  
! 
!EOP
    integer                 :: n
    integer                 :: rc
    real                    :: gridDesci(20)

    n = 1

    lsmsmobs%run_dd             = LDT_rc%udef

    lsmsmobs%security_class     = ''
    lsmsmobs%distribution_class = ''
    lsmsmobs%data_category      = ''
    lsmsmobs%area_of_data       = ''
    lsmsmobs%write_interval     = ''

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%format, &
         label="LIS soil moisture output format:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%wopt, &
         label="LIS soil moisture output methodology:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%wstyle, &
         label="LIS soil moisture output naming style:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%map_proj, &
         label="LIS soil moisture output map projection:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%nest, &
         label="LIS soil moisture output nest index:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%odir, &
         label="LIS soil moisture output directory:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output directory: not defined')

    ! WMO-convention specific identifiers
    if ( lsmsmobs%wstyle == "WMO convention") then 
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%security_class, &
       label="LIS soil moisture security class:",rc=rc)
       call LDT_verify(rc,'LIS soil moisture security class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%distribution_class, &
       label="LIS soil moisture distribution class:",rc=rc)
       call LDT_verify(rc,'LIS soil moisture distribution class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%data_category, &
       label="LIS soil moisture data category:",rc=rc)
       call LDT_verify(rc,'LIS soil moisture data category: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%area_of_data, &
       label="LIS soil moisture area of data:",rc=rc)
       call LDT_verify(rc,'LIS soil moisture area of data: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%write_interval, &
       label="LIS soil moisture write interval:",rc=rc)
       call LDT_verify(rc,'LIS soil moisture write interval: not defined')
    endif

    if(lsmsmobs%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(6),rc=rc)       
       
       lsmsmobs%datares = min(lsmsmobs%run_dd(5),lsmsmobs%run_dd(6))

       lsmsmobs%nc    = (nint((lsmsmobs%run_dd(4)-lsmsmobs%run_dd(2))/&
            lsmsmobs%run_dd(5))) + 1
       lsmsmobs%nr    = (nint((lsmsmobs%run_dd(3)-lsmsmobs%run_dd(1))/&
            lsmsmobs%run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = lsmsmobs%nc
       gridDesci(3) = lsmsmobs%nr
       gridDesci(4) = lsmsmobs%run_dd(1)
       gridDesci(5) = lsmsmobs%run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = lsmsmobs%run_dd(3)
       gridDesci(8) = lsmsmobs%run_dd(4)
       gridDesci(9) = lsmsmobs%run_dd(5)
       gridDesci(10) = lsmsmobs%run_dd(6)
       gridDesci(20) = 64

    elseif(lsmsmobs%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain y-dimension size: not defined')

       lsmsmobs%datares = lsmsmobs%run_dd(6)/100.0

       lsmsmobs%nc    = lsmsmobs%run_dd(7)
       lsmsmobs%nr    = lsmsmobs%run_dd(8)

       gridDesci = 0
       gridDesci(1) = 3 
       gridDesci(2) = lsmsmobs%nc
       gridDesci(3) = lsmsmobs%nr
       gridDesci(4) = lsmsmobs%run_dd(1)
       gridDesci(5) = lsmsmobs%run_dd(2)
       gridDesci(6) = 8
       gridDesci(7) = lsmsmobs%run_dd(4)
       gridDesci(8) = lsmsmobs%run_dd(6)
       gridDesci(9) = lsmsmobs%run_dd(6)
       gridDesci(10) = lsmsmobs%run_dd(3)
       gridDesci(11) = lsmsmobs%run_dd(5)
       gridDesci(20) = 64

    elseif(lsmsmobs%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil moisture domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmsmobs%run_dd(8),rc=rc)

    endif

!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,lsmsmobs%datares)) then 

       allocate(lsmsmobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmsmobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmsmobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmsmobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(lsmsmobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmsmobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmsmobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmsmobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            lsmsmobs%n11, &
            lsmsmobs%n12, lsmsmobs%n21, &
            lsmsmobs%n22, lsmsmobs%w11, &
            lsmsmobs%w12, lsmsmobs%w21, &
            lsmsmobs%w22)

    else

       allocate(lsmsmobs%n11(&
            lsmsmobs%nc*&
            lsmsmobs%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            lsmsmobs%nc*lsmsmobs%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            lsmsmobs%n11)
       
    endif
    
!  which variable we want in the DA obs computations. 
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%soilmoist_obs, "m3/m3",1,1)
    LDT_DAobsData(1)%soilmoist_obs%selectStats = 1

  end subroutine LISlsmSM_obsInit
  
end module LISlsmSM_obsMod
