!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISlsmPrecip_obsMod
!BOP
! 
! !MODULE: LISlsmPrecip_obsMod
! 
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation 
!  output as "observations" for data assimilation. This
!  plugin is typically used to handle the computations of 
!  scaling factors such as cumulative distribution function
!  (CDF) for use in DA
! 
! !REVISION HISTORY: 
!  16 Nov 2021    Mahdi Navari  Initial Specification (based on LISlsmSM_obsMod)
!

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISlsmPrecip_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lsmprecipobs
!
!EOP
  
  type, public :: lsmprecipobsdec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     real          :: datares
     real          :: run_dd(8)
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character*100 :: odir
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

  end type lsmprecipobsdec

  type(lsmprecipobsdec)  :: lsmprecipobs

contains

!BOP
! !ROUTINE: LISlsmPrecip_obsInit
! \label{LISlsmPrecip_obsInit}
! 
! !INTERFACE: 
  subroutine LISlsmPrecip_obsInit()
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

    lsmprecipobs%run_dd             = LDT_rc%udef

    lsmprecipobs%security_class     = ''
    lsmprecipobs%distribution_class = ''
    lsmprecipobs%data_category      = ''
    lsmprecipobs%area_of_data       = ''
    lsmprecipobs%write_interval     = ''

    call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%format, &
         label="LIS precipitation output format:",rc=rc)
    call LDT_verify(rc,'LIS precipitation output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%wopt, &
         label="LIS precipitation output methodology:",rc=rc)
    call LDT_verify(rc,'LIS precipitation output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%wstyle, &
         label="LIS precipitation output naming style:",rc=rc)
    call LDT_verify(rc,'LIS precipitation output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%map_proj, &
         label="LIS precipitation output map projection:",rc=rc)
    call LDT_verify(rc,'LIS precipitation output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%nest, &
         label="LIS precipitation output nest index:",rc=rc)
    call LDT_verify(rc,'LIS precipitation output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%odir, &
         label="LIS precipitation output directory:",rc=rc)
    call LDT_verify(rc,'LIS precipitation output directory: not defined')

    ! WMO-convention specific identifiers
    if ( lsmprecipobs%wstyle == "WMO convention") then 
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%security_class, &
       label="LIS precipitation security class:",rc=rc)
       call LDT_verify(rc,'LIS precipitation security class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%distribution_class, &
       label="LIS precipitation distribution class:",rc=rc)
       call LDT_verify(rc,'LIS precipitation distribution class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%data_category, &
       label="LIS precipitation data category:",rc=rc)
       call LDT_verify(rc,'LIS precipitation data category: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%area_of_data, &
       label="LIS precipitation area of data:",rc=rc)
       call LDT_verify(rc,'LIS precipitation area of data: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%write_interval, &
       label="LIS precipitation write interval:",rc=rc)
       call LDT_verify(rc,'LIS precipitation write interval: not defined')
    endif

    if(lsmprecipobs%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(6),rc=rc)       
       
       lsmprecipobs%datares = min(lsmprecipobs%run_dd(5),lsmprecipobs%run_dd(6))

       lsmprecipobs%nc    = (nint((lsmprecipobs%run_dd(4)-lsmprecipobs%run_dd(2))/&
            lsmprecipobs%run_dd(5))) + 1
       lsmprecipobs%nr    = (nint((lsmprecipobs%run_dd(3)-lsmprecipobs%run_dd(1))/&
            lsmprecipobs%run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = lsmprecipobs%nc
       gridDesci(3) = lsmprecipobs%nr
       gridDesci(4) = lsmprecipobs%run_dd(1)
       gridDesci(5) = lsmprecipobs%run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = lsmprecipobs%run_dd(3)
       gridDesci(8) = lsmprecipobs%run_dd(4)
       gridDesci(9) = lsmprecipobs%run_dd(5)
       gridDesci(10) = lsmprecipobs%run_dd(6)
       gridDesci(20) = 64

    elseif(lsmprecipobs%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS precipitation domain y-dimension size: not defined')

       lsmprecipobs%datares = lsmprecipobs%run_dd(6)/100.0

       lsmprecipobs%nc    = lsmprecipobs%run_dd(7)
       lsmprecipobs%nr    = lsmprecipobs%run_dd(8)

       gridDesci = 0
       gridDesci(1) = 3 
       gridDesci(2) = lsmprecipobs%nc
       gridDesci(3) = lsmprecipobs%nr
       gridDesci(4) = lsmprecipobs%run_dd(1)
       gridDesci(5) = lsmprecipobs%run_dd(2)
       gridDesci(6) = 8
       gridDesci(7) = lsmprecipobs%run_dd(4)
       gridDesci(8) = lsmprecipobs%run_dd(6)
       gridDesci(9) = lsmprecipobs%run_dd(6)
       gridDesci(10) = lsmprecipobs%run_dd(3)
       gridDesci(11) = lsmprecipobs%run_dd(5)
       gridDesci(20) = 64

    elseif(lsmprecipobs%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS precipitation domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmprecipobs%run_dd(8),rc=rc)

    endif

!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,lsmprecipobs%datares)) then 

       allocate(lsmprecipobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmprecipobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmprecipobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmprecipobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(lsmprecipobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmprecipobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmprecipobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmprecipobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            lsmprecipobs%n11, &
            lsmprecipobs%n12, lsmprecipobs%n21, &
            lsmprecipobs%n22, lsmprecipobs%w11, &
            lsmprecipobs%w12, lsmprecipobs%w21, &
            lsmprecipobs%w22)

    else

       allocate(lsmprecipobs%n11(&
            lsmprecipobs%nc*&
            lsmprecipobs%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            lsmprecipobs%nc*lsmprecipobs%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            lsmprecipobs%n11)
       
    endif
    
!  which variable we want in the DA obs computations. 
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%totalprecip_obs, "kg m-2",1,1)
    LDT_DAobsData(1)%totalprecip_obs%selectStats = 1

  end subroutine LISlsmPrecip_obsInit
  
end module LISlsmPrecip_obsMod
