!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module LISoutNatureRun_Mod
!BOP
! 
! !MODULE: LISoutNatureRun_Mod
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
  PUBLIC :: LISoutNatureRun_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LISoutNatureRunData
!
!EOP
  
  type, public :: LISoutNatureRunDatadec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     real          :: datares
     real          :: run_dd(8)
     character*50  :: map_proj
     character*50  :: mclass
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character(len=LDT_CONST_PATH_LEN) :: odir

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

  end type LISoutNatureRunDatadec

  type(LISoutNatureRunDatadec)  :: LISoutNatureRunData

contains

!BOP
! !ROUTINE: LISoutNatureRun_init
! \label{LISoutNatureRun_init}
! 
! !INTERFACE: 
  subroutine LISoutNatureRun_init()
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

    LISoutNatureRunData%run_dd             = LDT_rc%udef

    call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%mClass, &
         label="LIS Nature run output model class:",rc=rc)
    call LDT_verify(rc,'LIS Nature run output model class: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%format, &
         label="LIS Nature run output format:",rc=rc)
    call LDT_verify(rc,'LIS Nature run output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%wopt, &
         label="LIS Nature run output methodology:",rc=rc)
    call LDT_verify(rc,'LIS Nature run output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%wstyle, &
         label="LIS Nature run output naming style:",rc=rc)
    call LDT_verify(rc,'LIS Nature run output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%map_proj, &
         label="LIS Nature run output map projection:",rc=rc)
    call LDT_verify(rc,'LIS Nature run output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%nest, &
         label="LIS Nature run output nest index:",rc=rc)
    call LDT_verify(rc,'LIS Nature run output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%odir, &
         label="LIS Nature run output directory:",rc=rc)
    call LDT_verify(rc,'LIS Nature run output directory: not defined')

    if(LISoutNatureRunData%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(6),rc=rc)       
       
       LISoutNatureRunData%datares = min(LISoutNatureRunData%run_dd(5),LISoutNatureRunData%run_dd(6))

       LISoutNatureRunData%nc    = (nint((LISoutNatureRunData%run_dd(4)-LISoutNatureRunData%run_dd(2))/&
            LISoutNatureRunData%run_dd(5))) + 1
       LISoutNatureRunData%nr    = (nint((LISoutNatureRunData%run_dd(3)-LISoutNatureRunData%run_dd(1))/&
            LISoutNatureRunData%run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = LISoutNatureRunData%nc
       gridDesci(3) = LISoutNatureRunData%nr
       gridDesci(4) = LISoutNatureRunData%run_dd(1)
       gridDesci(5) = LISoutNatureRunData%run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = LISoutNatureRunData%run_dd(3)
       gridDesci(8) = LISoutNatureRunData%run_dd(4)
       gridDesci(9) = LISoutNatureRunData%run_dd(5)
       gridDesci(10) = LISoutNatureRunData%run_dd(6)
       gridDesci(20) = 64

    elseif(LISoutNatureRunData%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS Nature run domain y-dimension size: not defined')

       LISoutNatureRunData%datares = LISoutNatureRunData%run_dd(6)/100.0

       LISoutNatureRunData%nc    = LISoutNatureRunData%run_dd(7)
       LISoutNatureRunData%nr    = LISoutNatureRunData%run_dd(8)

       gridDesci = 0
       gridDesci(1) = 3 
       gridDesci(2) = LISoutNatureRunData%nc
       gridDesci(3) = LISoutNatureRunData%nr
       gridDesci(4) = LISoutNatureRunData%run_dd(1)
       gridDesci(5) = LISoutNatureRunData%run_dd(2)
       gridDesci(6) = 8
       gridDesci(7) = LISoutNatureRunData%run_dd(4)
       gridDesci(8) = LISoutNatureRunData%run_dd(6)
       gridDesci(9) = LISoutNatureRunData%run_dd(6)
       gridDesci(10) = LISoutNatureRunData%run_dd(3)
       gridDesci(11) = LISoutNatureRunData%run_dd(5)
       gridDesci(20) = 64

    elseif(LISoutNatureRunData%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS Nature run domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutNatureRunData%run_dd(8),rc=rc)

    endif

!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,LISoutNatureRunData%datares)) then 

       allocate(LISoutNatureRunData%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutNatureRunData%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutNatureRunData%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutNatureRunData%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(LISoutNatureRunData%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutNatureRunData%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutNatureRunData%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutNatureRunData%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            LISoutNatureRunData%n11, &
            LISoutNatureRunData%n12, LISoutNatureRunData%n21, &
            LISoutNatureRunData%n22, LISoutNatureRunData%w11, &
            LISoutNatureRunData%w12, LISoutNatureRunData%w21, &
            LISoutNatureRunData%w22)

    else

       allocate(LISoutNatureRunData%n11(&
            LISoutNatureRunData%nc*&
            LISoutNatureRunData%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            LISoutNatureRunData%nc*LISoutNatureRunData%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            LISoutNatureRunData%n11)
       
    endif
    
  end subroutine LISoutNatureRun_init
  
end module LISoutNatureRun_Mod
