!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module LISoutOSSEmask_Mod
!BOP
! 
! !MODULE: LISoutOSSEmask_Mod
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
  PUBLIC :: LISoutOSSEmask_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LISoutOSSEmaskData
!
!EOP
  
  type, public :: LISoutOSSEmaskDatadec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     real          :: datares
     real          :: run_dd(8)
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character*50  :: type
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

  end type LISoutOSSEmaskDatadec

  type(LISoutOSSEmaskDatadec)  :: LISoutOSSEmaskData

contains

!BOP
! !ROUTINE: LISoutOSSEmask_init
! \label{LISoutOSSEmask_init}
! 
! !INTERFACE: 
  subroutine LISoutOSSEmask_init()
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

    LISoutOSSEmaskData%run_dd             = LDT_rc%udef

    call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%format, &
         label="LIS output-based OSSE mask format:",rc=rc)
    call LDT_verify(rc,'LIS output-based OSSE mask format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%wopt, &
         label="LIS output-based OSSE mask methodology:",rc=rc)
    call LDT_verify(rc,'LIS output-based OSSE mask methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%wstyle, &
         label="LIS output-based OSSE mask naming style:",rc=rc)
    call LDT_verify(rc,'LIS output-based OSSE mask naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%map_proj, &
         label="LIS output-based OSSE mask map projection:",rc=rc)
    call LDT_verify(rc,'LIS output-based OSSE mask map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%nest, &
         label="LIS output-based OSSE mask nest index:",rc=rc)
    call LDT_verify(rc,'LIS output-based OSSE mask nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%odir, &
         label="LIS output-based OSSE mask directory:",rc=rc)
    call LDT_verify(rc,'LIS output-based OSSE mask directory: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%type, &
         label="LIS output-based OSSE mask type:",rc=rc)
    if(rc.ne.0) then 
       write(LDT_logunit,*)'[ERR] LIS output-based OSSE mask type: not defined'
       write(LDT_logunit,*)'[ERR] supported options are...'
       write(LDT_logunit,*)"[ERR] 'PMW soil moisture'"
       call LDT_endrun()
    endif

    if(LISoutOSSEmaskData%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(6),rc=rc)       
       
       LISoutOSSEmaskData%datares = min(LISoutOSSEmaskData%run_dd(5),LISoutOSSEmaskData%run_dd(6))

       LISoutOSSEmaskData%nc    = (nint((LISoutOSSEmaskData%run_dd(4)-LISoutOSSEmaskData%run_dd(2))/&
            LISoutOSSEmaskData%run_dd(5))) + 1
       LISoutOSSEmaskData%nr    = (nint((LISoutOSSEmaskData%run_dd(3)-LISoutOSSEmaskData%run_dd(1))/&
            LISoutOSSEmaskData%run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = LISoutOSSEmaskData%nc
       gridDesci(3) = LISoutOSSEmaskData%nr
       gridDesci(4) = LISoutOSSEmaskData%run_dd(1)
       gridDesci(5) = LISoutOSSEmaskData%run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = LISoutOSSEmaskData%run_dd(3)
       gridDesci(8) = LISoutOSSEmaskData%run_dd(4)
       gridDesci(9) = LISoutOSSEmaskData%run_dd(5)
       gridDesci(10) = LISoutOSSEmaskData%run_dd(6)
       gridDesci(20) = 64

    elseif(LISoutOSSEmaskData%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS output-based OSSE mask domain y-dimension size: not defined')

       LISoutOSSEmaskData%datares = LISoutOSSEmaskData%run_dd(6)/100.0

       LISoutOSSEmaskData%nc    = LISoutOSSEmaskData%run_dd(7)
       LISoutOSSEmaskData%nr    = LISoutOSSEmaskData%run_dd(8)

       gridDesci = 0
       gridDesci(1) = 3 
       gridDesci(2) = LISoutOSSEmaskData%nc
       gridDesci(3) = LISoutOSSEmaskData%nr
       gridDesci(4) = LISoutOSSEmaskData%run_dd(1)
       gridDesci(5) = LISoutOSSEmaskData%run_dd(2)
       gridDesci(6) = 8
       gridDesci(7) = LISoutOSSEmaskData%run_dd(4)
       gridDesci(8) = LISoutOSSEmaskData%run_dd(6)
       gridDesci(9) = LISoutOSSEmaskData%run_dd(6)
       gridDesci(10) = LISoutOSSEmaskData%run_dd(3)
       gridDesci(11) = LISoutOSSEmaskData%run_dd(5)
       gridDesci(20) = 64

    elseif(LISoutOSSEmaskData%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS output-based OSSE mask domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,LISoutOSSEmaskData%run_dd(8),rc=rc)

    endif

!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,LISoutOSSEmaskData%datares)) then 

       allocate(LISoutOSSEmaskData%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutOSSEmaskData%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutOSSEmaskData%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutOSSEmaskData%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(LISoutOSSEmaskData%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutOSSEmaskData%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutOSSEmaskData%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(LISoutOSSEmaskData%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            LISoutOSSEmaskData%n11, &
            LISoutOSSEmaskData%n12, LISoutOSSEmaskData%n21, &
            LISoutOSSEmaskData%n22, LISoutOSSEmaskData%w11, &
            LISoutOSSEmaskData%w12, LISoutOSSEmaskData%w21, &
            LISoutOSSEmaskData%w22)

    else

       allocate(LISoutOSSEmaskData%n11(&
            LISoutOSSEmaskData%nc*&
            LISoutOSSEmaskData%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            LISoutOSSEmaskData%nc*LISoutOSSEmaskData%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            LISoutOSSEmaskData%n11)
       
    endif
    
  end subroutine LISoutOSSEmask_init
  
end module LISoutOSSEmask_Mod
