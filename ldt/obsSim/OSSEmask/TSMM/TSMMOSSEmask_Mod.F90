!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module TSMMOSSEmask_Mod
!BOP
! 
! !MODULE: TSMMOSSEmask_Mod
! 
! !DESCRIPTION: 
!  This module handles the use of external orbital masks from
!  TSMM instrument to be applied to the observations created
!  for OSSEs
!
! !REVISION HISTORY: 
!  02 Jul 2021    Sujay Kumar  Initial Specification
!
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: TSMMOSSEmask_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: TSMMOSSEmaskData
!
!EOP
  
  type, public :: TSMMOSSEmaskDatadec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     real          :: datares
     real          :: gridDesci(20)
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

  end type TSMMOSSEmaskDatadec

  type(TSMMOSSEmaskDatadec)  :: TSMMOSSEmaskData

contains

!BOP
! !ROUTINE: TSMMOSSEmask_init
! \label{TSMMOSSEmask_init}
! 
! !INTERFACE: 
  subroutine TSMMOSSEmask_init()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_logMod

    implicit none
! 
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of 
! TSMM observation masks
! 
!EOP

    integer                 :: n
    integer                 :: rc
    real                    :: datares
    real                    :: run_dd(6)
    real                    :: gridDesci(20)
    real                    :: cornerlat1,cornerlat2
    real                    :: cornerlon1,cornerlon2
    
    n = 1


    call ESMF_ConfigGetAttribute(LDT_config,TSMMOSSEmaskData%odir, &
         label="TSMM OSSE mask directory:",rc=rc)
    call LDT_verify(rc,'TSMM OSSE mask directory: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,datares, &
         label="TSMM OSSE mask data resolution (km):",rc=rc)
    if(rc.ne.0) then
       write(LDT_logunit,*) '[ERR] TSMM OSSE mask data resolution (km): not defined'
       write(LDT_logunit,*) "[ERR] supported options are '1,5,25'"
       call LDT_endrun()
    endif
    
    if(datares.eq.1) then
       cornerlat1 = max(-90.0,nint((LDT_rc%gridDesc(n,4)+90.0)/0.01)*0.01-&
            90.0-5.0*0.01)
       cornerlon1 = max(-180.0,nint((LDT_rc%gridDesc(n,5)+180.0)/0.01)*0.01-&
            180.0-5.0*0.01)
       cornerlat2 = min(89.99,nint((LDT_rc%gridDesc(n,7)+90.0)/0.01)*0.01-&
            90.0+5.0*0.01)
       cornerlon2 = min(179.99,nint((LDT_rc%gridDesc(n,8)+180.0)/0.01)*0.01-&
            179.99+5.0*0.01)
       run_dd(1) = cornerlat1
       run_dd(2) = cornerlon1
       run_dd(3) = cornerlat2
       run_dd(4) = cornerlon2
       run_dd(5) = 0.01
       run_dd(6) = 0.01

       TSMMOSSEmaskData%nc = nint((cornerlon2-cornerlon1)/0.01)+1
       TSMMOSSEmaskData%nr = nint((cornerlat2-cornerlat1)/0.01)+1
       TSMMOSSEmaskData%datares = 0.01
    elseif(datares.eq.5) then
       TSMMOSSEmaskData%nc = 7200
       TSMMOSSEmaskData%nr = 3600
       run_dd(1) = -90.0
       run_dd(2) = -180.0
       run_dd(3) = 89.95
       run_dd(4) = 179.95
       run_dd(5) = 0.05
       run_dd(6) = 0.05
       TSMMOSSEmaskData%datares = 0.05
    elseif(datares.eq.25) then
       TSMMOSSEmaskData%nc = 1440
       TSMMOSSEmaskData%nr = 720
       run_dd(1) = -89.75
       run_dd(2) = -180.0
       run_dd(3) = 90
       run_dd(4) = 179.75
       run_dd(5) = 0.25
       run_dd(6) = 0.25
       TSMMOSSEmaskData%datares = 0.25
    endif

    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = TSMMOSSEmaskData%nc
    gridDesci(3) = TSMMOSSEmaskData%nr
    gridDesci(4) = run_dd(1)
    gridDesci(5) = run_dd(2)
    gridDesci(6) = 128
    gridDesci(7) = run_dd(3)
    gridDesci(8) = run_dd(4)
    gridDesci(9) = run_dd(5)
    gridDesci(10) = run_dd(6)
    gridDesci(20) = 64

    TSMMOSSEmaskData%gridDesci = gridDesci

!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,TSMMOSSEmaskData%datares)) then 

       allocate(TSMMOSSEmaskData%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(TSMMOSSEmaskData%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(TSMMOSSEmaskData%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(TSMMOSSEmaskData%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(TSMMOSSEmaskData%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(TSMMOSSEmaskData%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(TSMMOSSEmaskData%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(TSMMOSSEmaskData%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            TSMMOSSEmaskData%n11, &
            TSMMOSSEmaskData%n12, TSMMOSSEmaskData%n21, &
            TSMMOSSEmaskData%n22, TSMMOSSEmaskData%w11, &
            TSMMOSSEmaskData%w12, TSMMOSSEmaskData%w21, &
            TSMMOSSEmaskData%w22)

    else

       allocate(TSMMOSSEmaskData%n11(&
            TSMMOSSEmaskData%nc*&
            TSMMOSSEmaskData%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            TSMMOSSEmaskData%n11)
       
    endif
    
  end subroutine TSMMOSSEmask_init
  
end module TSMMOSSEmask_Mod
