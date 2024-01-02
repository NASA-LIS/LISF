!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module AMSR2OSSEmask_Mod
!BOP
! 
! !MODULE: AMSR2OSSEmask_Mod
! 
! !DESCRIPTION: 
!  This module handles the use of external orbital masks from
!  AMSR2 instrument to be applied to the observations created
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
  PUBLIC :: AMSR2OSSEmask_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: AMSR2OSSEmaskData
!
!EOP
  
  type, public :: AMSR2OSSEmaskDatadec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     real          :: datares

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

  end type AMSR2OSSEmaskDatadec

  type(AMSR2OSSEmaskDatadec)  :: AMSR2OSSEmaskData

contains

!BOP
! !ROUTINE: AMSR2OSSEmask_init
! \label{AMSR2OSSEmask_init}
! 
! !INTERFACE: 
  subroutine AMSR2OSSEmask_init()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_logMod

    implicit none
! 
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of 
! AMSR2 observation masks
! 
!EOP

    integer                 :: n
    integer                 :: rc
    real                    :: datares
    real                    :: run_dd(6)
    real                    :: gridDesci(20)
    real                    :: cornerlat1, cornerlat2
    real                    :: cornerlon1, cornerlon2

    n = 1


    call ESMF_ConfigGetAttribute(LDT_config,AMSR2OSSEmaskData%odir, &
         label="AMSR2 OSSE mask directory:",rc=rc)
    call LDT_verify(rc,'AMSR2 OSSE mask directory: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,datares, &
         label="AMSR2 OSSE mask data resolution (km):",rc=rc)
    if(rc.ne.0) then
       write(LDT_logunit,*) '[ERR] AMSR2 OSSE mask data resolution (km): not defined'
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

       AMSR2OSSEmaskData%nc = nint((cornerlon2-cornerlon1)/0.01)+1
       AMSR2OSSEmaskData%nr = nint((cornerlat2-cornerlat1)/0.01)+1
       AMSR2OSSEmaskData%datares = 0.01
    elseif(datares.eq.5) then
       AMSR2OSSEmaskData%nc = 7200
       AMSR2OSSEmaskData%nr = 3600
       run_dd(1) = -90.0
       run_dd(2) = -180.0
       run_dd(3) = 89.95
       run_dd(4) = 179.95
       run_dd(5) = 0.05
       run_dd(6) = 0.05
       AMSR2OSSEmaskData%datares = 0.05
    elseif(datares.eq.25) then
       AMSR2OSSEmaskData%nc = 1440
       AMSR2OSSEmaskData%nr = 720
       run_dd(1) = -90.0
       run_dd(2) = -180.0
       run_dd(3) = 89.75
       run_dd(4) = 179.75
       run_dd(5) = 0.25
       run_dd(6) = 0.25
       AMSR2OSSEmaskData%datares = 0.25
    endif

    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = AMSR2OSSEmaskData%nc
    gridDesci(3) = AMSR2OSSEmaskData%nr
    gridDesci(4) = run_dd(1)
    gridDesci(5) = run_dd(2)
    gridDesci(6) = 128
    gridDesci(7) = run_dd(3)
    gridDesci(8) = run_dd(4)
    gridDesci(9) = run_dd(5)
    gridDesci(10) = run_dd(6)
    gridDesci(20) = 64


!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,AMSR2OSSEmaskData%datares)) then 

       allocate(AMSR2OSSEmaskData%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(AMSR2OSSEmaskData%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(AMSR2OSSEmaskData%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(AMSR2OSSEmaskData%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(AMSR2OSSEmaskData%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(AMSR2OSSEmaskData%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(AMSR2OSSEmaskData%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(AMSR2OSSEmaskData%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            AMSR2OSSEmaskData%n11, &
            AMSR2OSSEmaskData%n12, AMSR2OSSEmaskData%n21, &
            AMSR2OSSEmaskData%n22, AMSR2OSSEmaskData%w11, &
            AMSR2OSSEmaskData%w12, AMSR2OSSEmaskData%w21, &
            AMSR2OSSEmaskData%w22)

    else

       allocate(AMSR2OSSEmaskData%n11(&
            AMSR2OSSEmaskData%nc*&
            AMSR2OSSEmaskData%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            AMSR2OSSEmaskData%nc*AMSR2OSSEmaskData%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            AMSR2OSSEmaskData%n11)
       
    endif
    
  end subroutine AMSR2OSSEmask_init
  
end module AMSR2OSSEmask_Mod
