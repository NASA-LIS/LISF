!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISlsmTEFF_obsMod
!BOP
! 
! !MODULE: LISlsmTEFF_obsMod
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
!  01 Nov 2021    Yonghwan Kwon: Modified for effective soil temperature
!

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISlsmTEFF_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lsmteffobs
!
!EOP
  
  type, public :: lsmteffobsdec
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
     ! EMK
     integer :: num_ens
     integer :: num_tiles
     integer :: ntiles_pergrid
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

  end type lsmteffobsdec

  type(lsmteffobsdec)  :: lsmteffobs

contains

!BOP
! !ROUTINE: LISlsmTEFF_obsInit
! \label{LISlsmTEFF_obsInit}
! 
! !INTERFACE: 
  subroutine LISlsmTEFF_obsInit()
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

    lsmteffobs%run_dd             = LDT_rc%udef

    lsmteffobs%security_class     = ''
    lsmteffobs%distribution_class = ''
    lsmteffobs%data_category      = ''
    lsmteffobs%area_of_data       = ''
    lsmteffobs%write_interval     = ''

    call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%format, &
         label="LIS soil temperature output format:",rc=rc)
    call LDT_verify(rc,'LIS soil temperature output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%wopt, &
         label="LIS soil temperature output methodology:",rc=rc)
    call LDT_verify(rc,'LIS soil temperature output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%wstyle, &
         label="LIS soil temperature output naming style:",rc=rc)
    call LDT_verify(rc,'LIS soil temperature output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%map_proj, &
         label="LIS soil temperature output map projection:",rc=rc)
    call LDT_verify(rc,'LIS soil temperature output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%nest, &
         label="LIS soil temperature output nest index:",rc=rc)
    call LDT_verify(rc,'LIS soil temperature output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%odir, &
         label="LIS soil temperature output directory:",rc=rc)
    call LDT_verify(rc,'LIS soil temperature output directory: not defined')

    ! EMK Add LIS ensemble information
    call ESMF_ConfigGetAttribute(LDT_config, lsmteffobs%num_ens, &
         label="LIS ensemble size:", rc=rc)
    call LDT_verify(rc,'LIS ensemble size: not defined')
    if (lsmteffobs%num_ens < 1) then
       write(LDT_logunit,*)'[ERR] LIS ensemble size must be at least 1!'
       write(LDT_logunit,*)'[ERR] Read in ', lsmteffobs%num_ens
       call LDT_endrun()
    end if

    call ESMF_ConfigGetAttribute(LDT_config, lsmteffobs%num_tiles, &
         label="LIS total number of tiles (including ensembles):", rc=rc)
    call LDT_verify(rc,'LIS total number of tiles (including ensembles): ' // &
         'not defined')
    if (lsmteffobs%num_tiles < 1) then
       write(LDT_logunit,*) &
            '[ERR] LIS total number of tiles (including ensembles) must be'  &
            //'at least 1!'
       write(LDT_logunit,*)'[ERR] Read in ', lsmteffobs%num_tiles
       call LDT_endrun()
    end if

    call ESMF_ConfigGetAttribute(LDT_config, lsmteffobs%ntiles_pergrid, &
         label="LIS number of tiles per grid point:", rc=rc)
    call LDT_verify(rc,'LIS number of tiles per grid point: not defined')
    if (lsmteffobs%num_tiles < 1) then
       write(LDT_logunit,*) &
            '[ERR] LIS number of tiles per grid point must be at least 1!'
       write(LDT_logunit,*)'[ERR] Read in ', lsmteffobs%ntiles_pergrid
       call LDT_endrun()
    end if

    ! WMO-convention specific identifiers
    if ( lsmteffobs%wstyle == "WMO convention") then 
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%security_class, &
       label="LIS soil temperature security class:",rc=rc)
       call LDT_verify(rc,'LIS soil temperature security class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%distribution_class, &
       label="LIS soil temperature distribution class:",rc=rc)
       call LDT_verify(rc,'LIS soil temperature distribution class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%data_category, &
       label="LIS soil temperature data category:",rc=rc)
       call LDT_verify(rc,'LIS soil temperature data category: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%area_of_data, &
       label="LIS soil temperature area of data:",rc=rc)
       call LDT_verify(rc,'LIS soil temperature area of data: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%write_interval, &
       label="LIS soil temperature write interval:",rc=rc)
       call LDT_verify(rc,'LIS soil temperature write interval: not defined')
    endif

    if(lsmteffobs%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(6),rc=rc)       
       
       lsmteffobs%datares = min(lsmteffobs%run_dd(5),lsmteffobs%run_dd(6))

       lsmteffobs%nc    = (nint((lsmteffobs%run_dd(4)-lsmteffobs%run_dd(2))/&
            lsmteffobs%run_dd(5))) + 1
       lsmteffobs%nr    = (nint((lsmteffobs%run_dd(3)-lsmteffobs%run_dd(1))/&
            lsmteffobs%run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = lsmteffobs%nc
       gridDesci(3) = lsmteffobs%nr
       gridDesci(4) = lsmteffobs%run_dd(1)
       gridDesci(5) = lsmteffobs%run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = lsmteffobs%run_dd(3)
       gridDesci(8) = lsmteffobs%run_dd(4)
       gridDesci(9) = lsmteffobs%run_dd(5)
       gridDesci(10) = lsmteffobs%run_dd(6)
       gridDesci(20) = 64

    elseif(lsmteffobs%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS soil temperature domain y-dimension size: not defined')

       lsmteffobs%datares = lsmteffobs%run_dd(6)/100.0

       lsmteffobs%nc    = lsmteffobs%run_dd(7)
       lsmteffobs%nr    = lsmteffobs%run_dd(8)

       gridDesci = 0
       gridDesci(1) = 3 
       gridDesci(2) = lsmteffobs%nc
       gridDesci(3) = lsmteffobs%nr
       gridDesci(4) = lsmteffobs%run_dd(1)
       gridDesci(5) = lsmteffobs%run_dd(2)
       gridDesci(6) = 8
       gridDesci(7) = lsmteffobs%run_dd(4)
       gridDesci(8) = lsmteffobs%run_dd(6)
       gridDesci(9) = lsmteffobs%run_dd(6)
       gridDesci(10) = lsmteffobs%run_dd(3)
       gridDesci(11) = lsmteffobs%run_dd(5)
       gridDesci(20) = 64

    elseif(lsmteffobs%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(1),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS soil temperature domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,lsmteffobs%run_dd(8),rc=rc)

    endif

!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,lsmteffobs%datares)) then 

       allocate(lsmteffobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmteffobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmteffobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmteffobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(lsmteffobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmteffobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmteffobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(lsmteffobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            lsmteffobs%n11, &
            lsmteffobs%n12, lsmteffobs%n21, &
            lsmteffobs%n22, lsmteffobs%w11, &
            lsmteffobs%w12, lsmteffobs%w21, &
            lsmteffobs%w22)

    else

       allocate(lsmteffobs%n11(&
            lsmteffobs%nc*&
            lsmteffobs%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            lsmteffobs%nc*lsmteffobs%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            lsmteffobs%n11)
       
    endif
    
!  which variable we want in the DA obs computations. 
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%teff_obs, "K",1,1)
    LDT_DAobsData(1)%teff_obs%selectStats = 1

  end subroutine LISlsmTEFF_obsInit
  
end module LISlsmTEFF_obsMod
