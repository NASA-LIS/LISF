!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: GRACEQLtws_obsMod
! 
! !DESCRIPTION: 
!  This module provides the observation plugin for the 
!  GRACE terrestrial storage observations. 
!   
! !REVISION HISTORY: 
!  14 Feb 2018: Sujay Kumar, Initial Specification
!
!
module GRACEQLtws_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACEQLtws_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACEQLtwsobs
!EOP
  type, public :: gracetwsobsdec
     integer :: nest    
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character(len=LDT_CONST_PATH_LEN) :: odir
     character(len=LDT_CONST_PATH_LEN) :: datasource 

     integer       :: yr
     integer       :: reftime
     integer       :: tdims
     character(len=LDT_CONST_PATH_LEN) :: gracefile
     logical       :: startMode
     integer       :: gracenc, gracenr
     real, allocatable :: tvals(:)
     real, allocatable :: time_bounds(:,:)
     real, allocatable :: lwe_thickness(:,:,:)
     real, allocatable :: scalefactor(:,:)  !BZ
     real, allocatable :: merror(:,:)  
     real, allocatable :: lerror(:,:)  
     real, allocatable :: twsavg(:,:)
     real, allocatable :: lisavg(:,:)
     integer, allocatable :: nlisavg(:,:)

     integer           :: b_syr
     integer           :: b_eyr
     
     integer           :: nc,nr
     real              :: datares

     integer, allocatable :: n111(:)

     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)

     real, allocatable    :: w11(:)
     real, allocatable    :: w12(:)
     real, allocatable    :: w21(:)
     real, allocatable    :: w22(:)

     integer              :: process_basin_scale
     integer              :: basin_cat_max
     character(len=LDT_CONST_PATH_LEN)        :: basinmapfile
     real, allocatable    :: basin_cat(:,:)
     type(ESMF_Time)      :: startTime

  end type gracetwsobsdec

  type(gracetwsobsdec) :: GRACEQLtwsobs

contains
  
!BOP
! 
! !ROUTINE: GRACEQLtws_obsInit
! \label{GRACEQLtws_obsInit}
! 
! !INTERFACE: 
  subroutine GRACEQLtws_obsinit()
! !USES: 
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures 
!   required GRACE for terrestrial water storage (TWS) data. 
! 
!EOP
    real          :: run_dd(8)
    integer       :: rc
    integer       :: n 
    integer       :: ftn,c,r
    real          :: gridDesci(20)

    n = 1

    ! Config file entries:
    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%gracefile, &
         label="GRACE QL data directory:",rc=rc)
    call LDT_verify(rc,'GRACE raw data filename: not defined')

    ! LIS TWS output directory entries:
    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%odir, &
         label="LIS TWS output directory:",rc=rc)
    call LDT_verify(rc,'LIS TWS output directory: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%format, &
         label="LIS TWS output format:",rc=rc)
    call LDT_verify(rc,'LIS TWS output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%wopt, &
         label="LIS TWS output methodology:",rc=rc)
    call LDT_verify(rc,'LIS TWS output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%wstyle, &
         label="LIS TWS output naming style:",rc=rc)
    call LDT_verify(rc,'LIS TWS output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%map_proj, &
         label="LIS TWS output map projection:",rc=rc)
    call LDT_verify(rc,'LIS TWS output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%nest, &
         label="LIS TWS output nest index:",rc=rc)
    call LDT_verify(rc,'LIS TWS output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%b_syr, &
         label="GRACE baseline starting year:",rc=rc)
    call LDT_verify(rc,'GRACE baseline start year: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%b_eyr, &
         label="GRACE baseline ending year:",rc=rc)
    call LDT_verify(rc,'GRACE baseline ending year: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%datasource, &
         label="GRACE data source:",rc=rc)
    if(rc.ne.0) then 
       write(LDT_logunit,*) '[ERR] GRACE data source: not defined'
       write(LDT_logunit,*) '[ERR] The options are ..'
       write(LDT_logunit,*) "[ERR] 'GRACE TWS Mascon 0.5 deg' or "
       write(LDT_logunit,*) "[ERR] 'GRACE TWS Original 1 deg'"
       call LDT_endrun()
    endif
    ! Read in LIS output - latlon grid description array entries:
    if(GRACEQLtwsobs%map_proj.eq."latlon") then 

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
       
       GRACEQLtwsobs%datares = min(run_dd(5),run_dd(6))

       GRACEQLtwsobs%nc  = (nint((run_dd(4)-run_dd(2))/&
            run_dd(5))) + 1
       GRACEQLtwsobs%nr  = (nint((run_dd(3)-run_dd(1))/&
            run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = GRACEQLtwsobs%nc
       gridDesci(3) = GRACEQLtwsobs%nr
       gridDesci(4) = run_dd(1)
       gridDesci(5) = run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = run_dd(3)
       gridDesci(8) = run_dd(4)
       gridDesci(9) = run_dd(5)
       gridDesci(10) = run_dd(6)
       gridDesci(20) = 64

    ! Read in LIS output - Lambert grid description array entries:
    elseif(GRACEQLtwsobs%map_proj.eq."lambert") then 

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
       
    ! Read in LIS output - Polar grid description array entries:
    elseif(GRACEQLtwsobs%map_proj.eq."polar") then 
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

    ! Read Config entries if processing GRACE obs at basin average:
    GRACEQLtwsobs%process_basin_scale = 0 

    call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%process_basin_scale, &
         label="GRACE process basin averaged observations:",rc=rc)
    call LDT_verify(rc,&
         'GRACE process basin averaged observations: not defined')

    if(GRACEQLtwsobs%process_basin_scale.eq.1) then 
       call ESMF_ConfigGetAttribute(LDT_config,GRACEQLtwsobs%basinmapfile, &
            label="GRACE basin map file:",rc=rc)
       call LDT_verify(rc,&
            'GRACE basin map file: not defined')
!------------------------------------------------------------------------------
! It is assumed that the basin categories are consecutively ordered from 1. 
!------------------------------------------------------------------------------

       allocate(GRACEQLtwsobs%basin_cat(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       ftn = LDT_getNextUnitNumber()
       open(ftn,file=GRACEQLtwsobs%basinmapfile, &
            form='unformatted',access='direct',&
            recl=LDT_rc%lnc(n)*LDT_rc%lnr(n)*4)
       read(ftn,rec=1) GRACEQLtwsobs%basin_cat
       call LDT_releaseUnitNumber(ftn)

       GRACEQLtwsobs%basin_cat_max = -1

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(GRACEQLtwsobs%basin_cat(c,r).ne.LDT_rc%udef) then 
                if(GRACEQLtwsobs%basin_cat(c,r).gt.GRACEQLtwsobs%basin_cat_max) then 
                   GRACEQLtwsobs%basin_cat_max = nint(GRACEQLtwsobs%basin_cat(c,r))
                endif
             endif
          enddo
       enddo

    endif  ! End basin-average processing mode

    GRACEQLtwsobs%startMode = .true. 
    GRACEQLtwsobs%yr = -1

    ! Selecting which variable wanted in the DA obs computations:
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%TWS_obs, "mm",1,1)
    LDT_DAobsData(1)%tws_obs%selectStats = 1

    call LDT_get_julhr(2002,1,1,0,0,0,GRACEQLtwsobs%reftime)
    
    ! Check if GRACE Obs grid finer than the LIS output TWS files::
    if(LDT_isLDTatAfinerResolution(n,GRACEQLtwsobs%datares)) then 

       allocate(GRACEQLtwsobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEQLtwsobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEQLtwsobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEQLtwsobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       allocate(GRACEQLtwsobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEQLtwsobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEQLtwsobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEQLtwsobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       ! Downscale LIS grid to GRACE TWS obs:
       call bilinear_interp_input(n, gridDesci, &
            GRACEQLtwsobs%n11, &
            GRACEQLtwsobs%n12, GRACEQLtwsobs%n21, &
            GRACEQLtwsobs%n22, GRACEQLtwsobs%w11, &
            GRACEQLtwsobs%w12, GRACEQLtwsobs%w21, &
            GRACEQLtwsobs%w22)

    ! GRACE res. > LIS file res.:
    else
       allocate(GRACEQLtwsobs%n11(&
            GRACEQLtwsobs%nc*&
            GRACEQLtwsobs%nr))

       ! Average LIS file output to GRACE resolution:
       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            GRACEQLtwsobs%nc*GRACEQLtwsobs%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            GRACEQLtwsobs%n11)
    endif

    GRACEQLtwsobs%gracenc = 720
    GRACEQLtwsobs%gracenr = 360
    
    gridDesci = 0
    gridDesci(1) = 0
    gridDesci(2) = GRACEQLtwsobs%gracenc
    gridDesci(3) = GRACEQLtwsobs%gracenr
    gridDesci(4) = -89.75
    gridDesci(5) = -179.75
    gridDesci(6) = 128
    gridDesci(7) = 89.75
    gridDesci(8) = 179.75
    gridDesci(9) = 0.5
    gridDesci(10) = 0.5
    gridDesci(20) = 64

    ! Use nearest neighbor to have product on GRACE product grid:
    allocate(GRACEQLtwsobs%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    call neighbor_interp_input(n, gridDesci,&
         GRACEQLtwsobs%n111 )

    allocate(GRACEQLtwsobs%lisavg(LDT_rc%lnc(n),LDT_rc%lnr(n)))
    allocate(GRACEQLtwsobs%nlisavg(LDT_rc%lnc(n),LDT_rc%lnr(n)))

    GRACEQLtwsobs%lisavg = 0 
    GRACEQLtwsobs%nlisavg = 0 


    call system('mkdir -p '//(LDT_rc%odir))

    call ESMF_TimeSet(GRACEQLtwsobs%startTime, yy=1858,&
         mm = 11, dd=17, h = 0,&
         m = 0, s = 0, calendar = LDT_calendar,&
         rc=rc)
    call LDT_verify(rc, 'ESMF_TimeSet failed in GRACEQLtws_obsinit')
    

  end subroutine GRACEQLtws_obsinit
     
end module GRACEQLtws_obsMod
