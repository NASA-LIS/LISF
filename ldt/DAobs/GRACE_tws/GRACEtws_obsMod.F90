!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: GRACEtws_obsMod
! 
! !DESCRIPTION: 
!  This module provides the observation plugin for the 
!  GRACE terrestrial storage observations. 
!   
! !REVISION HISTORY: 
!  14 Mar 2013: Sujay Kumar, Initial Specification
!
!  21 DEc 2017: Augusto Getirana, Account for surface water storage (SWS) by subtracting the SWS signal from GRACE TWS. 
!  This process results in the land water storage (LWS), which is consistent with the LSM water storage. 
!  Getirana, A., et al., 2017b. Rivers and floodplains as key components of global terrestrial water storage variability. 
!  Geophysical Research Letters, 44. DOI: 10.1002/2017GL074684.
!
module GRACEtws_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACEtws_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GRACEtwsobs
!EOP
  type, public :: gracetwsobsdec
     integer :: nest    
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character(len=LDT_CONST_PATH_LEN) :: odir
     character*100 :: datasource 


     integer       :: reftime
     integer       :: tdims
     character(len=LDT_CONST_PATH_LEN) :: gracefile, gracescalefile
     character(len=LDT_CONST_PATH_LEN) :: graceerrfile
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

     !ag (21 DEc 2017)
     !Surface water storage variables
     integer              :: swsflag
     integer, allocatable :: swsdate(:)
     real,    allocatable :: swsavg(:,:,:)
     real,    allocatable :: nswsavg(:,:,:)

  end type gracetwsobsdec

  type(gracetwsobsdec) :: GRACEtwsobs

contains
  
!BOP
! 
! !ROUTINE: GRACEtws_obsInit
! \label{GRACEtws_obsInit}
! 
! !INTERFACE: 
  subroutine GRACEtws_obsinit()
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
    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%gracefile, &
         label="GRACE raw data filename:",rc=rc)
    call LDT_verify(rc,'GRACE raw data filename: not defined')

    ! Call to read scale factor file (BZ):
    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%gracescalefile, &
         label="GRACE scale factor filename:",rc=rc)
    call LDT_verify(rc,'GRACE scale factor filename: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%graceerrfile, &
         label="GRACE measurement error filename:",rc=rc)
    call LDT_verify(rc,'GRACE measurement error filename: not defined')

    ! LIS TWS output directory entries:
    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%odir, &
         label="LIS TWS output directory:",rc=rc)
    call LDT_verify(rc,'LIS TWS output directory: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%format, &
         label="LIS TWS output format:",rc=rc)
    call LDT_verify(rc,'LIS TWS output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%wopt, &
         label="LIS TWS output methodology:",rc=rc)
    call LDT_verify(rc,'LIS TWS output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%wstyle, &
         label="LIS TWS output naming style:",rc=rc)
    call LDT_verify(rc,'LIS TWS output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%map_proj, &
         label="LIS TWS output map projection:",rc=rc)
    call LDT_verify(rc,'LIS TWS output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%nest, &
         label="LIS TWS output nest index:",rc=rc)
    call LDT_verify(rc,'LIS TWS output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%b_syr, &
         label="GRACE baseline starting year:",rc=rc)
    call LDT_verify(rc,'GRACE baseline start year: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%b_eyr, &
         label="GRACE baseline ending year:",rc=rc)
    call LDT_verify(rc,'GRACE baseline ending year: not defined')

    !ag (21 Dec 2017)
    !Options accounting for surface water storage 
    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%swsflag, &
         label="LIS SWS output processing:",rc=rc)
    call LDT_verify(rc,'LIS SWS output processing: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%datasource, &
         label="GRACE data source:",rc=rc)
    if(rc.ne.0) then 
       write(LDT_logunit,*) '[ERR] GRACE data source: not defined'
       write(LDT_logunit,*) '[ERR] The options are ..'
       write(LDT_logunit,*) "[ERR] 'GRACE TWS Mascon 0.5 deg' or "
       write(LDT_logunit,*) "[ERR] 'GRACE TWS Original 1 deg'"
       call LDT_endrun()
    endif
    ! Read in LIS output - latlon grid description array entries:
    if(GRACEtwsobs%map_proj.eq."latlon") then 

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
       
       GRACEtwsobs%datares = min(run_dd(5),run_dd(6))

       GRACEtwsobs%nc  = (nint((run_dd(4)-run_dd(2))/&
            run_dd(5))) + 1
       GRACEtwsobs%nr  = (nint((run_dd(3)-run_dd(1))/&
            run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = GRACEtwsobs%nc
       gridDesci(3) = GRACEtwsobs%nr
       gridDesci(4) = run_dd(1)
       gridDesci(5) = run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = run_dd(3)
       gridDesci(8) = run_dd(4)
       gridDesci(9) = run_dd(5)
       gridDesci(10) = run_dd(6)
       gridDesci(20) = 64

    ! Read in LIS output - Lambert grid description array entries:
    elseif(GRACEtwsobs%map_proj.eq."lambert") then 

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
    elseif(GRACEtwsobs%map_proj.eq."polar") then 
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
    GRACEtwsobs%process_basin_scale = 0 

    call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%process_basin_scale, &
         label="GRACE process basin averaged observations:",rc=rc)
    call LDT_verify(rc,&
         'GRACE process basin averaged observations: not defined')

    if(GRACEtwsobs%process_basin_scale.eq.1) then 
       call ESMF_ConfigGetAttribute(LDT_config,GRACEtwsobs%basinmapfile, &
            label="GRACE basin map file:",rc=rc)
       call LDT_verify(rc,&
            'GRACE basin map file: not defined')
!------------------------------------------------------------------------------
! It is assumed that the basin categories are consecutively ordered from 1. 
!------------------------------------------------------------------------------

       allocate(GRACEtwsobs%basin_cat(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       ftn = LDT_getNextUnitNumber()
       open(ftn,file=GRACEtwsobs%basinmapfile, &
            form='unformatted',access='direct',&
            recl=LDT_rc%lnc(n)*LDT_rc%lnr(n)*4)
       read(ftn,rec=1) GRACEtwsobs%basin_cat
       call LDT_releaseUnitNumber(ftn)

       GRACEtwsobs%basin_cat_max = -1

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             if(GRACEtwsobs%basin_cat(c,r).ne.LDT_rc%udef) then 
                if(GRACEtwsobs%basin_cat(c,r).gt.GRACEtwsobs%basin_cat_max) then 
                   GRACEtwsobs%basin_cat_max = nint(GRACEtwsobs%basin_cat(c,r))
                endif
             endif
          enddo
       enddo

    endif  ! End basin-average processing mode

    GRACEtwsobs%startMode = .true. 

    ! Selecting which variable wanted in the DA obs computations:
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%TWS_obs, "mm",1,1)
    LDT_DAobsData(1)%tws_obs%selectStats = 1

    call LDT_get_julhr(2002,1,1,0,0,0,GRACEtwsobs%reftime)
    
    ! Check if GRACE Obs grid finer than the LIS output TWS files::
    if(LDT_isLDTatAfinerResolution(n,GRACEtwsobs%datares)) then 

       allocate(GRACEtwsobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEtwsobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEtwsobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEtwsobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       allocate(GRACEtwsobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEtwsobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEtwsobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(GRACEtwsobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       ! Downscale LIS grid to GRACE TWS obs:
       call bilinear_interp_input(n, gridDesci, &
            GRACEtwsobs%n11, &
            GRACEtwsobs%n12, GRACEtwsobs%n21, &
            GRACEtwsobs%n22, GRACEtwsobs%w11, &
            GRACEtwsobs%w12, GRACEtwsobs%w21, &
            GRACEtwsobs%w22)

    ! GRACE res. > LIS file res.:
    else
       allocate(GRACEtwsobs%n11(&
            GRACEtwsobs%nc*&
            GRACEtwsobs%nr))

       ! Average LIS file output to GRACE resolution:
       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            GRACEtwsobs%nc*GRACEtwsobs%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            GRACEtwsobs%n11)
    endif
!Bailing: These are GRACE product descriptions which are diff 
! currently, I use LDT_rc%gridDesc(n,9) to determine if it is 
   if (LDT_rc%gridDesc(n,9).lt.0.5)then  !0.25 GRACE                     
    gridDesci = 0                                            
    gridDesci(1) = 0                                         
    gridDesci(2) = 1440                                       
    gridDesci(3) = 720                                      
    gridDesci(4) = -89.875                                  
    gridDesci(5) = -179.875                                
    gridDesci(6) = 128                                   
    gridDesci(7) = 89.875                                
    gridDesci(8) = 179.875                              
    gridDesci(9) = 0.25                                
    gridDesci(10) = 0.25                              
    gridDesci(20) = 64                              
   elseif (LDT_rc%gridDesc(n,9).lt.1.0 .and. LDT_rc%gridDesc(n,9).gt.0.25)then  !0.5 GRACE                     
    gridDesci = 0                                            
    gridDesci(1) = 0                                         
    gridDesci(2) = 720                                       
    gridDesci(3) = 360                                      
    gridDesci(4) = -89.75                                  
    gridDesci(5) = -179.75                                
    gridDesci(6) = 128                                   
    gridDesci(7) = 89.75                                
    gridDesci(8) = 179.75                              
    gridDesci(9) = 0.5                                
    gridDesci(10) = 0.5                              
    gridDesci(20) = 64                              
   else                                            ! 1.0 GRACE
    gridDesci = 0                                 
    gridDesci(1) = 0                             
    gridDesci(2) = 360                          
    gridDesci(3) = 180                         
    gridDesci(4) = -89.5                      
    gridDesci(5) = -179.5                    
    gridDesci(6) = 128                      
    gridDesci(7) = 89.5                    
    gridDesci(8) = 179.5                  
    gridDesci(9) =  1.0                  
    gridDesci(10) = 1.0                 
    gridDesci(20) = 64                 
   end if 


    ! Define ASCAT Obs grid domain and resolution
    ! **NOTE: ASCAT domain longs. are actually 0 to 360 deg
    !    and not from -180W to 180E deg. This shift is accounted
    !    for in the reading and assingment of the values. 
  
    !  Latest GRACE TWS Mascon 0.5 deg dataset (B. Li):
!    elseif( LDT_rc%gridDesc(n,9) == 0.5 ) then
    if(GRACEtwsobs%datasource.eq."GRACE TWS Mascon 0.25 deg") then 
       GRACEtwsobs%gracenc = 1440
       GRACEtwsobs%gracenr = 720
       
       gridDesci = 0
       gridDesci(1) = 0
       gridDesci(2) = GRACEtwsobs%gracenc
       gridDesci(3) = GRACEtwsobs%gracenr
       gridDesci(4) = -89.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
    ! -- Original GRACE TWS 1.0 deg datasets:
    elseif(GRACEtwsobs%datasource.eq."GRACE TWS Mascon 0.5 deg") then 
       GRACEtwsobs%gracenc = 720
       GRACEtwsobs%gracenr = 360
       
       gridDesci = 0
       gridDesci(1) = 0
       gridDesci(2) = GRACEtwsobs%gracenc
       gridDesci(3) = GRACEtwsobs%gracenr
       gridDesci(4) = -89.75
       gridDesci(5) = -179.75
       gridDesci(6) = 128
       gridDesci(7) = 89.75
       gridDesci(8) = 179.75
       gridDesci(9) = 0.5
       gridDesci(10) = 0.5
       gridDesci(20) = 64
    elseif(GRACEtwsobs%datasource.eq."GRACE TWS Original 1 deg") then 
       GRACEtwsobs%gracenc = 360
       GRACEtwsobs%gracenr = 180

       gridDesci = 0
       gridDesci(1) = 0
       gridDesci(2) = GRACEtwsobs%gracenc
       gridDesci(3) = GRACEtwsobs%gracenr
       gridDesci(4) = -89.5
       gridDesci(5) = -179.5
       gridDesci(6) = 128
       gridDesci(7) = 89.5
       gridDesci(8) = 179.5
       gridDesci(9) = 1.0
       gridDesci(10) = 1.0
       gridDesci(20) = 64
    endif

    ! Use nearest neighbor to have product on GRACE product grid:
    allocate(GRACEtwsobs%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    call neighbor_interp_input(n, gridDesci,&
         GRACEtwsobs%n111 )

    allocate(GRACEtwsobs%lisavg(LDT_rc%lnc(n),LDT_rc%lnr(n)))
    allocate(GRACEtwsobs%nlisavg(LDT_rc%lnc(n),LDT_rc%lnr(n)))

    GRACEtwsobs%lisavg = 0 
    GRACEtwsobs%nlisavg = 0 


    call system('mkdir -p '//(LDT_rc%odir))

  end subroutine GRACEtws_obsinit
     
end module GRACEtws_obsMod
