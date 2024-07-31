!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module HYMAP2_routingMod
!BOP
! 
! !MODULE: HYMAP2_routingMod
! 
! !DESCRIPTION: 
!  This module provides the definition of data structures used to control
!  the operation of the HYMAP routing scheme. It also provides the entry
!  method for the initialization of HYMAP routing variables. 
! 
!  Reference: Augusto C.V. Getirana, A. Boone, D. Yamazaki, B. Decharme,
!             F. Papa, and N. Mognard, 2012: The Hydrological Modeling
!             and Analysis Platform (HyMAP): Evaluation in the Amazon
!             Basin.  Journal of Hydrometeorology, 13, 1641-1665.
!             doi: http://dx.doi.org/10.1175/JHM-D-12-021.1
!             Augusto C.V. Getirana, and Rodrigo C.D. Paiva, 2013: Mapping 
!             large-scale river flow hydraulics in the Amazon Basin,
!             Water Resour. Res., 49, doi:10.1002/wrcr.20212
!
! !REVISION HISTORY: 
! 08 Nov 2011: Augusto Getirana, Initial implementation in LIS based on the 
!                                HYMAP offline routing code. 
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia formulation, 
!                                adaptive time step and reservoir operation. 
! 13 Apr 2016: Augusto Getirana, Inclusion of option for hybrid runs with a 
!                                river flow map. 
!  7 Sep 2019: Augusto Getirana,  Added support for 2-way coupling
! 27 Apr 2020: Augusto Getirana,  Added support for urban drainage
! 
! !USES: 
  use ESMF
  use LIS_topoMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: HYMAP2_routingInit
  public :: HYMAP2_logunit ! file unit number used for diagnostic logging

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  !public :: LIS_evap_state
  public :: HYMAP2_routing_struc
  
  type, public :: HYMAP2_routing_dec
  type(ESMF_State), allocatable :: LIS_runoff_state(:)
     
  real             :: dt
  integer          :: useens
  integer          :: nz            !number of stages in the sub-grid discretization
! === Undefined Values ==================================
  integer          :: imis           ! undefined value
! === River sequence ====================================
  integer, allocatable :: seqx(:)       !1D sequence horizontal
  integer, allocatable :: seqy(:)       !1D sequence vertical
  integer, allocatable :: seqx_glb(:)
  integer, allocatable :: seqy_glb(:)

  integer              :: nseqall       !length of 1D sequnece for river and mouth  
! === Map ===============================================
  integer, allocatable :: sindex(:,:)     !2-D sequence index
  integer, allocatable :: outlet(:)       !outlet flag: 0 - river; 1 - ocean
  integer, allocatable :: outlet_glb(:)       !outlet flag: 0 - river; 1 - ocean
  integer, allocatable :: next(:)        !downstream grid cell
  integer, allocatable :: next_glb(:)        !downstream grid cell
  integer, allocatable :: nextx(:,:)        
  integer, allocatable :: nexty(:,:)        
  real,    allocatable :: elevtn(:)      !river bed elevation [m]
  real,    allocatable :: nxtdst(:)      !distance to the next grid [m]
  real,    allocatable :: grarea(:)      !area of the grid [m^2] 
  real,    allocatable :: fldgrd(:,:,:)    !floodplain gradient [-]
  real,    allocatable :: fldman(:)      !maning coefficient for floodplains [-]
  real,    allocatable :: fldhgt(:,:)    !floodplain height [m]
  real,    allocatable :: trnoff(:)      !runoff   concentration time parameter [day]
  real,    allocatable :: tbsflw(:)      !baseflow concentration time parameter [day]
  real,    allocatable :: cntime(:)      !concentration time [s]
  real,    allocatable :: fldstomax(:,:,:) !maximum floodplain storage [m3]
  real,    allocatable :: rivman(:,:)      !maning coefficient for rivers [-]
  real,    allocatable :: rivelv(:)      !elevation of river bed [m]
  real,    allocatable :: rivstomax(:,:)   !maximum river storage [m3]
  real,    allocatable :: rivlen(:)      !river length [m] 
  real,    allocatable :: rivwth(:,:)      !river width [m]
  real,    allocatable :: rivhgt(:,:)        !river heihgt [m]
  real                 :: rslpmin          !minimum slope [-]
  real,    allocatable :: rivare(:,:)        !river surface area [m2]
  real,    allocatable :: runoff0(:)       !input runoff [mm.dt-1]
  real,    allocatable :: basflw0(:)       !input baseflow [mm.dt-1]
  real,    allocatable :: flowmap(:)       !river flow type map [-]
  real,    allocatable :: bsfdwi_ratio(:) !deep water infiltration ratio from baseflow [-]
  real,    allocatable :: rnfdwi_ratio(:) !deep water infiltration ratio from surface runoff [-]
 ! === Outputs ==========================================
  real,    allocatable :: rivsto(:,:)    !river storage [m3]
  real,    allocatable :: rivdph(:,:)    !river depth [m]
  real,    allocatable :: rivvel(:,:)    !river flow velocity [m/s]
  real,    allocatable :: rivout(:,:)    !river outflow  [m3/s]
  real,    allocatable :: evpout(:,:)    !actual evaporation from open waters [m3]
  real,    allocatable :: fldout(:,:)    !floodplain outflow  [m3/s]
  real,    allocatable :: fldsto(:,:)    !floodplain storage   [m3]
  real,    allocatable :: flddph(:,:)    !floodplain depth [m]
  real,    allocatable :: fldvel(:,:)    !floodplain flow velocity [m/s] 
  real,    allocatable :: fldfrc(:,:)    !area fraction [-]
  real,    allocatable :: fldare(:,:)    !flooded area [m2]
  real,    allocatable :: sfcelv(:,:)    !water surface elevation [m]
  real,    allocatable :: bsfdwi(:,:)   !deep water infiltration from baseflow [mm.idt-1]
  real,    allocatable :: rnfdwi(:,:)   !deep water infiltration from surface runoff [mm.idt-1]
  real,    allocatable :: surfws(:,:)   !surface water storage [mm]
! === Linear reservoir components ================
  real,    allocatable :: rnfsto(:,:)    !runoff   reservoir storage [m3]
  real,    allocatable :: bsfsto(:,:)    !baseflow reservoir storage [m3]
  character(len=LIS_CONST_PATH_LEN) :: rstfile
  integer              :: numout
  integer              :: fileopen
  real                 :: outInterval 
  real                 :: rstInterval
  character*20         :: startMode
  integer              :: mo
     !ag (04Jun2012)
  integer              :: flowtype     !flow type flag
  integer              :: linresflag   !linear reservoir (time delay) flag
  integer              :: evapflag     !evaporation from open water flag

  !ag (03May2017)
  real,   allocatable  :: ewat(:,:)     !potential evaporation from open waters [km m-2 s-1]
  real,   allocatable  :: edif(:,:)     !differential evapotranspiration (evaporation from open waters - total evapotranspiration) used as input in HyMAP [km m-2 s-1]
  !integer              :: pevap_comp_method !Potential evaporation from open water internal computation method
  !character*200        :: pevap_source      !Potential evaporation readin source directory
  
  integer              :: dwiflag      !deep water infiltration flag
  !character*10         :: evapsrc
     !ag (11Sep2015)
  integer              :: steptype     !time step type flag

  character(len=LIS_CONST_PATH_LEN) :: LISdir     !if LIS output is being read
     
! === Local inertia variables =====================
  real,    allocatable  :: rivout_pre(:,:)
  real,    allocatable  :: rivdph_pre(:,:)
  real,    allocatable  :: fldout_pre(:,:)
  real,    allocatable  :: flddph_pre(:,:)
  real,    allocatable  :: fldelv1(:,:)
  real                  :: cadp
  real                  :: grv
  
  real,    allocatable  :: dtaout(:,:)
                                    
!ag (30Jan2016)
! === timer variables ====================
  real                 :: dt_proc

! === reservoir operation variables/parameters ===
  integer, allocatable :: resopaltloc(:)
  real,    allocatable :: resopoutmin(:)
  real*8,  allocatable :: tresopalt(:,:)
  real,    allocatable :: resopalt(:,:)
  integer              :: nresop        !number of reservoirs
  integer              :: ntresop       !time series length (number of time steps in the input files)
  integer              :: resopflag
  integer              :: resoptype
  character(len=LIS_CONST_PATH_LEN) :: resopdir
  character*100         :: resopheader
  !ag (29Jun2016)
  integer              :: floodflag
  character(len=LIS_CONST_PATH_LEN) :: HYMAP_dfile      

! === 2-way coupling variables/parameters ===
  real,   allocatable  :: rivstotmp(:,:)     !River Storage [m3]
  real,   allocatable  :: fldstotmp(:,:)     !Flood Storage [m3]
  real,   allocatable  :: fldfrctmp(:,:)     !Flooded Fraction [m3]
  integer                   :: enable2waycpl
  real                        :: fldfrc2waycpl

!ag(27Apr2020)
! === urban drainage and flood modeling ===
  real,   allocatable  :: drsto(:,:)      !urban drainage network water storage [m3]
  real,   allocatable  :: drout(:,:)      !urban drainage network outflow [m3/s]
  real,   allocatable  :: drstomax(:)     !urban drainage network water storage capacity [m3]
  real,   allocatable  :: droutlet(:)     !outlet flag: 0 - drainage network; 1 - river
  real,   allocatable  :: droutlet_glb(:)     !global outlet flag: 0 - drainage network; 1 - river
  real,   allocatable  :: drtotwth(:)     !sum of gutter width within HyMAP grid cell [m]
  real,   allocatable  :: drnoutlet(:)    ! average number of drainage outlets within a grid cell [-]
  real,   allocatable  :: drtotlgh(:)     ! total urban drainage network length within a grid cell [m]
  character(len=LIS_CONST_PATH_LEN) :: drfile          !urban drainage parametere file name
  !gutter parameters
  real                        :: drwth     ! gutter width [m]
  real                        :: drhgt      ! gutter height [m]
  real                        :: drden     ! gutter density [units/m2]
  real                        :: drvel      ! gutter water intake velocity [m/s]
  !drainage network parameters
  real                        :: drblk     ! street block length [m]
  real                        :: drrad     ! underground squared pipe width [m]
  real                        :: drlgh     ! drainage length density [m/m2]
  real                        :: drman   ! Roughness coefficient for ciment pipes [-]
  real                        :: drslp      ! drainage system slope [m/m]

!ag(8Aug2020)
! === direct streamflow insertion variables/parameters ===
  integer, allocatable :: insertloc(:)
  real*8,  allocatable :: tinsert(:,:)
  real,    allocatable :: insertdis(:,:)
  integer              :: ninsert        !number of gauges
  integer              :: ntinsert       !time series length (number of time steps in the input files)
  integer              :: insertflag
  character(len=LIS_CONST_PATH_LEN) :: insertdir
  character*100         :: insertheader

!ag(30Mar2021)
! === ocean tides variables/parameters ===
  integer, allocatable  :: tidesloc(:)
  real*8,  allocatable  :: ttides(:,:)
  real,    allocatable  :: tides(:,:)
  integer               :: ntides        !number of outlets
  integer               :: nttides       !time series length (number of time steps in the input files)
  integer               :: tidesflag
  character(len=LIS_CONST_PATH_LEN) :: tidesdir
  character*100         :: tidesheader

  end type HYMAP2_routing_dec

  type(HYMAP2_routing_dec), allocatable :: HYMAP2_routing_struc(:)
  
  integer :: HYMAP2_logunit ! file unit number used for diagnostic logging

contains
 
!BOP
!
! !ROUTINE: HYMAP2_routingInit
! \label{HYMAP2_routingInit}
! 
  subroutine HYMAP2_routingInit
    !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_routingMod    
    use HYMAP2_initMod
    use LIS_perturbMod
    use LIS_mpiMod
   !ag(27Apr2020)
    use HYMAP2_urbanMod
    
    integer              :: n 
    integer              :: ftn 
    integer              :: status
    integer              :: i,j,k
    type(ESMF_Grid)      :: global_grid
    type(ESMF_DistGrid)  :: global_grid_dg
    type(ESMF_ArraySpec) :: realarrspec
    type(ESMF_Field)     :: sf_runoff_field
    type(ESMF_Field)     :: baseflow_field
    type(ESMF_Field)     :: evapotranspiration_field
    !type(ESMF_Field)     :: potential_evaporation_field
    real, pointer        :: sfrunoff(:)
    real, pointer        :: baseflow(:)
    real, pointer        :: evapotranspiration(:)
    !real, pointer        :: potevap(:)
    character*100        :: ctitle
    integer              :: iseq,ix,iy,jx,jy
    character*10         :: time
    character*20         :: flowtype
    !ag (11Sep2015)
    character*20         :: steptype 
    !ag (31Jan2016)
    character*20         :: resoptype
    !ag (8Aaug2020)
    character*20         :: inserttype
        
    !ag (19Feb2016)
    real,    allocatable :: tmp_real(:,:),tmp_real_nz(:,:,:)
    integer, allocatable :: nextx(:,:),nexty(:,:),mask(:,:),maskg(:,:)
    real,    allocatable :: elevtn(:,:),uparea(:,:),basin(:,:)
    
    !ag (11Mar2016)
    character*100  :: temp1
    integer        :: rc
    character*1    :: fproc(4)
    integer        :: ios
    integer        :: final_dirpos
    character(len=LIS_CONST_PATH_LEN) :: diag_fname
    character(len=LIS_CONST_PATH_LEN) :: diag_dir
    integer, external  :: LIS_create_subdirs

    !ag (03Apr2017)
    character*20 :: evapflag0
    !character*20 :: pevap_comp_method
    integer       :: ic, c,r,ic_down
    integer       :: m
    integer       :: gdeltas

    !ag (12Sep2019)
    type(ESMF_Field)     :: rivsto_field
    type(ESMF_Field)     :: fldsto_field
    type(ESMF_Field)     :: fldfrc_field
    real, pointer        :: rivstotmp(:)
    real, pointer        :: fldstotmp(:)
    real, pointer        :: fldfrctmp(:)

    type(ESMF_DistGrid)  :: patchDG
    type(ESMF_DistGrid)  :: gridDG
    integer, allocatable :: deblklist(:,:,:)
    integer              :: stid,enid
    character*1   :: caseid(3)
    character*40, allocatable:: vname(:)
    character*40, allocatable:: pertobjs(:)
    integer     , allocatable:: order(:)
    real        ,allocatable :: ccorr(:,:)
    real        ,allocatable :: stmin(:)
    real        ,allocatable :: stmax(:)
    real        ,allocatable :: ssdev(:)
    type(ESMF_ArraySpec) :: arrspec1, arrspec2
    type(ESMF_Field)     :: varField
    type(ESMF_Field)     :: varIncrField
    type(ESMF_Field)     :: pertField
    type(pert_dec_type)  :: routing_pert
    character*100        :: temp
    character*1          :: nestid(2)
    integer              :: max_index
    logical              :: name_found
    character*20         :: alglist(10)
    logical              :: Routing_DAvalid
    
    allocate(HYMAP2_routing_struc(LIS_rc%nnest))
    
    HYMAP2_logunit=10001
    
    do n=1, LIS_rc%nnest
       HYMAP2_routing_struc(n)%rslpmin  = 1e-5  !minimum slope
       HYMAP2_routing_struc(n)%nz      = 10    !number of stages in the sub-grid discretization
       HYMAP2_routing_struc(n)%imis     = -9999 !undefined integer value
       !HYMAP2_routing_struc(n)%cadp     = 0.7   !alfa coefficient for adaptative time step as described in Bates et al., (2010) [-]
       !HYMAP2_routing_struc(n)%dstend   = 25000 !river length to the ocean [m]
       HYMAP2_routing_struc(n)%grv      = 9.81  !gravity accerelation [m/s2]
       HYMAP2_routing_struc(n)%numout   = 0 
       HYMAP2_routing_struc(n)%fileopen = 0 
       HYMAP2_routing_struc(n)%dt_proc  = 0.
    enddo
       
    !ag (12Sep2019)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 enable 2-way coupling:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP2_routing_struc(n)%enable2waycpl=0
       call ESMF_ConfigGetAttribute(LIS_config,HYMAP2_routing_struc(n)%enable2waycpl,&
            default=0, rc=status)

       if (HYMAP2_routing_struc(n)%enable2waycpl==1) then
         write(LIS_logunit,*) '[INFO] HYMAP2 2-way coupling: activated'
         call ESMF_ConfigFindLabel(LIS_config,&
              "HYMAP2 2-way coupling flooded fraction threshold:",rc=status)
         call ESMF_ConfigGetAttribute(LIS_config,HYMAP2_routing_struc(n)%fldfrc2waycpl,rc=status)
         call LIS_verify(status,&
              "HYMAP2 2-way coupling flooded fraction threshold: not defined")
       else
          write(LIS_logunit,*) '[INFO] HYMAP2 2-way coupling: deactivated'
       endif
    enddo


    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model time step:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model time step: not defined")

       call LIS_parseTimeString(time,HYMAP2_routing_struc(n)%dt)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model output interval: not defined")

       call LIS_parseTimeString(time,HYMAP2_routing_struc(n)%outInterval)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 run in ensemble mode:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP2_routing_struc(n)%useens,rc=status)
       call LIS_verify(status,&
            "HYMAP2 run in ensemble mode: not defined")
    enddo

    !ag (13Apr2016)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing method:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,& 
            flowtype,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing method: not defined")
       if(flowtype.eq."kinematic") then
          HYMAP2_routing_struc(n)%flowtype = 1
       elseif(flowtype.eq."diffusive") then 
          HYMAP2_routing_struc(n)%flowtype = 2
       elseif(flowtype.eq."local inertia") then 
          HYMAP2_routing_struc(n)%flowtype = 3
       elseif(flowtype.eq."hybrid") then 
          HYMAP2_routing_struc(n)%flowtype = 0
        !ag(27Apr2020)
       elseif(flowtype.eq."urban") then 
          HYMAP2_routing_struc(n)%flowtype = 4
       else
         write(LIS_logunit,*)"[ERR] HYMAP2 routing method: unknown value ",trim(flowtype)
         call LIS_endrun()
       endif
    enddo
    
    !ag (11Sep2015)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model time step method:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            steptype,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model time step method: not defined")
       if(steptype.eq."constant") then 
          HYMAP2_routing_struc(n)%steptype = 1
       elseif(steptype.eq."adaptive") then 
          HYMAP2_routing_struc(n)%steptype = 2
       else
         write(LIS_logunit,*)"[ERR] HYMAP2 routing model time step method: unknown value"
         call LIS_endrun()
       endif
    enddo

    !ag (29Jan2016)
    if(steptype.eq."adaptive")then 
      call ESMF_ConfigFindLabel(LIS_config,&
           "HYMAP2 routing model adaptive time step alfa coefficient:",rc=status)
      do n=1, LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,&
              HYMAP2_routing_struc(n)%cadp,rc=status)
         call LIS_verify(status,&
              "HYMAP2 routing model adaptive time step alfa coefficient: not defined")
      enddo
    endif
    
    !ag (29Jun2016)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 floodplain dynamics:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP2_routing_struc(n)%floodflag,rc=status)
       call LIS_verify(status,&
            "HYMAP2 floodplain dynamics: not defined")
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model linear reservoir flag:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP2_routing_struc(n)%linresflag,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model linear reservoir flag: not defined")
    enddo
    
    !ag (24Apr2017)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model evaporation option:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,& 
            evapflag0,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model evaporation option: not defined")
       if(evapflag0.eq."none") then
          HYMAP2_routing_struc(n)%evapflag = 0
       elseif(evapflag0.eq."penman") then 
          HYMAP2_routing_struc(n)%evapflag = 1
       else
         write(LIS_logunit,*)"[ERR] HYMAP2 routing model evaporation option: ",trim(evapflag0)
         call LIS_endrun()
       endif
    enddo 

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model dwi flag:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP2_routing_struc(n)%dwiflag,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model dwi flag: not defined")    
    enddo

    !ag (8Aug2020)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 discharge direct insertion:",rc=status)
    do n=1, LIS_rc%nnest
      HYMAP2_routing_struc(n)%insertflag=0
      call ESMF_ConfigGetAttribute(LIS_config,&
           HYMAP2_routing_struc(n)%insertflag,rc=status)
      !call LIS_verify(status,&
      !     "HYMAP2 discharge direct insertion: not defined")

      if(HYMAP2_routing_struc(n)%insertflag==1)then
        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 number of gauges:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%ninsert,rc=status)
        call LIS_verify(status,&
             "HYMAP2 number of gauges: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 discharge direct insertion time series length:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%ntinsert,rc=status)
        call LIS_verify(status,&
             "HYMAP2 discharge direct insertion time series length: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 discharge direct insertion input directory:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%insertdir,rc=status)
        call LIS_verify(status,&
             "HYMAP2 discharge direct insertion input directory: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
           "HYMAP2 discharge direct insertion header filename:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%insertheader,rc=status)
        call LIS_verify(status,&
             "HYMAP2 discharge direct insertion header filename: not defined")
      endif
    enddo

    !ag (30Mar2021)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 ocean tides:",rc=status)
    do n=1, LIS_rc%nnest
      HYMAP2_routing_struc(n)%tidesflag=0
      call ESMF_ConfigGetAttribute(LIS_config,&
           HYMAP2_routing_struc(n)%tidesflag,rc=status)
      if(HYMAP2_routing_struc(n)%tidesflag==1)then
        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 ocean tides number of outlets:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%ntides,rc=status)
        call LIS_verify(status,&
             "HYMAP2 ocean tides number of outlets: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 ocean tides time series length:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%nttides,rc=status)
        call LIS_verify(status,&
             "HYMAP2 ocean tides time series length: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 ocean tides input directory:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%tidesdir,rc=status)
        call LIS_verify(status,&
             "HYMAP2 ocean tides input directory: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
           "HYMAP2 ocean tides header filename:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%tidesheader,rc=status)
        call LIS_verify(status,&
             "HYMAP2 ocean tides header filename: not defined")
      endif
    enddo
 
    !ag (4Feb2016)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 reservoir operation option:",rc=status)
    do n=1, LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,&
           HYMAP2_routing_struc(n)%resopflag,default=0,rc=status)

      if(HYMAP2_routing_struc(n)%resopflag==1)then
        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 number of reservoirs:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%nresop,rc=status)
        call LIS_verify(status,&
             "HYMAP2 number of reservoirs: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 reservoir operation input time series length:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%ntresop,rc=status)
        call LIS_verify(status,&
             "HYMAP2 reservoir operation input time series length: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 reservoir operation input directory:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%resopdir,rc=status)
        call LIS_verify(status,&
             "HYMAP2 reservoir operation input directory: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
           "HYMAP2 reservoir operation header filename:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%resopheader,rc=status)
        call LIS_verify(status,&
             "HYMAP2 reservoir operation header filename: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 reservoir operation input data type:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
              resoptype,rc=status)
        call LIS_verify(status,&
             "HYMAP2 reservoir operation input data type: not defined")
        if(resoptype.eq."water level") then 
          HYMAP2_routing_struc(n)%resoptype = 1
        elseif(resoptype.eq."streamflow") then 
          HYMAP2_routing_struc(n)%resoptype = 2
          write(LIS_logunit,*)"[ERR] HYMAP2 reservoir operation input data type: unknown value"
          call LIS_endrun()
        endif
      endif
    enddo

    write(LIS_logunit,*) '[INFO] Initializing HYMAP....'
    !allocate matrixes
    do n=1, LIS_rc%nnest
       write(LIS_logunit,*)'[INFO] columns and rows', &
            LIS_rc%lnc(n),LIS_rc%lnr(n)
       
       !ag (19Feb2016)
       allocate(HYMAP2_routing_struc(n)%nextx(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       ctitle = 'HYMAP_flow_direction_x'
       call HYMAP2_read_param_int_2d_global(ctitle,n,HYMAP2_routing_struc(n)%nextx)

       allocate(HYMAP2_routing_struc(n)%nexty(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       ctitle = 'HYMAP_flow_direction_y'
       call HYMAP2_read_param_int_2d_global(ctitle,n,HYMAP2_routing_struc(n)%nexty)
       
       allocate(elevtn(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       ctitle = 'HYMAP_grid_elevation'
       call HYMAP2_read_param_real_2d(ctitle,n,elevtn)
       
       allocate(uparea(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       ctitle = 'HYMAP_drain_area'
       call HYMAP2_read_param_real_2d(ctitle,n,uparea)
       
       allocate(basin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       ctitle = 'HYMAP_basin'
       call HYMAP2_read_param_real_2d(ctitle,n,basin)  
       
       allocate(mask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       ctitle = 'HYMAP_basin_mask'
       call HYMAP2_read_param_int_2d(ctitle,n,mask)  

       allocate(maskg(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       ctitle = 'HYMAP_basin_mask'
       call HYMAP2_read_param_int_2d_global(ctitle,n,maskg)  

       !Assign the mask to the routing data structure
       LIS_routing(n)%dommask = mask
       LIS_routing(n)%nextx = HYMAP2_routing_struc(n)%nextx

       write(LIS_logunit,*) '[INFO] Get number of HYMAP2 grid cells'
       call HYMAP2_get_vector_size(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            LIS_rc%gnc(n),LIS_rc%gnr(n),&
            LIS_ews_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1), &
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%nextx,&
            mask,HYMAP2_routing_struc(n)%nseqall)

       LIS_rc%nroutinggrid(n) = HYMAP2_routing_struc(n)%nseqall

       gdeltas = HYMAP2_routing_struc(n)%nseqall
       
#if (defined SPMD)
       call MPI_ALLREDUCE(LIS_rc%nroutinggrid(n),&
            LIS_rc%glbnroutinggrid(n),1,&
            MPI_INTEGER,MPI_SUM,&
            LIS_mpi_comm,status)

       call MPI_ALLGATHER(gdeltas,1,MPI_INTEGER,&
            LIS_routing_gdeltas(n,:),1,MPI_INTEGER,&
            LIS_mpi_comm,status)

#else
       LIS_rc%glbnroutinggrid(n) = LIS_rc%nroutinggrid(n)
#endif

       if(LIS_masterproc) then 
          LIS_routing_goffsets(n,:) = 0 
          do i=1,LIS_npes-1
             LIS_routing_goffsets(n,i) = &
                  LIS_routing_goffsets(n,i-1) +& 
                  LIS_routing_gdeltas(n,i-1)
          enddo
       endif
#if (defined SPMD)
       call MPI_BCAST(LIS_routing_goffsets(n,:), &
            LIS_npes, MPI_INTEGER,0, &
            LIS_mpi_comm, status)
#endif
       allocate(HYMAP2_routing_struc(n)%seqx(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%seqy(HYMAP2_routing_struc(n)%nseqall))

       allocate(HYMAP2_routing_struc(n)%seqx_glb(LIS_rc%glbnroutinggrid(n)))
       allocate(HYMAP2_routing_struc(n)%seqy_glb(LIS_rc%glbnroutinggrid(n)))

       allocate(HYMAP2_routing_struc(n)%sindex(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(HYMAP2_routing_struc(n)%outlet(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%outlet_glb(LIS_rc%glbnroutinggrid(n)))
       allocate(HYMAP2_routing_struc(n)%next(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%next_glb(LIS_rc%glbnroutinggrid(n)))

       allocate(HYMAP2_routing_struc(n)%elevtn(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%nxtdst(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%grarea(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%fldman(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%fldhgt(HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%nz))
       allocate(HYMAP2_routing_struc(n)%cntime(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%trnoff(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%tbsflw(HYMAP2_routing_struc(n)%nseqall))
       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          allocate(HYMAP2_routing_struc(n)%rivman(HYMAP2_routing_struc(n)%nseqall,1))
          allocate(HYMAP2_routing_struc(n)%rivwth(HYMAP2_routing_struc(n)%nseqall,1))
          allocate(HYMAP2_routing_struc(n)%rivhgt(HYMAP2_routing_struc(n)%nseqall,1))
          allocate(HYMAP2_routing_struc(n)%rivstomax(HYMAP2_routing_struc(n)%nseqall,1))
          allocate(HYMAP2_routing_struc(n)%rivare(HYMAP2_routing_struc(n)%nseqall,1))
          allocate(HYMAP2_routing_struc(n)%fldstomax(HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%nz,1))
          allocate(HYMAP2_routing_struc(n)%fldgrd(HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%nz,1))

       else
          allocate(HYMAP2_routing_struc(n)%rivman(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivwth(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivhgt(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivstomax(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivare(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldstomax(HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%nz,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldgrd(HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%nz,LIS_rc%nensem(n)))
       endif
          
       allocate(HYMAP2_routing_struc(n)%rivelv(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rivlen(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%runoff0(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%basflw0(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%flowmap(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rnfdwi_ratio(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%bsfdwi_ratio(HYMAP2_routing_struc(n)%nseqall))
       !ag (4Feb2016)
       if(HYMAP2_routing_struc(n)%resopflag==1)then
          allocate(HYMAP2_routing_struc(n)%resopaltloc(HYMAP2_routing_struc(n)%nresop))
          allocate(HYMAP2_routing_struc(n)%resopoutmin(HYMAP2_routing_struc(n)%nresop))
          allocate(HYMAP2_routing_struc(n)%tresopalt(HYMAP2_routing_struc(n)%nresop,HYMAP2_routing_struc(n)%ntresop))
          allocate(HYMAP2_routing_struc(n)%resopalt(HYMAP2_routing_struc(n)%nresop,HYMAP2_routing_struc(n)%ntresop))
       endif
       
       !ag(27Apr2020)
!        if(HYMAP2_routing_struc(n)%flowtype ==4)then
          allocate(HYMAP2_routing_struc(n)%drstomax(HYMAP2_routing_struc(n)%nseqall))
          allocate(HYMAP2_routing_struc(n)%droutlet(HYMAP2_routing_struc(n)%nseqall))
          allocate(HYMAP2_routing_struc(n)%droutlet_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(HYMAP2_routing_struc(n)%drtotwth(HYMAP2_routing_struc(n)%nseqall))
          allocate(HYMAP2_routing_struc(n)%drnoutlet(HYMAP2_routing_struc(n)%nseqall))
          allocate(HYMAP2_routing_struc(n)%drtotlgh(HYMAP2_routing_struc(n)%nseqall))
!        endif
       
       !ag(8Aug2020)
       if(HYMAP2_routing_struc(n)%insertflag==1)then
          allocate(HYMAP2_routing_struc(n)%insertloc(HYMAP2_routing_struc(n)%ninsert))
          allocate(HYMAP2_routing_struc(n)%tinsert(HYMAP2_routing_struc(n)%ninsert,HYMAP2_routing_struc(n)%ntinsert))
          allocate(HYMAP2_routing_struc(n)%insertdis(HYMAP2_routing_struc(n)%ninsert,HYMAP2_routing_struc(n)%ntinsert))
       endif

       !ag(30Mar2021)
       if(HYMAP2_routing_struc(n)%tidesflag==1)then
          allocate(HYMAP2_routing_struc(n)%tidesloc(HYMAP2_routing_struc(n)%ntides))
          allocate(HYMAP2_routing_struc(n)%ttides(HYMAP2_routing_struc(n)%ntides,HYMAP2_routing_struc(n)%nttides))
          allocate(HYMAP2_routing_struc(n)%tides(HYMAP2_routing_struc(n)%ntides,HYMAP2_routing_struc(n)%nttides))
       endif

       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          allocate(HYMAP2_routing_struc(n)%rivsto(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%rivdph(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%rivvel(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%rivout(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%evpout(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldout(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldsto(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%flddph(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldvel(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldfrc(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldare(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%sfcelv(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%rnfsto(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%bsfsto(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%rnfdwi(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%bsfdwi(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%surfws(HYMAP2_routing_struc(n)%nseqall,&
               1))
          
          allocate(HYMAP2_routing_struc(n)%dtaout(HYMAP2_routing_struc(n)%nseqall,&
               1))
          !ag (19Jan2016)
          allocate(HYMAP2_routing_struc(n)%rivout_pre(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%rivdph_pre(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldout_pre(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%flddph_pre(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldelv1(HYMAP2_routing_struc(n)%nseqall,&
               1))         
          !ag (03May2017)
          allocate(HYMAP2_routing_struc(n)%ewat(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%edif(HYMAP2_routing_struc(n)%nseqall,&
               1))
          !ag (12Sep2019)
          allocate(HYMAP2_routing_struc(n)%rivstotmp(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldstotmp(HYMAP2_routing_struc(n)%nseqall,&
               1))
          allocate(HYMAP2_routing_struc(n)%fldfrctmp(HYMAP2_routing_struc(n)%nseqall,&
               1))
          !ag(27Apr2020)
!          if(HYMAP2_routing_struc(n)%flowtype ==4)then
             allocate(HYMAP2_routing_struc(n)%drsto(HYMAP2_routing_struc(n)%nseqall,1))
             allocate(HYMAP2_routing_struc(n)%drout(HYMAP2_routing_struc(n)%nseqall,1))
!          endif
       else
          allocate(HYMAP2_routing_struc(n)%rivsto(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivdph(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivvel(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivout(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%evpout(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldout(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldsto(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%flddph(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldvel(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldfrc(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldare(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%sfcelv(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rnfsto(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%bsfsto(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rnfdwi(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%bsfdwi(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%surfws(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))

          allocate(HYMAP2_routing_struc(n)%dtaout(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))

          !ag (19Jan2016)
          allocate(HYMAP2_routing_struc(n)%rivout_pre(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%rivdph_pre(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldout_pre(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%flddph_pre(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldelv1(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))

          !ag (03May2017)
          allocate(HYMAP2_routing_struc(n)%ewat(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%edif(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))

          !ag (12Sep2019)
          allocate(HYMAP2_routing_struc(n)%rivstotmp(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldstotmp(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          allocate(HYMAP2_routing_struc(n)%fldfrctmp(HYMAP2_routing_struc(n)%nseqall,&
               LIS_rc%nensem(n)))
          !ag(27Apr2020)
!          if(HYMAP2_routing_struc(n)%flowtype ==4)then
             allocate(HYMAP2_routing_struc(n)%drsto(HYMAP2_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
             allocate(HYMAP2_routing_struc(n)%drout(HYMAP2_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
!          endif
       endif

       HYMAP2_routing_struc(n)%seqx=0.0
       HYMAP2_routing_struc(n)%seqy=0.0

       HYMAP2_routing_struc(n)%sindex=0.0
       HYMAP2_routing_struc(n)%outlet=0.0
       HYMAP2_routing_struc(n)%outlet_glb=0.0
       HYMAP2_routing_struc(n)%next=0.0
       HYMAP2_routing_struc(n)%next_glb=0.0
       !HYMAP2_routing_struc(n)%nextx=0.0
       !HYMAP2_routing_struc(n)%nexty=0.0 
       !HYMAP2_routing_struc(n)%mask=0.0
       HYMAP2_routing_struc(n)%elevtn=0.0
       HYMAP2_routing_struc(n)%nxtdst=0.0
       HYMAP2_routing_struc(n)%grarea=0.0
       HYMAP2_routing_struc(n)%fldgrd=0.0
       HYMAP2_routing_struc(n)%fldman=0.0
       HYMAP2_routing_struc(n)%fldhgt=0.0
       HYMAP2_routing_struc(n)%cntime=0.0
       HYMAP2_routing_struc(n)%trnoff=0.0
       HYMAP2_routing_struc(n)%tbsflw=0.0
       HYMAP2_routing_struc(n)%fldstomax=0.0
       HYMAP2_routing_struc(n)%rivman=0.0
       HYMAP2_routing_struc(n)%rivelv=0.0
       HYMAP2_routing_struc(n)%rivstomax=0.0
       HYMAP2_routing_struc(n)%rivlen=0.0
       HYMAP2_routing_struc(n)%rivwth=0.0
       HYMAP2_routing_struc(n)%rivhgt=0.0
       HYMAP2_routing_struc(n)%rivare=0.0
       HYMAP2_routing_struc(n)%runoff0=0.0
       HYMAP2_routing_struc(n)%basflw0=0.0
       HYMAP2_routing_struc(n)%rivsto=0.0
       HYMAP2_routing_struc(n)%rivdph=0.0
       HYMAP2_routing_struc(n)%rivvel=0.0
       HYMAP2_routing_struc(n)%rivout=0.0
       HYMAP2_routing_struc(n)%evpout=0.0
       HYMAP2_routing_struc(n)%fldout=0.0
       HYMAP2_routing_struc(n)%fldsto=0.0
       HYMAP2_routing_struc(n)%flddph=0.0
       HYMAP2_routing_struc(n)%fldvel=0.0
       HYMAP2_routing_struc(n)%fldfrc=0.0
       HYMAP2_routing_struc(n)%fldare=0.0
       HYMAP2_routing_struc(n)%sfcelv=0.0
       HYMAP2_routing_struc(n)%rnfsto=0.0
       HYMAP2_routing_struc(n)%bsfsto=0.0
       HYMAP2_routing_struc(n)%flowmap=0.0
       HYMAP2_routing_struc(n)%rnfdwi_ratio=0.0
       HYMAP2_routing_struc(n)%bsfdwi_ratio=0.0
       HYMAP2_routing_struc(n)%rnfdwi=0.0
       HYMAP2_routing_struc(n)%bsfdwi=0.0
       HYMAP2_routing_struc(n)%surfws=0.0
       
       HYMAP2_routing_struc(n)%dtaout=1000000.
       
       !ag (19Jan2016)
       HYMAP2_routing_struc(n)%rivout_pre=0.0
       HYMAP2_routing_struc(n)%rivdph_pre=0.0
       HYMAP2_routing_struc(n)%fldout_pre=0.0
       HYMAP2_routing_struc(n)%flddph_pre=0.0
       HYMAP2_routing_struc(n)%fldelv1=0.0
       
       !ag (4Feb2016)
       if(HYMAP2_routing_struc(n)%resopflag==1)then
          HYMAP2_routing_struc(n)%resopaltloc=real(HYMAP2_routing_struc(n)%imis)
          HYMAP2_routing_struc(n)%resopoutmin=real(HYMAP2_routing_struc(n)%imis)
          HYMAP2_routing_struc(n)%tresopalt=dble(HYMAP2_routing_struc(n)%imis)
          HYMAP2_routing_struc(n)%resopalt=real(HYMAP2_routing_struc(n)%imis)
       endif

       !ag (03May2017)
       HYMAP2_routing_struc(n)%ewat=0.0
       HYMAP2_routing_struc(n)%edif=0.0
       
       !ag(27Apr2020)
!       if(HYMAP2_routing_struc(n)%flowtype ==4)then
         HYMAP2_routing_struc(n)%drsto=0.0
         HYMAP2_routing_struc(n)%drout=0.0
!       endif

       !ag(8Aug2020)
       if(HYMAP2_routing_struc(n)%insertflag==1)then
          HYMAP2_routing_struc(n)%insertloc=real(HYMAP2_routing_struc(n)%imis)
          HYMAP2_routing_struc(n)%tinsert=dble(HYMAP2_routing_struc(n)%imis)
          HYMAP2_routing_struc(n)%insertdis=real(HYMAP2_routing_struc(n)%imis)
       endif

       !ag(7Apr2021)
       if(HYMAP2_routing_struc(n)%tidesflag==1)then
         HYMAP2_routing_struc(n)%ttides=real(HYMAP2_routing_struc(n)%imis)
         HYMAP2_routing_struc(n)%tides=real(HYMAP2_routing_struc(n)%imis)
       endif
    enddo
   
 
    write(LIS_logunit,*) '[INFO] Get HYMAP2 cell vector sequence'
    do n=1, LIS_rc%nnest
       call HYMAP2_get_sindex(LIS_rc%gnc(n),&
            LIS_rc%gnr(n),&
            LIS_rc%glbnroutinggrid(n),&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%nextx,&
            HYMAP2_routing_struc(n)%nexty,maskg,&
            HYMAP2_routing_struc(n)%sindex,&
            HYMAP2_routing_struc(n)%outlet_glb,&
            HYMAP2_routing_struc(n)%next_glb)
#if 0 
       ic = 0 
       do c=1,LIS_rc%gnc(n)
          do r=1,LIS_rc%gnr(n)
             if(HYMAP2_routing_struc(n)%nextx(c,r)/=HYMAP2_routing_struc(n)%imis.and.maskg(c,r)>0) then 
                ic = ic + 1
                print*, ic, & !c,r,HYMAP2_routing_struc(n)%sindex(c,r)
                     -56.875 + (r-1)*0.25, -82.875 + (c-1)*0.25
             endif
          enddo
       enddo
       stop
#endif

       call HYMAP2_get_seq(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            LIS_rc%gnc(n),LIS_rc%gnr(n),&
            LIS_ews_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1), &
            HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%nextx,&
            HYMAP2_routing_struc(n)%nexty,maskg,&
            HYMAP2_routing_struc(n)%sindex,&
            HYMAP2_routing_struc(n)%outlet,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,HYMAP2_routing_struc(n)%next)

    enddo
    
#if (defined SPMD)
    do n=1,LIS_rc%nnest    
       call MPI_ALLGATHERV(HYMAP2_routing_struc(n)%seqx,&
            LIS_rc%nroutinggrid(n),MPI_INTEGER,&
            HYMAP2_routing_struc(n)%seqx_glb,&
            LIS_routing_gdeltas(n,:),&
            LIS_routing_goffsets(n,:),&
            MPI_INTEGER,&
            LIS_mpi_comm,status)

       call MPI_ALLGATHERV(HYMAP2_routing_struc(n)%seqy,&
            LIS_rc%nroutinggrid(n),MPI_INTEGER,&
            HYMAP2_routing_struc(n)%seqy_glb,&
            LIS_routing_gdeltas(n,:),&
            LIS_routing_goffsets(n,:),&
            MPI_INTEGER,&
            LIS_mpi_comm,status)
    enddo
#endif

    do n=1,LIS_rc%nnest
       allocate(tmp_real(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(tmp_real_nz(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            HYMAP2_routing_struc(n)%nz))

    ctitle = 'HYMAP_river_width'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,&
               HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,&
               tmp_real,HYMAP2_routing_struc(n)%rivwth(:,1))
       else
          do m=1,LIS_rc%nensem(n)
             call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP2_routing_struc(n)%nseqall,&
                  HYMAP2_routing_struc(n)%imis,&
                  HYMAP2_routing_struc(n)%seqx,&
                  HYMAP2_routing_struc(n)%seqy,&
                  tmp_real,HYMAP2_routing_struc(n)%rivwth(:,m))             
          enddo
       endif
!    enddo

    ctitle = 'HYMAP_river_length'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%rivlen)
!    enddo
    
    ctitle = 'HYMAP_river_height'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,&
               HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,&
               tmp_real,HYMAP2_routing_struc(n)%rivhgt(:,1))  
       else
          do m=1,LIS_rc%nensem(n)
             call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP2_routing_struc(n)%nseqall,&
                  HYMAP2_routing_struc(n)%imis,&
                  HYMAP2_routing_struc(n)%seqx,&
                  HYMAP2_routing_struc(n)%seqy,&
                  tmp_real,HYMAP2_routing_struc(n)%rivhgt(:,m))  
          enddo
       endif
!    enddo
    
    ctitle = 'HYMAP_river_roughness'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,&
               HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,&
               tmp_real,HYMAP2_routing_struc(n)%rivman(:,1))
       else
          do m=1,LIS_rc%nensem(n)
             call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP2_routing_struc(n)%nseqall,&
                  HYMAP2_routing_struc(n)%imis,&
                  HYMAP2_routing_struc(n)%seqx,&
                  HYMAP2_routing_struc(n)%seqy,&
                  tmp_real,HYMAP2_routing_struc(n)%rivman(:,m))
          enddo
       endif
!    enddo
    
    ctitle = 'HYMAP_floodplain_height'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real(ctitle,n, HYMAP2_routing_struc(n)%nz,&
            tmp_real_nz)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            HYMAP2_routing_struc(n)%nz,&
            HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real_nz,HYMAP2_routing_struc(n)%fldhgt)
!    enddo

    ctitle = 'HYMAP_floodplain_roughness'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%fldman)
!    enddo

    ctitle = 'HYMAP_grid_elevation'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       where(tmp_real<0)tmp_real=0
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%elevtn)
!    enddo

    ctitle = 'HYMAP_grid_distance'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%nxtdst)

!    enddo

    ctitle = 'HYMAP_grid_area'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%grarea)
!    enddo

    ctitle = 'HYMAP_runoff_time_delay'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%cntime)
       where(HYMAP2_routing_struc(n)%cntime==0)HYMAP2_routing_struc(n)%cntime=minval(HYMAP2_routing_struc(n)%cntime,HYMAP2_routing_struc(n)%cntime>0)
!    enddo


    ctitle = 'HYMAP_runoff_time_delay_multiplier'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%trnoff)
!    enddo

    ctitle = 'HYMAP_baseflow_time_delay'
!    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%tbsflw)
       write(LIS_logunit,*)'[INFO] Remove zeros from groundwater time delay'
       where(HYMAP2_routing_struc(n)%tbsflw<HYMAP2_routing_struc(n)%cntime/86400..and.&
            HYMAP2_routing_struc(n)%cntime>0)&
            HYMAP2_routing_struc(n)%tbsflw=HYMAP2_routing_struc(n)%cntime/86400.
!    enddo

    !ag (23Nov2016)
    ctitle = 'HYMAP_runoff_dwi_ratio'
!    do n=1, LIS_rc%nnest
       if(HYMAP2_routing_struc(n)%dwiflag==1)then
          call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,&
               HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,&
               tmp_real,HYMAP2_routing_struc(n)%rnfdwi_ratio)
       else
          HYMAP2_routing_struc(n)%rnfdwi_ratio=0.
       endif
!    enddo

    !ag (23Nov2016)
    ctitle = 'HYMAP_baseflow_dwi_ratio'
!    do n=1, LIS_rc%nnest
       if(HYMAP2_routing_struc(n)%dwiflag==1)then
          call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,&
               HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,&
               tmp_real,HYMAP2_routing_struc(n)%bsfdwi_ratio)
       else
          HYMAP2_routing_struc(n)%bsfdwi_ratio=0.
       endif
!    enddo

    !ag (13Apr2016)
    ctitle = 'HYMAP_river_flow_type'
!    do n=1, LIS_rc%nnest
       !ag(27Apr2020)
       !if(HYMAP2_routing_struc(n)%flowtype==0)then
       if(HYMAP2_routing_struc(n)%flowtype==0.or.HYMAP2_routing_struc(n)%flowtype==4)then
          call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,&
               HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,&
               tmp_real,HYMAP2_routing_struc(n)%flowmap) 
       else
          HYMAP2_routing_struc(n)%flowmap=HYMAP2_routing_struc(n)%flowtype
       endif
!    enddo

    !ag (7Dec2020)
    ctitle = 'HYMAP_urban_drainage_outlet'
!    do n=1, LIS_rc%nnest
       if(HYMAP2_routing_struc(n)%flowtype==4)then
          call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%imis,&
               HYMAP2_routing_struc(n)%seqx,&
               HYMAP2_routing_struc(n)%seqy,&
               tmp_real,HYMAP2_routing_struc(n)%droutlet)
       else
          HYMAP2_routing_struc(n)%droutlet=HYMAP2_routing_struc(n)%imis
       endif

#if (defined SPMD)
       call MPI_ALLGATHERV(HYMAP2_routing_struc(n)%droutlet,&
            LIS_rc%nroutinggrid(n),MPI_INTEGER,&
            HYMAP2_routing_struc(n)%droutlet_glb,&
            LIS_routing_gdeltas(n,:),&
            LIS_routing_goffsets(n,:),&
            MPI_INTEGER,&
            LIS_mpi_comm,status)
#endif

      deallocate(tmp_real)
      deallocate(tmp_real_nz)
    enddo


    !ag (20Sep2016) Correction for cases where parameter maps don't match
    do n=1, LIS_rc%nnest
       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          do i=1,HYMAP2_routing_struc(n)%nseqall
             if(HYMAP2_routing_struc(n)%rivwth(i,1)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%rivlen(i)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%rivhgt(i,1)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%rivman(i,1)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%fldhgt(i,1)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%elevtn(i)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%nxtdst(i)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%grarea(i)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%cntime(i)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%trnoff(i)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%tbsflw(i)==HYMAP2_routing_struc(n)%imis.or.&
                  HYMAP2_routing_struc(n)%flowmap(i)==HYMAP2_routing_struc(n)%imis)then
                
                HYMAP2_routing_struc(n)%rivwth(i,1)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%rivlen(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%rivhgt(i,1)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%rivman(i,1)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%fldhgt(i,:)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%elevtn(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%nxtdst(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%grarea(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%cntime(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%trnoff(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%tbsflw(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%flowmap(i)=HYMAP2_routing_struc(n)%imis
                HYMAP2_routing_struc(n)%outlet(i)=HYMAP2_routing_struc(n)%imis
             endif
          enddo
       else
          do i=1,HYMAP2_routing_struc(n)%nseqall
             do m=1,LIS_rc%nensem(n)
                if(HYMAP2_routing_struc(n)%rivwth(i,m)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%rivlen(i)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%rivhgt(i,m)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%rivman(i,m)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%fldhgt(i,1)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%elevtn(i)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%nxtdst(i)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%grarea(i)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%cntime(i)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%trnoff(i)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%tbsflw(i)==HYMAP2_routing_struc(n)%imis.or.&
                     HYMAP2_routing_struc(n)%flowmap(i)==HYMAP2_routing_struc(n)%imis)then
                   
                   HYMAP2_routing_struc(n)%rivwth(i,m)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%rivlen(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%rivhgt(i,m)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%rivman(i,m)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%fldhgt(i,:)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%elevtn(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%nxtdst(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%grarea(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%cntime(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%trnoff(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%tbsflw(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%flowmap(i)=HYMAP2_routing_struc(n)%imis
                   HYMAP2_routing_struc(n)%outlet(i)=HYMAP2_routing_struc(n)%imis
                endif
             enddo
          enddo
       endif
    enddo


    write(LIS_logunit,*) '[INFO] Processing data before running HYMAP'
    do n=1, LIS_rc%nnest
       write(LIS_logunit,*)'[INFO] Calculate maximum river storage'

       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          do i=1,HYMAP2_routing_struc(n)%nseqall
             HYMAP2_routing_struc(n)%rivstomax(i,1) = HYMAP2_routing_struc(n)%rivlen(i)* &
                  HYMAP2_routing_struc(n)%rivwth(i,1) * HYMAP2_routing_struc(n)%rivhgt(i,1)
          
!             write(LIS_logunit,*)'[INFO] Calculate river bed elevation'
             HYMAP2_routing_struc(n)%rivelv(i) = HYMAP2_routing_struc(n)%elevtn(i) -&
                  HYMAP2_routing_struc(n)%rivhgt(i,1)
!             write(LIS_logunit,*)'[INFO] Calculate river surface area'
             if(HYMAP2_routing_struc(n)%rivwth(i,1)>0) then 
                HYMAP2_routing_struc(n)%rivare(i,1) =&
                     min(HYMAP2_routing_struc(n)%grarea(i), &
                     HYMAP2_routing_struc(n)%rivlen(i) *&
                     HYMAP2_routing_struc(n)%rivwth(i,1))
!                write(LIS_logunit,*)'[INFO] Setting floodplain staging'
             endif
          enddo
          call HYMAP2_set_fldstg(HYMAP2_routing_struc(n)%nz,&
               HYMAP2_routing_struc(n)%nseqall,&
               HYMAP2_routing_struc(n)%fldhgt,&
               HYMAP2_routing_struc(n)%grarea,&
               HYMAP2_routing_struc(n)%rivlen,&
               HYMAP2_routing_struc(n)%rivwth(:,1),&
               HYMAP2_routing_struc(n)%rivstomax(:,1),&
               HYMAP2_routing_struc(n)%fldstomax(:,:,1),&
               HYMAP2_routing_struc(n)%fldgrd(:,:,1),&
               HYMAP2_routing_struc(n)%rivare(:,1))
       else
          do i=1,HYMAP2_routing_struc(n)%nseqall
             do m=1,LIS_rc%nensem(n)
                HYMAP2_routing_struc(n)%rivstomax(i,m) = HYMAP2_routing_struc(n)%rivlen(i)* &
                     HYMAP2_routing_struc(n)%rivwth(i,m) * HYMAP2_routing_struc(n)%rivhgt(i,m)
                
                HYMAP2_routing_struc(n)%rivelv(i) = HYMAP2_routing_struc(n)%elevtn(i) -&
                     HYMAP2_routing_struc(n)%rivhgt(i,m)
                if(HYMAP2_routing_struc(n)%rivwth(i,m)>0) then 
                   HYMAP2_routing_struc(n)%rivare(i,m) =&
                        min(HYMAP2_routing_struc(n)%grarea(i), &
                        HYMAP2_routing_struc(n)%rivlen(i) *&
                        HYMAP2_routing_struc(n)%rivwth(i,m))
                endif
             enddo
          enddo

          do m=1,LIS_rc%nensem(n)
             call HYMAP2_set_fldstg(HYMAP2_routing_struc(n)%nz,&
                  HYMAP2_routing_struc(n)%nseqall,&
                  HYMAP2_routing_struc(n)%fldhgt,&
                  HYMAP2_routing_struc(n)%grarea,&
                  HYMAP2_routing_struc(n)%rivlen,&
                  HYMAP2_routing_struc(n)%rivwth(:,m),&
                  HYMAP2_routing_struc(n)%rivstomax(:,m),&
                  HYMAP2_routing_struc(n)%fldstomax(:,:,m),&
                  HYMAP2_routing_struc(n)%fldgrd(:,:,m),&
                  HYMAP2_routing_struc(n)%rivare(:,m))
          enddo
       endif
    enddo

    !ag (4Feb2016) - read reservoir operation data
    do n=1, LIS_rc%nnest
!TBD: SVK - block below needs update for parallelism
       if(HYMAP2_routing_struc(n)%resopflag==1)then
          call HYMAP2_get_data_resop_alt(HYMAP2_routing_struc(n)%resopdir,&
               HYMAP2_routing_struc(n)%resopheader,&
               LIS_rc%gnc(n),LIS_rc%gnr(n),HYMAP2_routing_struc(n)%sindex,&
               HYMAP2_routing_struc(n)%nresop,HYMAP2_routing_struc(n)%ntresop,&
               HYMAP2_routing_struc(n)%resopaltloc,&
               HYMAP2_routing_struc(n)%resopoutmin,&
               HYMAP2_routing_struc(n)%tresopalt,&
               HYMAP2_routing_struc(n)%resopalt)
          where(HYMAP2_routing_struc(n)%resopoutmin==real(HYMAP2_routing_struc(n)%imis).or.HYMAP2_routing_struc(n)%resopoutmin<0.)&
               HYMAP2_routing_struc(n)%resopoutmin=0.
       endif
    enddo

    !ag(27Apr2020) - read urban flood data
    do n=1, LIS_rc%nnest
      if(HYMAP2_routing_struc(n)%flowtype==4)then
        write(LIS_logunit,*)'[INFO] Setting urban drainage parameters'
        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 urban drainage parameter file:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%drfile,rc=status)
        call LIS_verify(status,&
             "HYMAP2 urban drainage parameter file: not defined")

        call HYMAP2_get_urban_parameters(HYMAP2_routing_struc(n)%drfile,&
                 HYMAP2_routing_struc(n)%drwth,&
                 HYMAP2_routing_struc(n)%drhgt,&
                 HYMAP2_routing_struc(n)%drden,&
                 HYMAP2_routing_struc(n)%drvel,&
                 HYMAP2_routing_struc(n)%drblk,&
                 HYMAP2_routing_struc(n)%drrad,&
                 HYMAP2_routing_struc(n)%drlgh,&
                 HYMAP2_routing_struc(n)%drman,&
                 HYMAP2_routing_struc(n)%drslp)
                 
        call HYMAP2_gen_urban_drain_maps(HYMAP2_routing_struc(n)%nseqall,&
                  HYMAP2_routing_struc(n)%drrad,&
                  HYMAP2_routing_struc(n)%drlgh,&
                  HYMAP2_routing_struc(n)%drden,&
                  HYMAP2_routing_struc(n)%drwth,&
                  HYMAP2_routing_struc(n)%drblk,&
                  HYMAP2_routing_struc(n)%grarea,&
                  HYMAP2_routing_struc(n)%next,&
                  HYMAP2_routing_struc(n)%flowmap,&
                  HYMAP2_routing_struc(n)%drstomax,&
                  HYMAP2_routing_struc(n)%drtotwth,&
                  HYMAP2_routing_struc(n)%drnoutlet,&
                  HYMAP2_routing_struc(n)%drtotlgh)  
      endif
    enddo

    !ag(8Aug2020) - read discharge data for direct insertion
    do n=1, LIS_rc%nnest
!TBD: SVK - block below needs update for parallelism
       if(HYMAP2_routing_struc(n)%insertflag==1)then
          call HYMAP2_get_discharge_data(HYMAP2_routing_struc(n)%insertdir,&
               HYMAP2_routing_struc(n)%insertheader,&
               LIS_rc%gnc(n),LIS_rc%gnr(n),HYMAP2_routing_struc(n)%sindex,&
               HYMAP2_routing_struc(n)%ninsert,HYMAP2_routing_struc(n)%ntinsert,&
               HYMAP2_routing_struc(n)%insertloc,&
               HYMAP2_routing_struc(n)%tinsert,&
               HYMAP2_routing_struc(n)%insertdis)
       endif
    enddo

    !ag(30Mar2021) - read ocean tides data
    do n=1, LIS_rc%nnest
!TBD: SVK - block below needs update for parallelism
       if(HYMAP2_routing_struc(n)%tidesflag==1)then
          call HYMAP2_get_discharge_data(HYMAP2_routing_struc(n)%tidesdir,&
               HYMAP2_routing_struc(n)%tidesheader,&
               LIS_rc%gnc(n),LIS_rc%gnr(n),HYMAP2_routing_struc(n)%sindex,&
               HYMAP2_routing_struc(n)%ntides,HYMAP2_routing_struc(n)%nttides,&
               HYMAP2_routing_struc(n)%tidesloc,&
               HYMAP2_routing_struc(n)%ttides,&
               HYMAP2_routing_struc(n)%tides)
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model start mode:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP2_routing_struc(n)%startmode,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model start mode: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model restart interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model restart interval: not defined")
       call LIS_parseTimeString(time,HYMAP2_routing_struc(n)%rstInterval)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 routing model restart file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP2_routing_struc(n)%rstfile,rc=status)
       call LIS_verify(status,&
            "HYMAP2 routing model restart file: not defined")
    enddo

    if(LIS_rc%lsm.eq."none") then 
       call initrunoffdata(trim(LIS_rc%runoffdatasource)//char(0))
    endif

    do n=1, LIS_rc%nnest
       call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status)

!       global_grid_dg = ESMF_DistGridCreate(minIndex=(/1/),&
!            maxIndex=(/LIS_rc%glbntiles(n)/),&
!            regDecomp=(/1/),rc=status)
!       call LIS_verify(status,&
!            'ESMF_DistGridCreate failed in HYMAP2_routingInit')

!       global_grid = ESMF_GridCreate(name="Global Grid",&
!            coordTypeKind=ESMF_TYPEKIND_R4,&
!            distgrid = global_grid_dg, &
!            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
!       call LIS_verify(status,&
!            'ESMF_GridCreate failed in HYMAP2_routingInit')

       !create LSM interface objects to store runoff and baseflow
       sf_runoff_field =ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_vecTile(n), name="Surface Runoff",rc=status)
       call LIS_verify(status, 'ESMF_FieldCreate failed')

       baseflow_field =ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_vecTile(n), name="Subsurface Runoff",rc=status)
       call LIS_verify(status, 'ESMF_FieldCreate failed')

       call ESMF_FieldGet(sf_runoff_field,localDE=0,farrayPtr=sfrunoff,&
            rc=status)
       call LIS_verify(status)
       sfrunoff = 0.0 

       call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,&
            rc=status)
       call LIS_verify(status)
       baseflow = 0.0

       call ESMF_stateAdd(LIS_runoff_state(n),(/sf_runoff_field/),rc=status)
       call LIS_verify(status, 'ESMF_StateAdd failed for surface runoff')

       call ESMF_stateAdd(LIS_runoff_state(n),(/baseflow_field/),rc=status)
       call LIS_verify(status, 'ESMF_StateAdd failed for base flow')

       !hkb (4Mar2016)
       !create LSM interface objects to store evapotranspiration and 
       !potential evaporation only if source is readin
       if ( HYMAP2_routing_struc(n)%evapflag .ne. 0 ) then
          !if ( HYMAP2_routing_struc(n)%evapsrc .eq. "readin" ) then
          evapotranspiration_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=LIS_vecTile(n), name="Total Evapotranspiration",rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')

          !potential_evaporation_field =ESMF_FieldCreate(arrayspec=realarrspec,&
          !       grid=global_grid, name="Potential Evaporation",rc=status)
          !call LIS_verify(status, 'ESMF_FieldCreate failed')

          call ESMF_FieldGet(evapotranspiration_field,localDE=0,&
               farrayPtr=evapotranspiration,&
               rc=status)
          call LIS_verify(status)
          evapotranspiration = 0.0 

          !call ESMF_FieldGet(potential_evaporation_field,localDE=0, &
          !     farrayPtr=potevap,rc=status)
          !call LIS_verify(status)
          !potevap = 0.0 

          call ESMF_stateAdd(LIS_runoff_state(n),(/evapotranspiration_field/),rc=status)
          call LIS_verify(status, 'ESMF_StateAdd failed for Total Evapotranspiration')

          !call ESMF_stateAdd(LIS_runoff_state(n),(/potential_evaporation_field/),rc=status)
          !call LIS_verify(status, 'ESMF_StateAdd failed for potential evaporation')
          !endif
       endif

       HYMAP2_routing_struc(n)%mo = -1

       !ag (12Sep2019)
       call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling",&
            0, rc=status)
       call LIS_verify(status)

       if (HYMAP2_routing_struc(n)%enable2waycpl==1) then
           ! River Storage
           rivsto_field =ESMF_FieldCreate(arrayspec=realarrspec,&
                grid=LIS_vecTile(n), name="River Storage",rc=status)
           call LIS_verify(status, 'ESMF_FieldCreate failed')

           call ESMF_FieldGet(rivsto_field,localDE=0,farrayPtr=rivstotmp,&
                rc=status)
           call LIS_verify(status)
           rivstotmp = 0.0

           call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling",&
                HYMAP2_routing_struc(n)%enable2waycpl, rc=status)
           call LIS_verify(status)

           call ESMF_stateAdd(LIS_runoff_state(n),(/rivsto_field/),rc=status)
           call LIS_verify(status, 'ESMF_StateAdd failed for River Storage')

           ! Flood Storage
           fldsto_field =ESMF_FieldCreate(arrayspec=realarrspec,&
                grid=LIS_vecTile(n), name="Flood Storage",rc=status)
           call LIS_verify(status, 'ESMF_FieldCreate failed')

           call ESMF_FieldGet(fldsto_field,localDE=0,farrayPtr=fldstotmp,&
                rc=status)
           call LIS_verify(status)
           fldstotmp = 0.0

           call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling",&
                HYMAP2_routing_struc(n)%enable2waycpl, rc=status)
           call LIS_verify(status)

           call ESMF_stateAdd(LIS_runoff_state(n),(/fldsto_field/),rc=status)
           call LIS_verify(status, 'ESMF_StateAdd failed for Flood Storage')

           ! Flooded fraction
           fldfrc_field =ESMF_FieldCreate(arrayspec=realarrspec,&
                grid=LIS_vecTile(n), name="Flooded Fraction",rc=status)
           call LIS_verify(status, 'ESMF_FieldCreate failed')

           call ESMF_FieldGet(fldfrc_field,localDE=0,farrayPtr=fldfrctmp,&
                rc=status)
           call LIS_verify(status)
           fldfrctmp = 0.0

           call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling",&
                HYMAP2_routing_struc(n)%enable2waycpl, rc=status)
           call LIS_verify(status)

           call ESMF_stateAdd(LIS_runoff_state(n),(/fldfrc_field/),rc=status)
           call LIS_verify(status, 'ESMF_StateAdd failed for Flooded Fraction')
       endif

    enddo 

    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(LIS_runoff_state(n),"Routing model evaporation option",&
            HYMAP2_routing_struc(n)%evapflag, rc=status)
       call LIS_verify(status)
       
       call LIS_registerAlarm("HYMAP2 router model alarm",&
            HYMAP2_routing_struc(n)%dt,HYMAP2_routing_struc(n)%dt)
       
       call LIS_registerAlarm("HYMAP2 router output alarm",&
            HYMAP2_routing_struc(n)%dt,HYMAP2_routing_struc(n)%outInterval)
       
       call LIS_registerAlarm("HYMAP2 router restart alarm",&
            HYMAP2_routing_struc(n)%dt,HYMAP2_routing_struc(n)%rstInterval)          
    enddo

#if 0 
!test
    do n=1,LIS_rc%nnest
       do ic=1,HYMAP2_routing_struc(n)%nseqall 
          if(HYMAP2_routing_struc(n)%outlet(ic)==HYMAP2_routing_struc(n)%imis)cycle
          
          if(HYMAP2_routing_struc(n)%outlet(ic)==0)then
             ic_down=HYMAP2_routing_struc(n)%next(ic)
!             print*, LIS_localPet, ic, ic_down
          endif
       enddo
    enddo
#endif
    do n=1,LIS_rc%nnest
       allocate(deblklist(1,2,LIS_npes))
       do i=0,LIS_npes-1
          stid = LIS_routing_goffsets(n,i)*LIS_rc%nensem(n)+1
          enid = stid + LIS_routing_gdeltas(n,i)*LIS_rc%nensem(n)-1
          
          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)
       enddo

       patchDG = ESMF_DistGridCreate(minIndex=(/1/),&
            maxIndex=(/LIS_rc%glbnroutinggrid(n)*LIS_rc%nensem(n)/),&
            deBlockList=deblklist,rc=status)
       call LIS_verify(status)
       
       LIS_vecRoutingTile(n) = &
            ESMF_GridCreate(name="HYMAP2 Patch Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = patchDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LIS_verify(status,'ESMF_GridCreate failed in HYMAP2_routing_init')


       do i=0,LIS_npes-1
          stid = LIS_routing_goffsets(n,i)+1
          enid = stid + LIS_routing_gdeltas(n,i)-1
          
          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)
       enddo

       gridDG = ESMF_DistGridCreate(minIndex=(/1/),&
            maxIndex=(/LIS_rc%glbnroutinggrid(n)/),&
            deBlockList=deblklist,rc=status)
       call LIS_verify(status)
       
       LIS_vecRoutingGrid(n) = &
            ESMF_GridCreate(name="HYMAP2 Tile Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = gridDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LIS_verify(status,'ESMF_GridCreate failed in HYMAP2_routing_init')
       deallocate(deblklist)
    enddo

!---------------------------------------------------------------------
!  create the Routing state if data assimilation is being done, 
!  create the Routing perturbation state only if perturbation
!  option is turned on. 
!---------------------------------------------------------------------
    Routing_DAvalid = .false. 

    if(LIS_rc%ndas.gt.0.or.LIS_rc%nperts.gt.0) then 
       
       Routing_DAvalid = .true. 
       
       do i=1,LIS_rc%ndas
          Routing_DAvalid = Routing_DAvalid.and.LIS_rc%Routing_DAinst_valid(i)
       enddo

       allocate(LIS_Routing_State(LIS_rc%nnest, LIS_rc%nperts))
       allocate(LIS_Routing_Incr_State(LIS_rc%nnest, LIS_rc%nperts))
       
       do n=1,LIS_rc%nnest
          do k=1,LIS_rc%nperts
             write(LIS_logunit,*) &
                  '[INFO] Opening constraints for prognostic state variables: ',&
                  trim(LIS_rc%progattribFile(k))
             ftn = LIS_getNextUnitNumber()
             open(ftn, file = LIS_rc%progattribFile(k),status='old')
             read(ftn,*)
             read(ftn,*) LIS_rc%nstvars(k)
             read(ftn,*)
             
             allocate(vname(LIS_rc%nstvars(k)))
             allocate(stmin(LIS_rc%nstvars(k)))
             allocate(stmax(LIS_rc%nstvars(k)))
             
             call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
                  rc=status)
             call LIS_verify(status, &
                  "ESMF_ArraySpecSet failed in LIS_routing_init")
             
             write(unit=temp,fmt='(i2.2)') n
             read(unit=temp,fmt='(2a1)') nestid
             
             write(unit=temp,fmt='(i3.3)') k
             read(unit=temp,fmt='(3a1)') caseid

             LIS_Routing_State(n,k) = ESMF_StateCreate(name="Routing State"//&
                  nestid(1)//nestid(2)&
                  //'_'//caseid(1)//caseid(2)//caseid(3), rc=status)
             call LIS_verify(status, &
                  "ESMF_StateCreate failed in LIS_routing_init")

             LIS_Routing_Incr_State(n,k) = ESMF_StateCreate(name="Routing Incr State"//&
                  nestid(1)//nestid(2)// &
                  '_'//caseid(1)//caseid(2)//caseid(3), rc=status)
             call LIS_verify(status,&
                  "ESMF_StateCreate failed in LIS_routing_init")

             do i=1,LIS_rc%nstvars(k)
                read(ftn,fmt='(a40)') vname(i)
                read(ftn,*) stmin(i),stmax(i)
                write(LIS_logunit,*) '[INFO] ',vname(i),stmin(i),stmax(i)

                varField = ESMF_FieldCreate(&
                     grid=LIS_vecRoutingTile(n),&
                     arrayspec=arrspec1,name=trim(vname(i)), rc=status)
                call LIS_verify(status, &
                     "ESMF_FieldCreate failed in LIS_routing_init")

                varIncrField = ESMF_FieldCreate(&
                     grid=LIS_vecRoutingTile(n),&
                     arrayspec=arrspec1,name=trim(vname(i)), rc=status)
                call LIS_verify(status,&
                     "ESMF_FieldCreate failed in LIS_routing_init")

                call ESMF_AttributeSet(varField,"Max Value",stmax(i),rc=status)
                call LIS_verify(status,&
                     "ESMF_AttribteSet failed in LIS_routing_init")

                call ESMF_AttributeSet(varField,"Min Value",stmin(i),rc=status)
                call LIS_verify(status,&
                     "ESMF_AttributeSet failed in LIS_routing_init")


                call ESMF_AttributeSet(VarIncrField,"Max Value",stmax(i),rc=status)
                call LIS_verify(status,&
                     "ESMF_AttributeSet failed in LIS_routing_init")

                call ESMF_AttributeSet(VarIncrField,"Min Value",stmin(i),rc=status)
                call LIS_verify(status,&
                     "ESMF_AttributeSet failed in LIS_routing_init")

                call ESMF_StateAdd(LIS_Routing_State(n,k),(/varField/),rc=status)
                call LIS_verify(status,&
                     "ESMF_StateAdd failed in LIS_routing_init")

                call ESMF_StateAdd(LIS_Routing_Incr_State(n,k), &
                     (/VarIncrField/), rc=status)
                call LIS_verify(status,&
                     "ESMF_StateAdd failed in LIS_routing_init")
!----------------------------------------------------------------------------
! Initially set the fresh increments available status to false. 
!----------------------------------------------------------------------------
                call ESMF_AttributeSet(LIS_Routing_Incr_State(n,k), &
                     name="Fresh Increments Status", value=.false., &
                     rc=status)
                call LIS_verify(status,&
                     "ESMF_AttributeSet failed in LIS_routing_init")
             enddo
             deallocate(vname)
             deallocate(stmin)
             deallocate(stmax)
             call LIS_releaseUnitNumber(ftn)
          enddo
       enddo
    endif
    
    if(LIS_rc%nperts.gt.0) then 
       allocate(LIS_Routing_Pert_State(LIS_rc%nnest, LIS_rc%nperts))
       
       call ESMF_ArraySpecSet(arrspec2,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status,&
            "ESMF_ArraySpecSet failed in LIS_routing_init")

       do n=1,LIS_rc%nnest
          allocate(ssdev(LIS_rc%ngrid(n)))
          do k=1,LIS_rc%nperts
             if(LIS_rc%perturb_state(k).ne."none") then 
                allocate(routing_pert%vname(LIS_rc%nstvars(k)))
                allocate(routing_pert%perttype(LIS_rc%nstvars(k)))
                allocate(routing_pert%ssdev(LIS_rc%nstvars(k)))
                allocate(routing_pert%stdmax(LIS_rc%nstvars(k)))
                allocate(routing_pert%zeromean(LIS_rc%nstvars(k)))
                allocate(routing_pert%tcorr(LIS_rc%nstvars(k)))
                allocate(routing_pert%xcorr(LIS_rc%nstvars(k)))
                allocate(routing_pert%ycorr(LIS_rc%nstvars(k)))
                allocate(routing_pert%ccorr(LIS_rc%nstvars(k),LIS_rc%nstvars(k)))

                write(unit=temp,fmt='(i2.2)') n
                read(unit=temp,fmt='(2a1)') nestid

                LIS_Routing_Pert_State(n,k) = ESMF_StateCreate(&
                     name="Routing_Pert_State"//&
                     nestid(1)//nestid(2),&
                     rc=status)
                call LIS_verify(status,&
                     "ESMF_StateCreate: Routing_Pert_State failed in LIS_routing_init")

                call LIS_readPertAttributes(LIS_rc%nstvars(k),&
                     LIS_rc%progpertAttribfile(k),&
                     routing_pert)

                do i=1,LIS_rc%nstvars(k)
                   pertField = ESMF_FieldCreate(&
                        grid=LIS_vecRoutingTile(n),&
                        arrayspec=arrspec2,name=trim(routing_pert%vname(i)),&
                        rc=status)

                   call ESMF_StateAdd(LIS_Routing_Pert_State(n,k),(/pertField/),&
                        rc=status)
                   call LIS_verify(status,&
                        "ESMF_StateAdd failed in LIS_routing_init")
                enddo

                allocate(pertobjs(LIS_rc%nstvars(k)))
                allocate(order(LIS_rc%nstvars(k)))
                allocate(ccorr(LIS_rc%nstvars(k),LIS_rc%nstvars(k)))
                order = -1

                call ESMF_StateGet(LIS_Routing_Pert_State(n,k),&
                     itemNameList=pertobjs,rc=status)
                call LIS_verify(status,&
                     "ESMF_StateGet failed in LIS_routing_init")

                do i=1,LIS_rc%nstvars(k)
                   do j=1,LIS_rc%nstvars(k)
                      if(routing_pert%vname(j).eq.pertobjs(i)) then 
                         order(i) = j
                         exit;
                      endif
                   enddo
                enddo

                do i=1,LIS_rc%nstvars(k)
                   do j=1,LIS_rc%nstvars(k)
                      ccorr(i,j) = routing_pert%ccorr(order(i),order(j))
                   enddo
                enddo

                do i=1,LIS_rc%nstvars(k)
                   call ESMF_StateGet(LIS_Routing_Pert_State(n,k),&
                        pertobjs(i),pertField,rc=status)
                   call LIS_verify(status,&
                        "ESMF_StateGet failed in LIS_routing_init")

                   call ESMF_AttributeSet(pertField,"Perturbation Type",&
                        routing_pert%perttype(order(i)),&
                        rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet: Perturbation Type failed in LIS_routing_init")

                   if(LIS_rc%ngrid(n).gt.0) then 
                      ssdev = routing_pert%ssdev(order(i))

                      call ESMF_AttributeSet(pertField,"Standard Deviation",&
                           ssdev,itemCount=LIS_rc%ngrid(n),&
                           rc=status)
                      call LIS_verify(status,&
                           "ESMF_AttributeSet: Standard Deviation failed in LIS_routing_init")
                   endif
                   call ESMF_AttributeSet(pertField,"Std Normal Max",&
                        routing_pert%stdmax(order(i)),&
                        rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet: Std Normal Max failed in LIS_routing_init")

                   call ESMF_AttributeSet(pertField,"Ensure Zero Mean",&
                        routing_pert%zeromean(order(i)),&
                        rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet: Ensure Zero Mean failed in LIS_routing_init")

                   call ESMF_AttributeSet(pertField,&
                        "Temporal Correlation Scale",&
                        routing_pert%tcorr(order(i)), rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet: Temporal Correlation Scale failed in LIS_routing_init")

                   call ESMF_AttributeSet(pertField,"X Correlation Scale",&
                        routing_pert%xcorr(order(i)), rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet: X Correlation Scale failed in LIS_routing_init")

                   call ESMF_AttributeSet(pertField,"Y Correlation Scale",&
                        routing_pert%ycorr(order(i)), rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet: Y Correlation Scale failed in LIS_routing_init")

                   call ESMF_AttributeSet(pertField,&
                        "Cross Correlation Strength",&
                        ccorr(i,:), itemCount=LIS_rc%nstvars(k),&
                        rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet: Cross Correlation Strength failed in LIS_routing_init")

                enddo
                deallocate(pertobjs)
                deallocate(order)
                deallocate(ccorr)
                deallocate(routing_pert%vname)
                deallocate(routing_pert%perttype)
                deallocate(routing_pert%ssdev)
                deallocate(routing_pert%stdmax)
                deallocate(routing_pert%zeromean)
                deallocate(routing_pert%tcorr)
                deallocate(routing_pert%xcorr)
                deallocate(routing_pert%ycorr)
                deallocate(routing_pert%ccorr)
             endif
          enddo
          deallocate(ssdev)          
       enddo
    endif
    if(Routing_DAvalid) then 
       max_index = -1
       do i=1,LIS_rc%nperts
          if(max_index.eq.-1.and.LIS_rc%perturb_state(i).ne."none") then 
             max_index = 1
             alglist(max_index) = LIS_rc%perturb_state(i)
          else
             name_found = .false. 
             do k=1,max_index
                if(LIS_rc%perturb_state(i).ne."none".and.&
                     LIS_rc%perturb_state(i).eq.alglist(k)) then
                   name_found = .true. 
                endif
             enddo
             if(.not.name_found.and.max_index.ne.-1) then 
                max_index = max_index + 1
                alglist(max_index) = LIS_rc%perturb_state(i)
             endif
          endif
       enddo

       if(max_index.gt.0) then 
          do i=1,max_index
             !Call this only once for all instances of the algorithm
             call perturbinit(trim(alglist(i))//char(0), 4)
          enddo
       endif

       do i=1,LIS_rc%nperts   
          if(LIS_rc%perturb_state(i).ne."none") then 
             call perturbsetup(trim(LIS_rc%perturb_state(i))//char(0), 4, i, &
                  LIS_Routing_State(:,i), LIS_Routing_Pert_State(:,i))
          endif
       enddo
    end if
  end subroutine HYMAP2_routingInit
  !=============================================
  !=============================================
  subroutine HYMAP2_read_param_real(ctitle,n,z,array)
    
    !USES: 
    use LIS_coreMod
    use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    
    implicit none	  
    character(*), intent(in)    :: ctitle
    integer,      intent(in)    :: z
    integer,      intent(in)    :: n 
    real,         intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n),z)
    integer                     :: ftn 
    logical                     :: file_exists
    integer                     :: status
    integer                     :: l
    integer                     :: varid
    character*100               :: cfile
!    real                        :: array1(LIS_rc%gnc(n),LIS_rc%gnr(n),z)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)
    
    if(file_exists) then 

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_real in HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1),z/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP2_routingMod')
       
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routingMod')
       
!       do l=1,z
!          array(:,:,l) = nint(array1(&
!               LIS_ews_halo_ind(n,LIS_localPet+1):&         
!               LIS_ewe_halo_ind(n,LIS_localPet+1), &
!               LIS_nss_halo_ind(n,LIS_localPet+1): &
!               LIS_nse_halo_ind(n,LIS_localPet+1),l))
!       enddo
    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_real in HYMAP2_routingMod'
       call LIS_endrun()
    endif
    
#endif
  end subroutine HYMAP2_read_param_real

  !=============================================
  subroutine HYMAP2_read_param_real_2d(ctitle,n,array)
    
    !USES: 
    use LIS_coreMod
    use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    
    implicit none	  
    character(*), intent(in)    :: ctitle
    integer,      intent(in)    :: n 
    real,         intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
    integer                     :: ftn 
    logical                     :: file_exists
    integer                     :: status
    integer                     :: varid
    character*100               :: cfile

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)
    
    if(file_exists) then 

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_real in HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1)/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1)/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP2_routingMod')
       
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routingMod')
       
    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_real in HYMAP2_routingMod'
       call LIS_endrun()
    endif
    
#endif
  end subroutine HYMAP2_read_param_real_2d

  subroutine HYMAP2_read_param_int(ctitle,z,n,array)
    
    !USES: 
    use LIS_coreMod
    use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    
    implicit none	  
    character(*), intent(in)    :: ctitle
    integer,      intent(in)    :: z
    integer,      intent(in)    :: n 
    integer,      intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n),z)
    integer                     :: ftn 
    logical                     :: file_exists
    integer                     :: status
    integer                     :: c,r,l
    character*100               :: cfile
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)
    
    if(file_exists) then 

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_int in HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1),z/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP2_routingMod')
       
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_int in HYMAP2_routingMod'
       call LIS_endrun()
    endif
    
#endif
  end subroutine HYMAP2_read_param_int

  subroutine HYMAP2_read_param_int_2d(ctitle,n,array)
    
    !USES: 
    use LIS_coreMod
    use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    
    implicit none	  
    character(*), intent(in)    :: ctitle
    integer,      intent(in)    :: n 
    integer,      intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
    integer                     :: ftn 
    logical                     :: file_exists
    integer                     :: status
    integer                     :: c,r
    character*100               :: cfile
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)
    
    if(file_exists) then 

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_int in HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP2_routingMod')
       

       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1)/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1)/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP2_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_int in HYMAP2_routingMod'
       call LIS_endrun()
    endif
    
#endif
  end subroutine HYMAP2_read_param_int_2d

  subroutine HYMAP2_read_param_int_2d_global(ctitle,n,array)
    
    !USES: 
    use LIS_coreMod
    use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    
    implicit none	  
    character(*), intent(in)    :: ctitle
    integer,      intent(in)    :: n 
    integer,      intent(inout) :: array(LIS_rc%gnc(n),LIS_rc%gnr(n))
    integer                     :: ftn 
    logical                     :: file_exists
    integer                     :: status
    integer                     :: c,r
    character*100               :: cfile
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)
    
    if(file_exists) then 

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_int in HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array), &
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP2_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_int in HYMAP2_routingMod'
       call LIS_endrun()
    endif
    
#endif
  end subroutine HYMAP2_read_param_int_2d_global
  
  !=============================================
  !=============================================
  subroutine HYMAP2_get_data_resop_alt(resopdir,resopheader,nx,ny,sindex,nresop,ntresop,resoploc,resopoutmin,tresop,resop)
  
    implicit none
    character(*), intent(in)  :: resopdir,resopheader
    integer,      intent(in)  :: nx,ny
    integer,      intent(in)  :: sindex(nx,ny)
    integer,      intent(in)  :: nresop,ntresop
    integer,      intent(out) :: resoploc(nresop)
    real*8,       intent(out) :: tresop(nresop,ntresop)
    real,         intent(out) :: resop(nresop,ntresop),resopoutmin(nresop)
    integer                   :: res
    integer                   :: xresop(nresop),yresop(nresop)
    character(50)            :: resopname(nresop)
    character(len=LIS_CONST_PATH_LEN) :: yfile
    
    call HYMAP2_read_header_resop(trim(resopheader),nresop,resopname,xresop,yresop,resopoutmin)
    do res=1,nresop
      resoploc(res)=sindex(xresop(res),yresop(res))
    enddo
    do res=1,nresop
      yfile=trim(resopdir)//trim(resopname(res))//'.txt'
      call HYMAP2_read_time_series(ntresop,trim(yfile),tresop(res,:),resop(res,:))
    enddo
  end subroutine HYMAP2_get_data_resop_alt
  !=============================================
  !=============================================
  subroutine HYMAP2_get_discharge_data(insertdir,insertheader,nx,ny,sindex,ninsert,ntinsert,insertloc,tinsert,insertdis)

    implicit none
    character(*), intent(in)  :: insertdir,insertheader
    integer,      intent(in)  :: nx,ny
    integer,      intent(in)  :: sindex(nx,ny)
    integer,      intent(in)  :: ninsert,ntinsert
    integer,      intent(out) :: insertloc(ninsert)
    real*8,       intent(out) :: tinsert(ninsert,ntinsert)
    real,         intent(out) :: insertdis(ninsert,ntinsert)
    integer                   :: ins
    integer                   :: xinsert(ninsert),yinsert(ninsert)
    character(50)            :: insertname(ninsert)
    character(len=LIS_CONST_PATH_LEN) :: yfile

    call HYMAP2_read_header(trim(insertheader),ninsert,insertname,xinsert,yinsert)
    do ins=1,ninsert
      insertloc(ins)=sindex(xinsert(ins),yinsert(ins))
    enddo
    do ins=1,ninsert
      yfile=trim(insertdir)//trim(insertname(ins))//'.txt'
      call HYMAP2_read_time_series(ntinsert,trim(yfile),tinsert(ins,:),insertdis(ins,:))
    enddo
  end subroutine HYMAP2_get_discharge_data
  !=============================================
  !=============================================
!  subroutine HYMAP2_get_tides_data(tidesdir,tidesheader,nx,ny,sindex,ntides,nttides,tidesloc,ttides,tides)
  
!    implicit none
!    character(*), intent(in)  :: tidesdir,tidesheader
!    integer,      intent(in)  :: nx,ny
!    integer,      intent(in)  :: sindex(nx,ny)
!    integer,      intent(in)  :: ntides,nttides
!    integer,      intent(out) :: tidesloc(ntides)
!    real*8,       intent(out) :: ttides(ntides,nttides)
!    real,         intent(out) :: tides(ntides,nttides)
!    integer                   :: ins
!    integer                   :: xtides(ntides),ytides(ntides)
!    character(50)            :: itidesname(ntides)
!    character(500)            :: yfile
    
!    call HYMAP2_read_header(trim(insertheader),ninsert,insertname,xinsert,yinsert)
!    do ins=1,ninsert
!      insertloc(ins)=sindex(xinsert(ins),yinsert(ins))
!    enddo
!    do ins=1,ninsert
!      yfile=trim(insertdir)//trim(insertname(ins))//'.txt'
!      call HYMAP2_read_time_series(ntinsert,trim(yfile),tinsert(ins,:),insertdis(ins,:))
!    enddo
!  end subroutine HYMAP2_get_tides_data
  !=============================================
  !=============================================  
  subroutine HYMAP2_read_header_resop(yheader,inst,yqname,ix,iy,outmin)

    use LIS_logMod
    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(in)    :: inst
    character(*), intent(inout) :: yqname(inst)
    integer,      intent(inout) :: ix(inst),iy(inst)
    real,         intent(inout) :: outmin(inst)   
    logical                     :: file_exists
    integer                     :: ist

    inquire(file=yheader,exist=file_exists)
    
    if(file_exists) then 
      !get name of station files
      write(LIS_logunit,*)'[read_header] get stations info: name and coordinates'
      write(LIS_logunit,*)'[read_header] ',yheader,inst
      open(2,file=trim(yheader), status='old')
      !print*,inst
      do ist=1,inst
        read(2,*,end=10)yqname(ist),ix(ist),iy(ist),outmin(ist)
        write(LIS_logunit,'(a,i5,a,2i5,f10.2)')'[read_header] ',ist,trim(yqname(ist)),ix(ist),iy(ist),outmin(ist)
      enddo
      close(2)
    else
      write(LIS_logunit,*) 'header file '//trim(yheader)
      write(LIS_logunit,*) 'failed in read_header in HYMAP2_routingMod'
      call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) 'header file '//trim(yheader)
    write(LIS_logunit,*) 'failed in read_header in HYMAP2_routingMod'
    call LIS_endrun()

  end subroutine HYMAP2_read_header_resop
  !=============================================
  !=============================================  
  subroutine HYMAP2_read_header(yheader,inst,yqname,ix,iy)

    use LIS_logMod
    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(in)    :: inst
    character(*), intent(inout) :: yqname(inst)
    integer,      intent(inout) :: ix(inst),iy(inst)
    logical                     :: file_exists
    integer                     :: ist

    inquire(file=yheader,exist=file_exists)
    
    if(file_exists) then 
      !get name of station files
      write(LIS_logunit,*)'[read_header] get stations info: name and coordinates'
      write(LIS_logunit,*)'[read_header] ',yheader,inst
      open(2,file=trim(yheader), status='old')
      !print*,inst
      do ist=1,inst
        read(2,*,end=10)yqname(ist),ix(ist),iy(ist)
        write(LIS_logunit,*)'[read_header] ',ist,trim(yqname(ist)),ix(ist),iy(ist)
      enddo
      close(2)
    else
      write(LIS_logunit,*) 'header file '//trim(yheader)
      write(LIS_logunit,*) 'failed in read_header in HYMAP2_routingMod'
      call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) 'header file '//trim(yheader)
    write(LIS_logunit,*) 'failed in read_header in HYMAP2_routingMod'
    call LIS_endrun()

  end subroutine HYMAP2_read_header
  !=============================================
  !=============================================  
  subroutine HYMAP2_read_time_series(itmax,yfile,ztalt,zhalt)
  
    use LIS_logMod
    implicit none

    integer,      intent(in)    :: itmax
    character(*), intent(in)    :: yfile
    real*8,       intent(inout) :: ztalt(itmax)
    real,         intent(inout) :: zhalt(itmax)
    logical                     :: file_exists
    integer :: it, ftn
  
    inquire(file=yfile,exist=file_exists)
    
    if(file_exists)then 
      ftn = LIS_getNextUnitNumber()
      open(ftn,file=trim(yfile),status='old')
      do it=1,itmax
        read(ftn,*,end=10)ztalt(it),zhalt(it)
        write(LIS_logunit,*)ztalt(it),zhalt(it)
      enddo
10    continue
      call LIS_releaseUnitNumber(ftn)
    else
      write(LIS_logunit,*) '[ERR] time series file '//trim(yfile)
      write(LIS_logunit,*) '[ERR] failed in read_time_series in HYMAP2_routing_init'
      call LIS_endrun()
    endif  

  end subroutine HYMAP2_read_time_series
  !============================================= 
!ag(27Apr2020)
! ================================================
!BOP
!
! !ROUTINE: HYMAP2_gather_tiles
! \label{HYMAP2_gather_tiles}
! 
! !INTERFACE:
subroutine HYMAP2_gather_tiles_int(n,var,var_glb)
! !USES:
  use LIS_coreMod
  use LIS_routingMod
  use LIS_mpiMod
!
! !DESCRIPTION: 
!  This subroutine gathers an individual variable
!  across different processors into a global array
!EOP

  implicit none

  integer        :: n 
  integer           :: var(LIS_rc%nroutinggrid(n))
  integer           :: var_glb(LIS_rc%glbnroutinggrid(n))

  integer           :: tmpvar(LIS_rc%glbnroutinggrid(n))
  integer        :: i,l,ix,iy,ix1,iy1
  integer        :: status

#if (defined SPMD)
  call MPI_ALLGATHERV(var,&
       LIS_rc%nroutinggrid(n),&
       MPI_REAL,tmpvar,&
       LIS_routing_gdeltas(n,:),&
       LIS_routing_goffsets(n,:),&
       MPI_REAL,LIS_mpi_comm,status)
#endif
  !rearrange them to be in correct order.
  do l=1,LIS_npes
     do i=1,LIS_routing_gdeltas(n,l-1)
        ix = HYMAP2_routing_struc(n)%seqx_glb(i+&
             LIS_routing_goffsets(n,l-1))
        iy = HYMAP2_routing_struc(n)%seqy_glb(i+&
             LIS_routing_goffsets(n,l-1))
        ix1 = ix + LIS_ews_halo_ind(n,l) - 1
        iy1 = iy + LIS_nss_halo_ind(n,l)-1
        var_glb(HYMAP2_routing_struc(n)%sindex(ix1,iy1)) = &
             tmpvar(i+LIS_routing_goffsets(n,l-1))
     enddo
  enddo

end subroutine HYMAP2_gather_tiles_int
! ================================================ 
end module HYMAP2_routingMod
