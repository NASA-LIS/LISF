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
module HYMAP3_routingMod
!BOP
!
! !MODULE: HYMAP3_routingMod
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
! 08 Nov 2011: Augusto Getirana, Initial implementation in LIS based on
!   the HYMAP offline routing code.
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia
!   formulation, adaptive time step and reservoir operation.
! 13 Apr 2016: Augusto Getirana, Inclusion of option for hybrid runs with
!   a river flow map.
! 07 Sep 2019: Augusto Getirana,  Added support for 2-way coupling
! 27 Apr 2020: Augusto Getirana,  Added support for urban drainage
!
! !USES:
  use ESMF
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_topoMod

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: HYMAP3_routingInit
  public :: HYMAP3_logunit ! file unit number used for diagnostic logging
  public :: HYMAP3_gather_tiles_int
  public :: HYMAP3_gather_tiles
  public :: HYMAP3_map_g2l
  public :: HYMAP3_map_gxy2l_index
  public :: HYMAP3_map_l2g_index
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

  public :: HYMAP3_routing_struc

  type, public :: HYMAP3_routing_dec

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
     real,    allocatable :: elevtn(:)      !grid elevation [m]
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
     character(LIS_CONST_PATH_LEN) :: rstfile
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
     real,   allocatable  :: edif(:,:) !differential evapotranspiration
                                       ! (evaporation from open waters -
                                       ! total evapotranspiration) used
                                       ! as input in HyMAP [km m-2 s-1]

     integer              :: dwiflag      !deep water infiltration flag
     !ag (11Sep2015)
     integer              :: steptype     !time step type flag

     character(LIS_CONST_PATH_LEN) :: LISdir     !if LIS output is being read

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
     integer, allocatable :: resoploc(:)
     integer, allocatable :: resoploc_dwn(:)
     integer, allocatable :: resoploc_glb(:)
     integer, allocatable :: resoploc_dwn_glb(:)
     real,    allocatable :: resopoutmin(:)
     real*8,  allocatable :: tresopalt(:,:)
     real,    allocatable :: resopalt(:,:)
     integer              :: nresop   !number of reservoirs
     integer              :: ntresop  !time series length (number of time
                                      ! steps in the input files)
     integer              :: resopflag
     character(LIS_CONST_PATH_LEN) :: resopdir
     character(LIS_CONST_PATH_LEN) :: resopheader
     !ag (29Jun2016)
     integer              :: floodflag
     character(LIS_CONST_PATH_LEN) :: HYMAP_dfile
     !ag(17Apr2024)
     integer, allocatable :: resoptype(:)
     !ag(11Oct2024)
     real,   allocatable  :: elevtn_resop(:)
     real,   allocatable  :: fldhgt_resop(:,:)
     real,   allocatable  :: fldstomax_resop(:,:)
     real,   allocatable  :: grarea_resop(:)
     real,   allocatable  :: rivstomax_resop(:)
     real,   allocatable  :: rivelv_resop(:)
     real,   allocatable  :: rivlen_resop(:)
     real,   allocatable  :: rivwth_resop(:)
     ! === 2-way coupling variables/parameters ===
     real,   allocatable  :: rivstotmp(:,:)     !River Storage [m3]
     real,   allocatable  :: fldstotmp(:,:)     !Flood Storage [m3]
     real,   allocatable  :: fldfrctmp(:,:)     !Flooded Fraction [m3]
     integer              :: enable2waycpl
     real                 :: fldfrc2waycpl

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
     character(LIS_CONST_PATH_LEN) :: drfile          !urban drainage parametere file name
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
     character(LIS_CONST_PATH_LEN) :: insertdir
     character(LIS_CONST_PATH_LEN) :: insertheader

     !ag(30Mar2021)
     ! === sea level variables/parameters ===
     integer, allocatable  :: outletid(:)      !outlet identification
     real*8,  allocatable  :: tsealevel(:,:)   !date/time for sea level time series

     real,    allocatable  :: sealevel(:,:)    !sea level time series
     integer               :: nsealevel        !number of sealevel time series
     integer               :: ntsealevel       !time series length (number of time steps in the input files)
     integer               :: noutlet          !number of outlets
     integer               :: sealevelflag     !flag accounting for varying sea level
     character(LIS_CONST_PATH_LEN) :: sealeveldir      !directory containing files with sea level time series
     character(LIS_CONST_PATH_LEN) :: sealevelheader   !file containing list of sea level time series
     character(LIS_CONST_PATH_LEN) :: outletlist       !file containing list of outlet under varying sea level effect
     !ag
     ! === water management variables/parameters ===
     integer, allocatable  :: managact(:)
     integer, allocatable  :: managloc(:,:)
     integer, allocatable  :: managtype(:)
     real,    allocatable  :: managqmax(:)
     real,    allocatable  :: managcoef(:,:)
     integer               :: nmanagcoef        !number of management coefficients
     integer               :: nmanag            !number of locations with water management
     integer               :: managflag
     character(LIS_CONST_PATH_LEN) :: managheader
     !ag(30Mar2022)
     ! === bifurcation variables/parameters ===
     integer               :: bifflag          !bifurcation flag
     character(LIS_CONST_PATH_LEN) :: biffile          !bifurcation pathway input file
     integer               :: nbif             !number of bifurcations within the domain
     integer               :: nbifelv          !number of pathways in one bifurcation (defined by elevations)
     real,    allocatable  :: bifout(:)        !bifurcation streamflow [m3/s]
     real,    allocatable  :: bifout_pre(:)    !previousbifurcation streamflow [m3/s]
     real,    allocatable  :: bifwth(:,:)      !pathway width for each elevation [m]
     real,    allocatable  :: bifelv(:)        !pathway reference elevation [m]
     real,    allocatable  :: bifdelv(:)       !pathway delta elevations, having bifelv as reference [m]
     real,    allocatable  :: bifman(:)        !pathway roughness coefficients [-]
     real,    allocatable  :: biflen(:)        !pathway length [m]
     real,    allocatable  :: bifsto(:,:)      !bifurcation water storage for each elevation [m3]
     integer, allocatable  :: bifloc(:,:)      !location of upstream and downstream grids composing bifurcation

     !ag(23Feb2023)
     ! === levee variables/parameters ===
     integer               :: levflag          !levee flag
     real,    allocatable  :: levhgt(:)        !levee elevation above surface elevation [m]
     real,    allocatable  :: levstomax(:,:)   !levee elevation above surface elevation [m3]
     real,    allocatable  :: fldonlystomax(:,:,:) !maximum floodplain-only storage - excludes river storage [m3]
     real,    allocatable  :: fldstoatlev(:,:) !isolated floodplain storage at levee elevation - excludes river storage [m3]

     !ag(8Nov2024)
     ! === slope constraints ===
     real,    allocatable  :: rivslp(:)        !riverbed slope [m/m]
     real,    allocatable  :: maxsfcslp(:)     !upper threshold (max. allowed) surface slope [m/m]

     !ag(4Apr2025)
     ! === Yassin's reservoir operation scheme ===
     character(LIS_CONST_PATH_LEN) :: resopncfile
     integer,  allocatable :: ncloc_resop(:)
     real*8,   allocatable :: maxsto_resop(:)
     real*8,   allocatable :: inidis_resop(:)
     real*8,   allocatable :: inisto_resop(:)
     real*8,   allocatable :: dwndis_resop(:)
     real*8,   allocatable :: deadis_resop(:)

     real*8,   allocatable :: minsto_mo_resop(:,:)
     real*8,   allocatable :: nupsto_mo_resop(:,:)
     real*8,   allocatable :: uppsto_mo_resop(:,:)
     real*8,   allocatable :: mindis_mo_resop(:,:)
     real*8,   allocatable :: nupdis_mo_resop(:,:)
     real*8,   allocatable :: uppdis_mo_resop(:,:)

     real*8,   allocatable :: reg1_resop(:)
     real*8,   allocatable :: reg2_resop(:)
     real*8,   allocatable :: reg3_resop(:)

     !ag(14Apr2025)
     ! === vector-based HyMAP implementation ===
     character(LIS_CONST_PATH_LEN) :: vecfile  !vector input file
     integer               :: vecflag          !vector input flag

  end type HYMAP3_routing_dec

  type(HYMAP3_routing_dec), allocatable :: HYMAP3_routing_struc(:)
  integer :: HYMAP3_logunit ! file unit number used for diagnostic logging

contains

!BOP
!
! !ROUTINE: HYMAP3_routingInit
! \label{HYMAP3_routingInit}
!
  subroutine HYMAP3_routingInit
    !USES:
    !ag(6Apr2022)
    use HYMAP3_bifMod
    use HYMAP3_initMod
    !ag(1May2021)
    use HYMAP3_managMod
    use HYMAP3_modelMod
    !ag(27Apr2020)
    use HYMAP3_urbanMod
    use LIS_coreMod
    use LIS_logMod
    use LIS_mpiMod
    use LIS_perturbMod
    use LIS_routingMod
    use LIS_timeMgrMod

    implicit none

    integer              :: n
    integer              :: ftn
    integer              :: status
    integer              :: i,j,k
    type(ESMF_ArraySpec) :: realarrspec
    type(ESMF_Field)     :: sf_runoff_field
    type(ESMF_Field)     :: baseflow_field
    type(ESMF_Field)     :: evapotranspiration_field
    real, pointer        :: sfrunoff(:)
    real, pointer        :: baseflow(:)
    real, pointer        :: evapotranspiration(:)
    character*100        :: ctitle
    character*10         :: time
    character*20         :: flowtype
    !ag (11Sep2015)
    character*20         :: steptype
    !ag (31Jan2016)

    !ag (19Feb2016)
    real,    allocatable :: tmp_real(:,:),tmp_real_nz(:,:,:)
    integer, allocatable :: mask(:,:),maskg(:,:)
    real,    allocatable :: elevtn(:,:),uparea(:,:),basin(:,:)

    !ag (11Mar2016)
    integer, external  :: LIS_create_subdirs

    !ag (03Apr2017)
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

    !ag(7Sep2022)
    ! === bifurcation local variables/parameters ===
    real,    allocatable :: elevtn_glb(:)      !river bed elevation [m]
    real,    allocatable :: grarea_glb(:)      !area of the grid [m^2]
    real,    allocatable :: fldhgt_glb(:,:)    !floodplain height [m]
    real,    allocatable :: rivstomax_glb(:)   !maximum river storage [m3]
    real,    allocatable :: fldstomax_glb(:,:) !maximum floodplain storage [m3]
    real,    allocatable :: rivlen_glb(:)      !river length [m]
    real,    allocatable :: rivwth_glb(:)      !river width [m]
    real,    allocatable :: rivelv_glb(:)      !elevation of river bed [m]
    real*8               :: bifelv1,bifsto1
    integer              :: ibif,ielv,icg,iz

    !ag(11Oct2024)
    integer              :: iresop
    real                 :: cadp

    external :: initrunoffdata
    external :: perturbinit
    external :: perturbsetup

#if (defined SPMD)
    external :: MPI_ALLREDUCE
    external :: MPI_ALLGATHER
    external :: MPI_BCAST
    external :: MPI_ALLGATHERV
#endif

    allocate(HYMAP3_routing_struc(LIS_rc%nnest))

    HYMAP3_logunit = LIS_getNextUnitNumber()

    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%rslpmin  = 1e-5  !minimum slope
       HYMAP3_routing_struc(n)%nz      = 10    !number of stages in the sub-grid discretization
       HYMAP3_routing_struc(n)%imis     = -9999 !undefined integer value
       HYMAP3_routing_struc(n)%grv      = 9.81  !gravity accerelation [m/s2]
       HYMAP3_routing_struc(n)%numout   = 0
       HYMAP3_routing_struc(n)%fileopen = 0
       HYMAP3_routing_struc(n)%dt_proc  = 0.
    enddo

    !ag(14Apr2025)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 vector input file:",rc=status)
    do n=1, LIS_rc%nnest
       !by default, vecflag is zero, i.e., 2D domain
       HYMAP3_routing_struc(n)%vecflag=0
       call ESMF_ConfigGetAttribute(LIS_config, &
            HYMAP3_routing_struc(n)%vecfile,rc=status)
       if(status==0)then
          write(LIS_logunit,*) &
               '[INFO] HYMAP3 parameters in 1D vector mode'
          !if vector input file is provided, vecflag gets value 1
          HYMAP3_routing_struc(n)%vecflag=1
          !get HYMAP3 vector dimensions
          call HYMAP3_vector_read_dims(n)
          write(LIS_logunit,*) '[INFO] HYMAP3 vector input file: ',&
               trim(HYMAP3_routing_struc(n)%vecfile), &
               HYMAP3_routing_struc(n)%nz,HYMAP3_routing_struc(n)%nseqall
       else
          write(LIS_logunit,*) '[INFO] HYMAP3 parameters in 2D grid mode'
       endif
    enddo

    !ag (12Sep2019)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 enable 2-way coupling:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%enable2waycpl=0
       call ESMF_ConfigGetAttribute(LIS_config, &
            HYMAP3_routing_struc(n)%enable2waycpl,&
            default=0, rc=status)

       if (HYMAP3_routing_struc(n)%enable2waycpl==1) then
          write(LIS_logunit,*) '[INFO] HYMAP3 2-way coupling: activated'
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 2-way coupling flooded fraction threshold:", &
               rc=status)
          call ESMF_ConfigGetAttribute(LIS_config, &
               HYMAP3_routing_struc(n)%fldfrc2waycpl,rc=status)
          call LIS_verify(status,&
               "HYMAP3 2-way coupling flooded fraction threshold: not defined")
       else
          write(LIS_logunit,*) '[INFO] HYMAP3 2-way coupling: off'
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model time step:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP3 routing model time step: not defined")

       call LIS_parseTimeString(time,HYMAP3_routing_struc(n)%dt)

       call LIS_update_timestep(LIS_rc,n,HYMAP3_routing_struc(n)%dt)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP3 routing model output interval: not defined")

       call LIS_parseTimeString(time,HYMAP3_routing_struc(n)%outInterval)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 run in ensemble mode:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%useens,rc=status)
       call LIS_verify(status,&
            "HYMAP3 run in ensemble mode: not defined")
    enddo

    !Local inertia is the default routing method
    flowtype="local inertia"
    !ag (13Apr2016)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing method:",rc=status)
    if(status==0)call ESMF_ConfigGetAttribute(LIS_config,&
         flowtype,rc=status)
    do n=1, LIS_rc%nnest
       if(flowtype.eq."kinematic") then
          HYMAP3_routing_struc(n)%flowtype = 1
       elseif(flowtype.eq."diffusive") then
          HYMAP3_routing_struc(n)%flowtype = 2
       elseif(flowtype.eq."local inertia") then
          HYMAP3_routing_struc(n)%flowtype = 3
       elseif(flowtype.eq."hybrid") then
          HYMAP3_routing_struc(n)%flowtype = 0
          !ag(27Apr2020)
       elseif(flowtype.eq."urban") then
          HYMAP3_routing_struc(n)%flowtype = 4
       else
          write(LIS_logunit,*) &
               "[ERR] HYMAP3 routing method: unknown value ", &
               trim(flowtype)
          call LIS_endrun()
       endif
       write(LIS_logunit,*)"[INFO] HYMAP3 routing method: ",trim(flowtype)
    enddo

    !adaptive is the default time step method
    steptype="adaptive"
    !ag (11Sep2015)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model time step method:",rc=status)
    if(status==0)call ESMF_ConfigGetAttribute(LIS_config,&
         steptype,rc=status)
    do n=1, LIS_rc%nnest
       if(steptype.eq."constant") then
          HYMAP3_routing_struc(n)%steptype = 1
       elseif(steptype.eq."adaptive") then
          HYMAP3_routing_struc(n)%steptype = 2
       else
          write(LIS_logunit,*) &
               "[ERR] HYMAP3 routing model time step method: unknown value"
          call LIS_endrun()
       endif
       write(LIS_logunit,*) &
            "[INFO] HYMAP3 routing model time step method: ", &
            trim(steptype)
    enddo

    !ag (29Jan2016)
    if(steptype.eq."adaptive")then
       !0.5 is the default alpha coefficient
       cadp=0.5
       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP3 routing model adaptive time step alpha coefficient:",&
            rc=status)
       if(status==0)call ESMF_ConfigGetAttribute(LIS_config,&
            cadp,rc=status)
       do n=1, LIS_rc%nnest
          HYMAP3_routing_struc(n)%cadp=cadp
          write(LIS_logunit,*) &
               "[INFO] HYMAP3 adaptive alpha coefficient: ", &
               HYMAP3_routing_struc(n)%cadp
       enddo
    endif

    !ag (29Jun2016)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 floodplain dynamics:",rc=status)
    do n=1, LIS_rc%nnest
       !"1" is the default floodplain dynamics flag, i.e., "on"
       HYMAP3_routing_struc(n)%floodflag=1
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%floodflag,rc=status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model linear reservoir flag:",rc=status)
    do n=1, LIS_rc%nnest
       !"0" is the default linear reservoir flag, i.e., "off"
       HYMAP3_routing_struc(n)%linresflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%linresflag,rc=status)
    enddo

    !ag (24Apr2017)
    !"none" is the default evaporation option, i.e., "off"
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model evaporation option",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%evapflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%evapflag,rc=status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model dwi flag:",rc=status)
    do n=1, LIS_rc%nnest
       !"0" is the default linear dwi flag, i.e., "off"
       HYMAP3_routing_struc(n)%dwiflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%dwiflag,rc=status)
       write(LIS_logunit,*)"[INFO] HYMAP3 dwi flag: ", &
            HYMAP3_routing_struc(n)%dwiflag
    enddo

    !ag (24May2021)
    do n=1, LIS_rc%nnest
       if (HYMAP3_routing_struc(n)%enable2waycpl==1) then
          write(LIS_logunit,*) &
               "[INFO] 2-way coupling enabled. Disable linear reservoir and dwi"
          HYMAP3_routing_struc(n)%linresflag = 0
          HYMAP3_routing_struc(n)%dwiflag = 0
       endif
    enddo

    !ag (8Aug2020)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 discharge direct insertion:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%insertflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%insertflag,rc=status)
       if(HYMAP3_routing_struc(n)%insertflag==1)then
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 number of gauges:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%ninsert,rc=status)
          call LIS_verify(status,&
               "HYMAP3 number of gauges: not defined")

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 discharge direct insertion time series length:", &
               rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%ntinsert,rc=status)
          call LIS_verify(status,&
               "HYMAP3 discharge direct insertion time series length: not defined")

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 discharge direct insertion input directory:", &
               rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%insertdir,rc=status)
          call LIS_verify(status,&
               "HYMAP3 discharge direct insertion input directory: not defined")

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 discharge direct insertion header filename:", &
               rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%insertheader,rc=status)
          call LIS_verify(status,&
               "HYMAP3 discharge direct insertion header filename: not defined")
       endif
    enddo

    !ag (30Apr2021)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 water management:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%managflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%managflag,rc=status)
       if(HYMAP3_routing_struc(n)%managflag==1)then
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 water management header file:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%managheader,rc=status)
          call LIS_verify(status,&
               "HYMAP3 water management header file: not defined")
       endif
    enddo

    !ag (30Mar2021)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 sea level:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%sealevelflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%sealevelflag,rc=status)
       if(HYMAP3_routing_struc(n)%sealevelflag==1)then
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 sea level input directory:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%sealeveldir,rc=status)
          call LIS_verify(status,&
               "HYMAP3 sea level input directory: not defined")

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 sea level header filename:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%sealevelheader,rc=status)
          call LIS_verify(status,&
               "HYMAP3 sea level header filename: not defined")

          call HYMAP3_read_header_size1( &
               trim(HYMAP3_routing_struc(n)%sealevelheader),&
               HYMAP3_routing_struc(n)%nsealevel,&
               HYMAP3_routing_struc(n)%ntsealevel)

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 sea level outlet list:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%outletlist,rc=status)
          call LIS_verify(status,&
               "HYMAP3 sea level outlet list: not defined")

          call HYMAP3_read_header_size( &
               trim(HYMAP3_routing_struc(n)%outletlist),&
               HYMAP3_routing_struc(n)%noutlet)

          allocate(HYMAP3_routing_struc(n)%tsealevel( &
               HYMAP3_routing_struc(n)%nsealevel, &
               HYMAP3_routing_struc(n)%ntsealevel))
          allocate(HYMAP3_routing_struc(n)%sealevel( &
               HYMAP3_routing_struc(n)%nsealevel, &
               HYMAP3_routing_struc(n)%ntsealevel))

          HYMAP3_routing_struc(n)%tsealevel= &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%sealevel= &
               real(HYMAP3_routing_struc(n)%imis)
       endif
    enddo

    !ag(30Mar2022)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 river bifurcation:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%bifflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%bifflag,rc=status)
       if(HYMAP3_routing_struc(n)%bifflag==1)then
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 river bifurcation pathway file:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%biffile,rc=status)
          call LIS_verify(status,&
               "HYMAP3 river bifurcation pathway file: not defined")
          !ag(27Jul2025)
          !get number of bifurcations and elevations
          call HYMAP3_read_header_size1(HYMAP3_routing_struc(n)%biffile,&
               HYMAP3_routing_struc(n)%nbif, &
               HYMAP3_routing_struc(n)%nbifelv)
       endif
    enddo

    !ag(23Feb2023)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 levee:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%levflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%levflag,rc=status)
    enddo

    !ag (4Feb2016)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 reservoir operation option:",rc=status)
    do n=1, LIS_rc%nnest
       HYMAP3_routing_struc(n)%resopflag=0
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%resopflag,default=0,rc=status)

       if(HYMAP3_routing_struc(n)%resopflag==1)then

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 reservoir operation input directory:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%resopdir,rc=status)
          call LIS_verify(status,&
               "HYMAP3 reservoir operation input directory: not defined")

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 reservoir operation header filename:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%resopheader,rc=status)
          call LIS_verify(status,&
               "HYMAP3 reservoir operation header filename: not defined")
          !ag(27Jul2025)
          call HYMAP3_read_header_size1( &
               trim(HYMAP3_routing_struc(n)%resopheader),&
               HYMAP3_routing_struc(n)%nresop, &
               HYMAP3_routing_struc(n)%ntresop)
          !ag(1Apr2025)
          !Yassin's scheme
       elseif(HYMAP3_routing_struc(n)%resopflag==2)then
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 reservoir operation input file:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%resopncfile,rc=status)
          call LIS_verify(status,&
               "HYMAP3 reservoir operation input file: not defined")

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 reservoir list:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%resopheader,rc=status)
          call LIS_verify(status,&
               "HYMAP3 reservoir list: not defined")
          !ag(27Jul2025)
          call HYMAP3_read_header_size( &
               trim(HYMAP3_routing_struc(n)%resopheader),&
               HYMAP3_routing_struc(n)%nresop)
       endif
    enddo

    write(LIS_logunit,*) '[INFO] Initializing HYMAP3....'
    !allocate matrixes
    do n=1, LIS_rc%nnest
       !ag(14Apr2025)
       if(HYMAP3_routing_struc(n)%vecflag==0)then
          write(LIS_logunit,*)'[INFO] columns and rows', &
               LIS_rc%lnc(n),LIS_rc%lnr(n)

          !ag (19Feb2016)
          allocate(HYMAP3_routing_struc(n)%nextx( &
               LIS_rc%gnc(n),LIS_rc%gnr(n)))
          ctitle = 'HYMAP_flow_direction_x'
          call HYMAP3_read_param_int_2d_global(ctitle,n, &
               HYMAP3_routing_struc(n)%nextx)

          allocate(HYMAP3_routing_struc(n)%nexty( &
               LIS_rc%gnc(n),LIS_rc%gnr(n)))
          ctitle = 'HYMAP_flow_direction_y'
          call HYMAP3_read_param_int_2d_global(ctitle,n, &
               HYMAP3_routing_struc(n)%nexty)

          allocate(elevtn(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          ctitle = 'HYMAP_grid_elevation'
          call HYMAP3_read_param_real_2d(ctitle,n,elevtn)

          allocate(uparea(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          ctitle = 'HYMAP_drain_area'
          call HYMAP3_read_param_real_2d(ctitle,n,uparea)

          allocate(basin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          ctitle = 'HYMAP_basin'
          call HYMAP3_read_param_real_2d(ctitle,n,basin)

          allocate(mask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          ctitle = 'HYMAP_basin_mask'
          call HYMAP3_read_param_int_2d(ctitle,n,mask)

          allocate(maskg(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          ctitle = 'HYMAP_basin_mask'
          call HYMAP3_read_param_int_2d_global(ctitle,n,maskg)

          !Assign the mask to the routing data structure
          LIS_routing(n)%dommask = mask
          LIS_routing(n)%nextx = HYMAP3_routing_struc(n)%nextx

          write(LIS_logunit,*) '[INFO] Get number of HYMAP3 grid cells'
          call HYMAP3_get_vector_size(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               LIS_rc%gnc(n),LIS_rc%gnr(n),&
               LIS_ews_halo_ind(n,LIS_localPet+1), &
               LIS_nss_halo_ind(n,LIS_localPet+1), &
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%nextx,&
               mask,HYMAP3_routing_struc(n)%nseqall)

       elseif(HYMAP3_routing_struc(n)%vecflag==1)then
          !there's nothing to be done here
       endif

       LIS_rc%nroutinggrid(n) = HYMAP3_routing_struc(n)%nseqall
       gdeltas = HYMAP3_routing_struc(n)%nseqall

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
                  LIS_routing_goffsets(n,i-1) + &
                  LIS_routing_gdeltas(n,i-1)
          enddo
       endif
#if (defined SPMD)
       call MPI_BCAST(LIS_routing_goffsets(n,:), &
            LIS_npes, MPI_INTEGER,0, &
            LIS_mpi_comm, status)
#endif

       allocate(HYMAP3_routing_struc(n)%seqx( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%seqy( &
            HYMAP3_routing_struc(n)%nseqall))

       allocate(HYMAP3_routing_struc(n)%seqx_glb( &
            LIS_rc%glbnroutinggrid(n)))
       allocate(HYMAP3_routing_struc(n)%seqy_glb( &
            LIS_rc%glbnroutinggrid(n)))

       allocate(HYMAP3_routing_struc(n)%sindex( &
            LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(HYMAP3_routing_struc(n)%outlet( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%outlet_glb( &
            LIS_rc%glbnroutinggrid(n)))
       allocate(HYMAP3_routing_struc(n)%next( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%next_glb( &
            LIS_rc%glbnroutinggrid(n)))

       allocate(HYMAP3_routing_struc(n)%elevtn( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%nxtdst( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%grarea(&
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%fldman(&
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%fldhgt(&
            HYMAP3_routing_struc(n)%nseqall,&
            HYMAP3_routing_struc(n)%nz))
       allocate(HYMAP3_routing_struc(n)%cntime( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%trnoff(&
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%tbsflw(&
            HYMAP3_routing_struc(n)%nseqall))

       if(HYMAP3_routing_struc(n)%useens.eq.0) then
          allocate(HYMAP3_routing_struc(n)%rivman( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivwth( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivhgt( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivstomax( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivare( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldstomax( &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%nz,1))
          allocate(HYMAP3_routing_struc(n)%fldgrd( &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%nz,1))
       else
          allocate(HYMAP3_routing_struc(n)%rivman( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivwth( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivhgt( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivstomax( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivare( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldstomax( &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%nz,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldgrd( &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%nz,LIS_rc%nensem(n)))
       endif

       allocate(HYMAP3_routing_struc(n)%rivelv( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%rivlen( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%runoff0( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%basflw0( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%flowmap( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%rnfdwi_ratio( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%bsfdwi_ratio( &
            HYMAP3_routing_struc(n)%nseqall))

       !ag (4Feb2016)
       if(HYMAP3_routing_struc(n)%resopflag==1)then
          allocate(HYMAP3_routing_struc(n)%resoploc( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%resoploc_dwn( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%resoploc_glb( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%resoploc_dwn_glb( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%resopoutmin( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%tresopalt( &
               HYMAP3_routing_struc(n)%nresop, &
               HYMAP3_routing_struc(n)%ntresop))
          allocate(HYMAP3_routing_struc(n)%resopalt( &
               HYMAP3_routing_struc(n)%nresop, &
               HYMAP3_routing_struc(n)%ntresop))
          allocate(HYMAP3_routing_struc(n)%resoptype( &
               HYMAP3_routing_struc(n)%nresop))

          !ag(11Oct2024)
          allocate(HYMAP3_routing_struc(n)%elevtn_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%fldhgt_resop( &
               HYMAP3_routing_struc(n)%nresop,HYMAP3_routing_struc(n)%nz))
          allocate(HYMAP3_routing_struc(n)%fldstomax_resop( &
               HYMAP3_routing_struc(n)%nresop,HYMAP3_routing_struc(n)%nz))
          allocate(HYMAP3_routing_struc(n)%grarea_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%rivstomax_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%rivelv_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%rivlen_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%rivwth_resop( &
               HYMAP3_routing_struc(n)%nresop))

          !ag(4Apr2025)
          !Yassin's scheme
       elseif(HYMAP3_routing_struc(n)%resopflag==2)then
          allocate(HYMAP3_routing_struc(n)%resoploc( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%resoploc_dwn( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%resoploc_glb( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%resoploc_dwn_glb( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%ncloc_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%maxsto_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%inidis_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%inisto_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%dwndis_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%deadis_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%minsto_mo_resop( &
               HYMAP3_routing_struc(n)%nresop,12))
          allocate(HYMAP3_routing_struc(n)%nupsto_mo_resop( &
               HYMAP3_routing_struc(n)%nresop,12))
          allocate(HYMAP3_routing_struc(n)%uppsto_mo_resop( &
               HYMAP3_routing_struc(n)%nresop,12))
          allocate(HYMAP3_routing_struc(n)%mindis_mo_resop( &
               HYMAP3_routing_struc(n)%nresop,12))
          allocate(HYMAP3_routing_struc(n)%nupdis_mo_resop( &
               HYMAP3_routing_struc(n)%nresop,12))
          allocate(HYMAP3_routing_struc(n)%uppdis_mo_resop( &
               HYMAP3_routing_struc(n)%nresop,12))
          allocate(HYMAP3_routing_struc(n)%reg1_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%reg2_resop( &
               HYMAP3_routing_struc(n)%nresop))
          allocate(HYMAP3_routing_struc(n)%reg3_resop( &
               HYMAP3_routing_struc(n)%nresop))
       endif

       !ag(27Apr2020)
       allocate(HYMAP3_routing_struc(n)%drstomax( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%droutlet( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%droutlet_glb( &
            LIS_rc%glbnroutinggrid(n)))
       allocate(HYMAP3_routing_struc(n)%drtotwth( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%drnoutlet( &
            HYMAP3_routing_struc(n)%nseqall))
       allocate(HYMAP3_routing_struc(n)%drtotlgh( &
            HYMAP3_routing_struc(n)%nseqall))

       !ag(8Aug2020)
       if(HYMAP3_routing_struc(n)%insertflag==1)then
          allocate(HYMAP3_routing_struc(n)%insertloc( &
               HYMAP3_routing_struc(n)%ninsert))
          allocate(HYMAP3_routing_struc(n)%tinsert( &
               HYMAP3_routing_struc(n)%ninsert, &
               HYMAP3_routing_struc(n)%ntinsert))
          allocate(HYMAP3_routing_struc(n)%insertdis( &
               HYMAP3_routing_struc(n)%ninsert, &
               HYMAP3_routing_struc(n)%ntinsert))
       endif

       !ag(23Feb2023)
       allocate(HYMAP3_routing_struc(n)%levhgt( &
            HYMAP3_routing_struc(n)%nseqall))

       if(HYMAP3_routing_struc(n)%useens.eq.0) then
          allocate(HYMAP3_routing_struc(n)%rivsto( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivdph( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivvel( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivout( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%evpout( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldout( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldsto( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%flddph( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldvel( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldfrc( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldare( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%sfcelv( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rnfsto( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%bsfsto( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rnfdwi( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%bsfdwi( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%surfws( &
               HYMAP3_routing_struc(n)%nseqall,1))

          allocate(HYMAP3_routing_struc(n)%dtaout( &
               HYMAP3_routing_struc(n)%nseqall,1))
          !ag (19Jan2016)
          allocate(HYMAP3_routing_struc(n)%rivout_pre( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%rivdph_pre( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldout_pre( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%flddph_pre( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldelv1( &
               HYMAP3_routing_struc(n)%nseqall,1))
          !ag (03May2017)
          allocate(HYMAP3_routing_struc(n)%ewat( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%edif( &
               HYMAP3_routing_struc(n)%nseqall,1))
          !ag (12Sep2019)
          allocate(HYMAP3_routing_struc(n)%rivstotmp( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldstotmp( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldfrctmp( &
               HYMAP3_routing_struc(n)%nseqall,1))
          !ag(27Apr2020)
          allocate(HYMAP3_routing_struc(n)%drsto( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%drout( &
               HYMAP3_routing_struc(n)%nseqall,1))

          !ag(22May2024)
          allocate(HYMAP3_routing_struc(n)%levstomax( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldstoatlev( &
               HYMAP3_routing_struc(n)%nseqall,1))
          allocate(HYMAP3_routing_struc(n)%fldonlystomax( &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%nz,1))
       else
          allocate(HYMAP3_routing_struc(n)%rivsto( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivdph( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivvel( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivout( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%evpout( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldout( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldsto( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%flddph( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldvel( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldfrc( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldare( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%sfcelv( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rnfsto( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%bsfsto( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rnfdwi( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%bsfdwi( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%surfws( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%dtaout( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          !ag (19Jan2016)
          allocate(HYMAP3_routing_struc(n)%rivout_pre( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%rivdph_pre( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldout_pre( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%flddph_pre( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldelv1( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          !ag (03May2017)
          allocate(HYMAP3_routing_struc(n)%ewat( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%edif( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          !ag (12Sep2019)
          allocate(HYMAP3_routing_struc(n)%rivstotmp( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldstotmp( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldfrctmp( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          !ag(27Apr2020)
          allocate(HYMAP3_routing_struc(n)%drsto( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%drout( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          !ag(22May2024)
          allocate(HYMAP3_routing_struc(n)%levstomax( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldstoatlev( &
               HYMAP3_routing_struc(n)%nseqall,LIS_rc%nensem(n)))
          allocate(HYMAP3_routing_struc(n)%fldonlystomax( &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%nz,LIS_rc%nensem(n)))
       endif
       !ag(9Oct2024)
       allocate(HYMAP3_routing_struc(n)%outletid( &
            HYMAP3_routing_struc(n)%nseqall))
       HYMAP3_routing_struc(n)%outletid=HYMAP3_routing_struc(n)%imis

       HYMAP3_routing_struc(n)%seqx=0.0
       HYMAP3_routing_struc(n)%seqy=0.0
       HYMAP3_routing_struc(n)%sindex=0.0
       HYMAP3_routing_struc(n)%outlet=0.0
       HYMAP3_routing_struc(n)%outlet_glb=0.0
       HYMAP3_routing_struc(n)%next=0.0
       HYMAP3_routing_struc(n)%next_glb=0.0
       HYMAP3_routing_struc(n)%elevtn=0.0
       HYMAP3_routing_struc(n)%nxtdst=0.0
       HYMAP3_routing_struc(n)%grarea=0.0
       HYMAP3_routing_struc(n)%fldgrd=0.0
       HYMAP3_routing_struc(n)%fldman=0.0
       HYMAP3_routing_struc(n)%fldhgt=0.0
       HYMAP3_routing_struc(n)%cntime=0.0
       HYMAP3_routing_struc(n)%trnoff=0.0
       HYMAP3_routing_struc(n)%tbsflw=0.0
       HYMAP3_routing_struc(n)%fldstomax=0.0
       HYMAP3_routing_struc(n)%rivman=0.0
       HYMAP3_routing_struc(n)%rivelv=0.0
       HYMAP3_routing_struc(n)%rivstomax=0.0
       HYMAP3_routing_struc(n)%rivlen=0.0
       HYMAP3_routing_struc(n)%rivwth=0.0
       HYMAP3_routing_struc(n)%rivhgt=0.0
       HYMAP3_routing_struc(n)%rivare=0.0
       HYMAP3_routing_struc(n)%runoff0=0.0
       HYMAP3_routing_struc(n)%basflw0=0.0
       HYMAP3_routing_struc(n)%rivsto=0.0
       HYMAP3_routing_struc(n)%rivdph=0.0
       HYMAP3_routing_struc(n)%rivvel=0.0
       HYMAP3_routing_struc(n)%rivout=0.0
       HYMAP3_routing_struc(n)%evpout=0.0
       HYMAP3_routing_struc(n)%fldout=0.0
       HYMAP3_routing_struc(n)%fldsto=0.0
       HYMAP3_routing_struc(n)%flddph=0.0
       HYMAP3_routing_struc(n)%fldvel=0.0
       HYMAP3_routing_struc(n)%fldfrc=0.0
       HYMAP3_routing_struc(n)%fldare=0.0
       HYMAP3_routing_struc(n)%sfcelv=0.0
       HYMAP3_routing_struc(n)%rnfsto=0.0
       HYMAP3_routing_struc(n)%bsfsto=0.0
       HYMAP3_routing_struc(n)%flowmap=0.0
       HYMAP3_routing_struc(n)%rnfdwi_ratio=0.0
       HYMAP3_routing_struc(n)%bsfdwi_ratio=0.0
       HYMAP3_routing_struc(n)%rnfdwi=0.0
       HYMAP3_routing_struc(n)%bsfdwi=0.0
       HYMAP3_routing_struc(n)%surfws=0.0
       HYMAP3_routing_struc(n)%dtaout=1000000.

       !ag (19Jan2016)
       HYMAP3_routing_struc(n)%rivout_pre=0.0
       HYMAP3_routing_struc(n)%rivdph_pre=0.0
       HYMAP3_routing_struc(n)%fldout_pre=0.0
       HYMAP3_routing_struc(n)%flddph_pre=0.0
       HYMAP3_routing_struc(n)%fldelv1=0.0

       !ag (4Feb2016)
       if(HYMAP3_routing_struc(n)%resopflag==1)then
          HYMAP3_routing_struc(n)%resoploc = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%resoploc_dwn = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%resoploc_glb = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%resoploc_dwn_glb = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%resopoutmin = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%tresopalt = &
               dble(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%resopalt = &
               real(HYMAP3_routing_struc(n)%imis)

          !ag(11Oct2024)
          HYMAP3_routing_struc(n)%elevtn_resop = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%fldhgt_resop = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%fldstomax_resop = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%grarea_resop = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%rivstomax_resop = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%rivelv_resop = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%rivlen_resop = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%rivwth_resop = &
               real(HYMAP3_routing_struc(n)%imis)
       endif

       !ag (03May2017)
       HYMAP3_routing_struc(n)%ewat=0.0
       HYMAP3_routing_struc(n)%edif=0.0

       !ag(27Apr2020)
       HYMAP3_routing_struc(n)%drsto=0.0
       HYMAP3_routing_struc(n)%drout=0.0

       !ag(8Aug2020)
       if(HYMAP3_routing_struc(n)%insertflag==1)then
          HYMAP3_routing_struc(n)%insertloc = &
               real(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%tinsert = &
               dble(HYMAP3_routing_struc(n)%imis)
          HYMAP3_routing_struc(n)%insertdis = &
               real(HYMAP3_routing_struc(n)%imis)
       endif

       !ag(22May2024)
       HYMAP3_routing_struc(n)%levstomax=0.0
       HYMAP3_routing_struc(n)%fldstoatlev=0.0
       HYMAP3_routing_struc(n)%fldonlystomax=0.0
    enddo

    do n=1, LIS_rc%nnest
       !ag(15May2025)
       !read grid-based input parameters
       if(HYMAP3_routing_struc(n)%vecflag==0)then
          write(LIS_logunit,*) '[INFO] Get HYMAP3 cell vector sequence'
          call HYMAP3_get_sindex(LIS_rc%gnc(n),&
               LIS_rc%gnr(n),&
               LIS_rc%glbnroutinggrid(n),&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%nextx,&
               HYMAP3_routing_struc(n)%nexty,maskg,&
               HYMAP3_routing_struc(n)%sindex,&
               HYMAP3_routing_struc(n)%outlet_glb,&
               HYMAP3_routing_struc(n)%next_glb)

          call HYMAP3_get_seq(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               LIS_rc%gnc(n),LIS_rc%gnr(n),&
               LIS_ews_halo_ind(n,LIS_localPet+1), &
               LIS_nss_halo_ind(n,LIS_localPet+1), &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%nextx,&
               HYMAP3_routing_struc(n)%nexty,maskg,&
               HYMAP3_routing_struc(n)%sindex,&
               HYMAP3_routing_struc(n)%outlet,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,HYMAP3_routing_struc(n)%next)

#if (defined SPMD)
          call MPI_ALLGATHERV(HYMAP3_routing_struc(n)%seqx,&
               LIS_rc%nroutinggrid(n),MPI_INTEGER,&
               HYMAP3_routing_struc(n)%seqx_glb,&
               LIS_routing_gdeltas(n,:),&
               LIS_routing_goffsets(n,:),&
               MPI_INTEGER,&
               LIS_mpi_comm,status)

          call MPI_ALLGATHERV(HYMAP3_routing_struc(n)%seqy,&
               LIS_rc%nroutinggrid(n),MPI_INTEGER,&
               HYMAP3_routing_struc(n)%seqy_glb,&
               LIS_routing_gdeltas(n,:),&
               LIS_routing_goffsets(n,:),&
               MPI_INTEGER,&
               LIS_mpi_comm,status)
#endif

          allocate(tmp_real(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(tmp_real_nz(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               HYMAP3_routing_struc(n)%nz))

          ctitle = 'HYMAP_river_width'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          if(HYMAP3_routing_struc(n)%useens.eq.0) then
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%rivwth(:,1))
          else
             do m=1,LIS_rc%nensem(n)
                call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                     1,HYMAP3_routing_struc(n)%nseqall,&
                     HYMAP3_routing_struc(n)%imis,&
                     HYMAP3_routing_struc(n)%seqx,&
                     HYMAP3_routing_struc(n)%seqy,&
                     tmp_real,HYMAP3_routing_struc(n)%rivwth(:,m))
             enddo
          endif

          ctitle = 'HYMAP_river_length'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%rivlen)

          ctitle = 'HYMAP_river_height'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          if(HYMAP3_routing_struc(n)%useens.eq.0) then
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%rivhgt(:,1))
          else
             do m=1,LIS_rc%nensem(n)
                call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                     1,HYMAP3_routing_struc(n)%nseqall,&
                     HYMAP3_routing_struc(n)%imis,&
                     HYMAP3_routing_struc(n)%seqx,&
                     HYMAP3_routing_struc(n)%seqy,&
                     tmp_real,HYMAP3_routing_struc(n)%rivhgt(:,m))
             enddo
          endif

          ctitle = 'HYMAP_river_roughness'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          if(HYMAP3_routing_struc(n)%useens.eq.0) then
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%rivman(:,1))
          else
             do m=1,LIS_rc%nensem(n)
                call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                     1,HYMAP3_routing_struc(n)%nseqall,&
                     HYMAP3_routing_struc(n)%imis,&
                     HYMAP3_routing_struc(n)%seqx,&
                     HYMAP3_routing_struc(n)%seqy,&
                     tmp_real,HYMAP3_routing_struc(n)%rivman(:,m))
             enddo
          endif

          ctitle = 'HYMAP_floodplain_height'
          call HYMAP3_read_param_real(ctitle,n, &
               HYMAP3_routing_struc(n)%nz,&
               tmp_real_nz)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               HYMAP3_routing_struc(n)%nz,&
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real_nz,HYMAP3_routing_struc(n)%fldhgt)

          ctitle = 'HYMAP_floodplain_roughness'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%fldman)

          ctitle = 'HYMAP_grid_elevation'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%elevtn)

          ctitle = 'HYMAP_grid_distance'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%nxtdst)

          ctitle = 'HYMAP_grid_area'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%grarea)

          ctitle = 'HYMAP_runoff_time_delay'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%cntime)
          where(HYMAP3_routing_struc(n)%cntime==0) &
               HYMAP3_routing_struc(n)%cntime = &
               minval(HYMAP3_routing_struc(n)%cntime, &
               HYMAP3_routing_struc(n)%cntime>0)

          ctitle = 'HYMAP_runoff_time_delay_multiplier'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%trnoff)

          ctitle = 'HYMAP_baseflow_time_delay'
          call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
          call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
               1,HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%imis,&
               HYMAP3_routing_struc(n)%seqx,&
               HYMAP3_routing_struc(n)%seqy,&
               tmp_real,HYMAP3_routing_struc(n)%tbsflw)
          write(LIS_logunit,*) &
               '[INFO] Remove zeros from groundwater time delay'
          where(HYMAP3_routing_struc(n)%tbsflw < &
               HYMAP3_routing_struc(n)%cntime/86400..and.&
               HYMAP3_routing_struc(n)%cntime>0) &
               HYMAP3_routing_struc(n)%tbsflw = &
               HYMAP3_routing_struc(n)%cntime/86400.

          !ag (23Nov2016)
          ctitle = 'HYMAP_runoff_dwi_ratio'
          if(HYMAP3_routing_struc(n)%dwiflag==1)then
             call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%rnfdwi_ratio)
          else
             HYMAP3_routing_struc(n)%rnfdwi_ratio=0.
          endif

          !ag (23Nov2016)
          ctitle = 'HYMAP_baseflow_dwi_ratio'
          if(HYMAP3_routing_struc(n)%dwiflag==1)then
             call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%bsfdwi_ratio)
          else
             HYMAP3_routing_struc(n)%bsfdwi_ratio=0.
          endif

          !ag (13Apr2016)
          ctitle = 'HYMAP_river_flow_type'
          !ag(27Apr2020)
          if(HYMAP3_routing_struc(n)%flowtype==0 .or. &
               HYMAP3_routing_struc(n)%flowtype==4)then
             call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%flowmap)
          else
             HYMAP3_routing_struc(n)%flowmap = &
                  HYMAP3_routing_struc(n)%flowtype
          endif

          !ag (7Dec2020)
          ctitle = 'HYMAP_urban_drainage_outlet'
          if(HYMAP3_routing_struc(n)%flowtype==4)then
             call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%droutlet)
          else
             HYMAP3_routing_struc(n)%droutlet = &
                  HYMAP3_routing_struc(n)%imis
          endif

          !ag(23Feb2023)
          ctitle = 'HYMAP_levee_height'
          if(HYMAP3_routing_struc(n)%levflag==1)then
             call HYMAP3_read_param_real_2d(ctitle,n,tmp_real)
             call HYMAP3_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
                  1,HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%imis,&
                  HYMAP3_routing_struc(n)%seqx,&
                  HYMAP3_routing_struc(n)%seqy,&
                  tmp_real,HYMAP3_routing_struc(n)%levhgt)
          else
             HYMAP3_routing_struc(n)%levhgt=0.0
          endif

#if (defined SPMD)
          call MPI_ALLGATHERV(HYMAP3_routing_struc(n)%droutlet,&
               LIS_rc%nroutinggrid(n),MPI_INTEGER,&
               HYMAP3_routing_struc(n)%droutlet_glb,&
               LIS_routing_gdeltas(n,:),&
               LIS_routing_goffsets(n,:),&
               MPI_INTEGER,&
               LIS_mpi_comm,status)
#endif
          deallocate(tmp_real)
          deallocate(tmp_real_nz)

          !read vectorized input parameters
       elseif(HYMAP3_routing_struc(n)%vecflag==1)then
          ctitle = 'HYMAP_river_width'
          if(HYMAP3_routing_struc(n)%useens.eq.0) then
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%rivwth(:,1))
          else
             do m=1,LIS_rc%nensem(n)
                call HYMAP3_vector_read_param(ctitle,n,1, &
                     HYMAP3_routing_struc(n)%rivwth(:,m))
             enddo
          endif
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%rivwth)

          ctitle = 'HYMAP_river_height'
          if(HYMAP3_routing_struc(n)%useens.eq.0) then
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%rivhgt(:,1))
          else
             do m=1,LIS_rc%nensem(n)
                call HYMAP3_vector_read_param(ctitle,n,1, &
                     HYMAP3_routing_struc(n)%rivhgt(:,m))
             enddo
          endif
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%rivhgt)

          ctitle = 'HYMAP_river_roughness'
          if(HYMAP3_routing_struc(n)%useens.eq.0) then
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%rivman(:,1))
          else
             do m=1,LIS_rc%nensem(n)
                call HYMAP3_vector_read_param(ctitle,n,1, &
                     HYMAP3_routing_struc(n)%rivman(:,m))
             enddo
          endif
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%rivman)

          ctitle = 'HYMAP_river_length'
          call HYMAP3_vector_read_param(ctitle,n,1, &
               HYMAP3_routing_struc(n)%rivlen)
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%rivlen)

          ctitle = 'HYMAP_floodplain_height'
          call HYMAP3_vector_read_param(ctitle,n, &
               HYMAP3_routing_struc(n)%nz,&
               HYMAP3_routing_struc(n)%fldhgt)
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%fldhgt)

          ctitle = 'HYMAP_floodplain_roughness'
          call HYMAP3_vector_read_param(ctitle,n,1, &
               HYMAP3_routing_struc(n)%fldman)
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%fldman)

          ctitle = 'HYMAP_grid_elevation'
          call HYMAP3_vector_read_param(ctitle,n,1, &
               HYMAP3_routing_struc(n)%elevtn)
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%elevtn)

          ctitle = 'HYMAP_grid_distance'
          call HYMAP3_vector_read_param(ctitle,n,1, &
               HYMAP3_routing_struc(n)%nxtdst)
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%nxtdst)

          ctitle = 'HYMAP_grid_area'
          call HYMAP3_vector_read_param(ctitle,n,1, &
               HYMAP3_routing_struc(n)%grarea)
          write(LIS_logunit,*) trim(ctitle), &
               maxval(HYMAP3_routing_struc(n)%grarea)

          if(HYMAP3_routing_struc(n)%linresflag==1)then
             ctitle = 'HYMAP_runoff_time_delay'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%cntime)
             where(HYMAP3_routing_struc(n)%cntime==0) &
                  HYMAP3_routing_struc(n)%cntime = &
                  minval(HYMAP3_routing_struc(n)%cntime, &
                  HYMAP3_routing_struc(n)%cntime>0)

             ctitle = 'HYMAP_runoff_time_delay_multiplier'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%trnoff)

             ctitle = 'HYMAP_baseflow_time_delay'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%tbsflw)
             write(LIS_logunit,*) &
                  '[INFO] Remove zeros from groundwater time delay'
             where(HYMAP3_routing_struc(n)%tbsflw < &
                  HYMAP3_routing_struc(n)%cntime/86400. .and. &
                  HYMAP3_routing_struc(n)%cntime>0) &
                  HYMAP3_routing_struc(n)%tbsflw = &
                  HYMAP3_routing_struc(n)%cntime/86400.
          endif

          if(HYMAP3_routing_struc(n)%dwiflag==1)then
             ctitle = 'HYMAP_baseflow_dwi_ratio'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%rnfdwi_ratio)
             ctitle = 'HYMAP_runoff_dwi_ratio'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%bsfdwi_ratio)
          else
             HYMAP3_routing_struc(n)%rnfdwi_ratio=0.
             HYMAP3_routing_struc(n)%bsfdwi_ratio=0.
          endif

          if(HYMAP3_routing_struc(n)%flowtype==0.or. &
               HYMAP3_routing_struc(n)%flowtype==4)then
             ctitle = 'HYMAP_river_flow_type'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%flowmap)
          else
             HYMAP3_routing_struc(n)%flowmap = &
                  HYMAP3_routing_struc(n)%flowtype
          endif

          if(HYMAP3_routing_struc(n)%flowtype==4)then
             ctitle = 'HYMAP_urban_drainage_outlet'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%droutlet)
          else
             HYMAP3_routing_struc(n)%droutlet = &
                  HYMAP3_routing_struc(n)%imis
          endif

          if(HYMAP3_routing_struc(n)%levflag==1)then
             ctitle = 'HYMAP_levee_height'
             call HYMAP3_vector_read_param(ctitle,n,1, &
                  HYMAP3_routing_struc(n)%levhgt)
          else
             HYMAP3_routing_struc(n)%levhgt=0.
          endif
       endif
    enddo

    !ag (20Sep2016) Correction for cases where parameter maps don't match
    do n=1, LIS_rc%nnest
       if(HYMAP3_routing_struc(n)%useens.eq.0) then
          do i=1,HYMAP3_routing_struc(n)%nseqall
             if(HYMAP3_routing_struc(n)%rivwth(i,1) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%rivlen(i) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%rivhgt(i,1) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%rivman(i,1) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%fldhgt(i,1) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%elevtn(i) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%nxtdst(i) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%grarea(i) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%cntime(i) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%trnoff(i) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%tbsflw(i) == &
                  HYMAP3_routing_struc(n)%imis.or. &
                  HYMAP3_routing_struc(n)%flowmap(i) == &
                  HYMAP3_routing_struc(n)%imis) then

                HYMAP3_routing_struc(n)%rivwth(i,1) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%rivlen(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%rivhgt(i,1) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%rivman(i,1) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%fldhgt(i,:) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%elevtn(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%nxtdst(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%grarea(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%cntime(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%trnoff(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%tbsflw(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%flowmap(i) = &
                     HYMAP3_routing_struc(n)%imis
                HYMAP3_routing_struc(n)%outlet(i) = &
                     HYMAP3_routing_struc(n)%imis
             endif
          enddo
       else
          do i=1,HYMAP3_routing_struc(n)%nseqall
             do m=1,LIS_rc%nensem(n)
                if(HYMAP3_routing_struc(n)%rivwth(i,m) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%rivlen(i) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%rivhgt(i,m) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%rivman(i,m) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%fldhgt(i,1) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%elevtn(i) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%nxtdst(i) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%grarea(i) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%cntime(i) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%trnoff(i) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%tbsflw(i) == &
                     HYMAP3_routing_struc(n)%imis.or. &
                     HYMAP3_routing_struc(n)%flowmap(i) == &
                     HYMAP3_routing_struc(n)%imis) then

                   HYMAP3_routing_struc(n)%rivwth(i,m) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%rivlen(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%rivhgt(i,m) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%rivman(i,m) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%fldhgt(i,:) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%elevtn(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%nxtdst(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%grarea(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%cntime(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%trnoff(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%tbsflw(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%flowmap(i) = &
                        HYMAP3_routing_struc(n)%imis
                   HYMAP3_routing_struc(n)%outlet(i) = &
                        HYMAP3_routing_struc(n)%imis
                endif
             enddo
          enddo
       endif
    enddo

    write(LIS_logunit,*) '[INFO] Processing data before running HYMAP'
    do n=1, LIS_rc%nnest
       write(LIS_logunit,*)'[INFO] Calculate maximum river storage'

       if(HYMAP3_routing_struc(n)%useens.eq.0) then
          do i=1,HYMAP3_routing_struc(n)%nseqall
             HYMAP3_routing_struc(n)%rivstomax(i,1) = &
                  HYMAP3_routing_struc(n)%rivlen(i) * &
                  HYMAP3_routing_struc(n)%rivwth(i,1) * &
                  HYMAP3_routing_struc(n)%rivhgt(i,1)
             !ag(23Feb2023) (22May2024)
             HYMAP3_routing_struc(n)%levstomax(i,1) = &
                  HYMAP3_routing_struc(n)%rivlen(i) * &
                  HYMAP3_routing_struc(n)%rivwth(i,1) * &
                  HYMAP3_routing_struc(n)%levhgt(i)
             HYMAP3_routing_struc(n)%rivelv(i) = &
                  HYMAP3_routing_struc(n)%elevtn(i) - &
                  HYMAP3_routing_struc(n)%rivhgt(i,1)
             if(HYMAP3_routing_struc(n)%rivwth(i,1)>0) then
                HYMAP3_routing_struc(n)%rivare(i,1) = &
                     min(HYMAP3_routing_struc(n)%grarea(i), &
                     HYMAP3_routing_struc(n)%rivlen(i) * &
                     HYMAP3_routing_struc(n)%rivwth(i,1))
             endif
          enddo

          call HYMAP3_set_fldstg(HYMAP3_routing_struc(n)%nz,&
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%fldhgt,&
               HYMAP3_routing_struc(n)%grarea,&
               HYMAP3_routing_struc(n)%rivlen,&
               HYMAP3_routing_struc(n)%rivwth(:,1),&
               HYMAP3_routing_struc(n)%rivstomax(:,1),&
               HYMAP3_routing_struc(n)%levhgt,&
               HYMAP3_routing_struc(n)%fldstomax(:,:,1),&
               HYMAP3_routing_struc(n)%fldgrd(:,:,1),&
               HYMAP3_routing_struc(n)%rivare(:,1),&
               HYMAP3_routing_struc(n)%fldonlystomax(:,:,1),&
               HYMAP3_routing_struc(n)%fldstoatlev(:,1))
       else
          do i=1,HYMAP3_routing_struc(n)%nseqall
             do m=1,LIS_rc%nensem(n)
                HYMAP3_routing_struc(n)%rivstomax(i,m) = &
                     HYMAP3_routing_struc(n)%rivlen(i) * &
                     HYMAP3_routing_struc(n)%rivwth(i,m) * &
                     HYMAP3_routing_struc(n)%rivhgt(i,m)
                !ag(23Feb2023) (22May2024)
                HYMAP3_routing_struc(n)%levstomax(i,m) = &
                     HYMAP3_routing_struc(n)%rivlen(i) * &
                     HYMAP3_routing_struc(n)%rivwth(i,m) * &
                     HYMAP3_routing_struc(n)%levhgt(i)
                HYMAP3_routing_struc(n)%rivelv(i) = &
                     HYMAP3_routing_struc(n)%elevtn(i) - &
                     HYMAP3_routing_struc(n)%rivhgt(i,m)
                if(HYMAP3_routing_struc(n)%rivwth(i,m)>0) then
                   HYMAP3_routing_struc(n)%rivare(i,m) = &
                        min(HYMAP3_routing_struc(n)%grarea(i), &
                        HYMAP3_routing_struc(n)%rivlen(i) * &
                        HYMAP3_routing_struc(n)%rivwth(i,m))
                endif
             enddo
          enddo

          do m=1,LIS_rc%nensem(n)
             call HYMAP3_set_fldstg(HYMAP3_routing_struc(n)%nz,&
                  HYMAP3_routing_struc(n)%nseqall,&
                  HYMAP3_routing_struc(n)%fldhgt,&
                  HYMAP3_routing_struc(n)%grarea,&
                  HYMAP3_routing_struc(n)%rivlen,&
                  HYMAP3_routing_struc(n)%rivwth(:,m),&
                  HYMAP3_routing_struc(n)%rivstomax(:,m),&
                  HYMAP3_routing_struc(n)%levhgt,&
                  HYMAP3_routing_struc(n)%fldstomax(:,:,m),&
                  HYMAP3_routing_struc(n)%fldgrd(:,:,m),&
                  HYMAP3_routing_struc(n)%rivare(:,m),&
                  HYMAP3_routing_struc(n)%fldstoatlev(:,m),&
                  HYMAP3_routing_struc(n)%fldonlystomax(:,:,m))
          enddo
       endif
    enddo

    !ag (4Feb2016) - read reservoir operation data
    do n=1, LIS_rc%nnest
       if(HYMAP3_routing_struc(n)%resopflag==1)then
          call HYMAP3_get_data_resop_alt(n, &
               HYMAP3_routing_struc(n)%resopdir,&
               HYMAP3_routing_struc(n)%resopheader,&
               HYMAP3_routing_struc(n)%nresop,&
               HYMAP3_routing_struc(n)%ntresop,&
               HYMAP3_routing_struc(n)%resoploc,&
               HYMAP3_routing_struc(n)%resoploc_dwn,&
               HYMAP3_routing_struc(n)%resoploc_glb,&
               HYMAP3_routing_struc(n)%resoploc_dwn_glb,&
               HYMAP3_routing_struc(n)%resopoutmin,&
               HYMAP3_routing_struc(n)%tresopalt,&
               HYMAP3_routing_struc(n)%resopalt,&
               HYMAP3_routing_struc(n)%resoptype)
          where(HYMAP3_routing_struc(n)%resopoutmin == &
               real(HYMAP3_routing_struc(n)%imis).or. &
               HYMAP3_routing_struc(n)%resopoutmin<0.) &
               HYMAP3_routing_struc(n)%resopoutmin=0.

          !ag(11Oct2024)
          allocate(elevtn_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(grarea_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivelv_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivlen_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivwth_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivstomax_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(fldhgt_glb(LIS_rc%glbnroutinggrid(n), &
               HYMAP3_routing_struc(n)%nz))
          allocate(fldstomax_glb(LIS_rc%glbnroutinggrid(n), &
               HYMAP3_routing_struc(n)%nz))

          call HYMAP3_gather_tiles(n,HYMAP3_routing_struc(n)%elevtn, &
               elevtn_glb)
          call HYMAP3_gather_tiles(n,HYMAP3_routing_struc(n)%grarea, &
               grarea_glb)
          call HYMAP3_gather_tiles(n,HYMAP3_routing_struc(n)%rivelv, &
               rivelv_glb)
          call HYMAP3_gather_tiles(n,HYMAP3_routing_struc(n)%rivlen, &
               rivlen_glb)
          call HYMAP3_gather_tiles(n, &
               HYMAP3_routing_struc(n)%rivwth(:,1),&
               rivwth_glb)
          call HYMAP3_gather_tiles(n, &
               HYMAP3_routing_struc(n)%rivstomax(:,1),rivstomax_glb)
          do iz=1,HYMAP3_routing_struc(n)%nz
             call HYMAP3_gather_tiles(n, &
                  HYMAP3_routing_struc(n)%fldhgt(:,iz),fldhgt_glb(:,iz))
             call HYMAP3_gather_tiles(n, &
                  HYMAP3_routing_struc(n)%fldstomax(:,iz,1), &
                  fldstomax_glb(:,iz))
          enddo

          do iresop=1,HYMAP3_routing_struc(n)%nresop
             icg=HYMAP3_routing_struc(n)%resoploc_glb(iresop)
             HYMAP3_routing_struc(n)%elevtn_resop(iresop)=elevtn_glb(icg)
             HYMAP3_routing_struc(n)%fldhgt_resop(iresop,:) = &
                  fldhgt_glb(icg,:)
             HYMAP3_routing_struc(n)%fldstomax_resop(iresop,:) = &
                  fldstomax_glb(icg,:)
             HYMAP3_routing_struc(n)%grarea_resop(iresop) = &
                  grarea_glb(icg)
             HYMAP3_routing_struc(n)%rivstomax_resop(iresop) = &
                  rivstomax_glb(icg)
             HYMAP3_routing_struc(n)%rivelv_resop(iresop) = &
                  rivelv_glb(icg)
             HYMAP3_routing_struc(n)%rivlen_resop(iresop)=rivlen_glb(icg)
             HYMAP3_routing_struc(n)%rivwth_resop(iresop)=rivwth_glb(icg)
          enddo
          deallocate(elevtn_glb)
          deallocate(fldhgt_glb)
          deallocate(grarea_glb)
          deallocate(rivelv_glb)
          deallocate(rivlen_glb)
          deallocate(rivwth_glb)
          deallocate(rivstomax_glb)
          deallocate(fldstomax_glb)

       !Yassin's scheme
       elseif(HYMAP3_routing_struc(n)%resopflag==2)then
          call HYMAP3_read_header_resop_yassin(n, &
               HYMAP3_routing_struc(n)%resopheader,&
               HYMAP3_routing_struc(n)%nresop,&
               HYMAP3_routing_struc(n)%resoploc,&
               HYMAP3_routing_struc(n)%resoploc_dwn,&
               HYMAP3_routing_struc(n)%resoploc_glb,&
               HYMAP3_routing_struc(n)%resoploc_dwn_glb,&
               HYMAP3_routing_struc(n)%ncloc_resop)
          do iresop=1,HYMAP3_routing_struc(n)%nresop
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'max_sto',&
                  HYMAP3_routing_struc(n)%maxsto_resop(iresop))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'init_dis',&
                  HYMAP3_routing_struc(n)%inidis_resop(iresop))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'init_sto',&
                  HYMAP3_routing_struc(n)%inisto_resop(iresop))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'down_dis',&
                  HYMAP3_routing_struc(n)%dwndis_resop(iresop))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile, &
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'dead_sto',&
                  HYMAP3_routing_struc(n)%deadis_resop(iresop))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile, &
                  12,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'mon_min_sto',&
                  HYMAP3_routing_struc(n)%minsto_mo_resop(iresop,:))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile, &
                  12,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'mon_nupper_sto',&
                  HYMAP3_routing_struc(n)%nupsto_mo_resop(iresop,:))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile, &
                  12,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'mon_upper_sto',&
                  HYMAP3_routing_struc(n)%uppsto_mo_resop(iresop,:))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  12,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'mon_min_dis',&
                  HYMAP3_routing_struc(n)%mindis_mo_resop(iresop,:))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile, &
                  12,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'mon_nupper_dis',&
                  HYMAP3_routing_struc(n)%nupdis_mo_resop(iresop,:))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  12,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'mon_upper_dis',&
                  HYMAP3_routing_struc(n)%uppdis_mo_resop(iresop,:))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'poly_a',&
                  HYMAP3_routing_struc(n)%reg1_resop(iresop))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'poly_b',&
                  HYMAP3_routing_struc(n)%reg2_resop(iresop))
             call read_netcdf_1d_real( &
                  HYMAP3_routing_struc(n)%resopncfile,&
                  1,1,HYMAP3_routing_struc(n)%ncloc_resop(iresop),1, &
                  'poly_c',&
                  HYMAP3_routing_struc(n)%reg3_resop(iresop))
          enddo
       endif
    enddo

    !ag(27Apr2020) - read urban flood data
    do n=1, LIS_rc%nnest
       if(HYMAP3_routing_struc(n)%flowtype==4)then
          write(LIS_logunit,*)'[INFO] Setting urban drainage parameters'
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP3 urban drainage parameter file:",rc=status)
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP3_routing_struc(n)%drfile,rc=status)
          call LIS_verify(status,&
               "HYMAP3 urban drainage parameter file: not defined")

          call HYMAP3_get_urban_parameters( &
               HYMAP3_routing_struc(n)%drfile,&
               HYMAP3_routing_struc(n)%drwth,&
               HYMAP3_routing_struc(n)%drhgt,&
               HYMAP3_routing_struc(n)%drden,&
               HYMAP3_routing_struc(n)%drvel,&
               HYMAP3_routing_struc(n)%drblk,&
               HYMAP3_routing_struc(n)%drrad,&
               HYMAP3_routing_struc(n)%drlgh,&
               HYMAP3_routing_struc(n)%drman,&
               HYMAP3_routing_struc(n)%drslp)

          call HYMAP3_gen_urban_drain_maps( &
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%drrad,&
               HYMAP3_routing_struc(n)%drlgh,&
               HYMAP3_routing_struc(n)%drden,&
               HYMAP3_routing_struc(n)%drwth,&
               HYMAP3_routing_struc(n)%drblk,&
               HYMAP3_routing_struc(n)%grarea,&
               HYMAP3_routing_struc(n)%next,&
               HYMAP3_routing_struc(n)%flowmap,&
               HYMAP3_routing_struc(n)%drstomax,&
               HYMAP3_routing_struc(n)%drtotwth,&
               HYMAP3_routing_struc(n)%drnoutlet,&
               HYMAP3_routing_struc(n)%drtotlgh)
       endif
    enddo

    !ag(8Aug2020) - read discharge data for direct insertion
    do n=1, LIS_rc%nnest
       if(HYMAP3_routing_struc(n)%insertflag==1)then
          call HYMAP3_get_discharge_data( &
               HYMAP3_routing_struc(n)%insertdir,&
               HYMAP3_routing_struc(n)%insertheader,&
               LIS_rc%gnc(n),LIS_rc%gnr(n), &
               HYMAP3_routing_struc(n)%sindex,&
               HYMAP3_routing_struc(n)%ninsert, &
               HYMAP3_routing_struc(n)%ntinsert,&
               HYMAP3_routing_struc(n)%insertloc,&
               HYMAP3_routing_struc(n)%tinsert,&
               HYMAP3_routing_struc(n)%insertdis)
       endif
    enddo

    !ag(30Mar2021) - read sea level data
    do n=1, LIS_rc%nnest
       if(HYMAP3_routing_struc(n)%sealevelflag==1)then
          call HYMAP3_get_sea_level_data(n, &
               HYMAP3_routing_struc(n)%sealeveldir,&
               HYMAP3_routing_struc(n)%sealevelheader,&
               HYMAP3_routing_struc(n)%outletlist,&
               HYMAP3_routing_struc(n)%nsealevel,&
               HYMAP3_routing_struc(n)%ntsealevel,&
               HYMAP3_routing_struc(n)%noutlet,&
               HYMAP3_routing_struc(n)%nseqall,&
               HYMAP3_routing_struc(n)%outletid,&
               HYMAP3_routing_struc(n)%tsealevel,&
               HYMAP3_routing_struc(n)%sealevel)
       endif
    enddo

    !ag(30Apr2021) - read water management header file
    do n=1, LIS_rc%nnest
       if(HYMAP3_routing_struc(n)%managflag==1)then
          !get number of water management locations
          call HYMAP3_read_header_size( &
               HYMAP3_routing_struc(n)%managheader, &
               HYMAP3_routing_struc(n)%nmanag)
          HYMAP3_routing_struc(n)%nmanagcoef=12
          allocate(HYMAP3_routing_struc(n)%managqmax( &
               HYMAP3_routing_struc(n)%nmanag))
          allocate(HYMAP3_routing_struc(n)%managact( &
               HYMAP3_routing_struc(n)%nmanag))
          allocate(HYMAP3_routing_struc(n)%managloc( &
               HYMAP3_routing_struc(n)%nmanag,2))
          allocate(HYMAP3_routing_struc(n)%managtype( &
               HYMAP3_routing_struc(n)%nmanag))
          allocate(HYMAP3_routing_struc(n)%managcoef( &
               HYMAP3_routing_struc(n)%nmanag, &
               HYMAP3_routing_struc(n)%nmanagcoef))

          call HYMAP3_get_management_rules( &
               HYMAP3_routing_struc(n)%managheader,&
               LIS_rc%gnc(n),LIS_rc%gnr(n), &
               HYMAP3_routing_struc(n)%sindex, &
               HYMAP3_routing_struc(n)%nmanag, &
               HYMAP3_routing_struc(n)%nmanagcoef,&
               HYMAP3_routing_struc(n)%managloc,&
               HYMAP3_routing_struc(n)%managqmax,&
               HYMAP3_routing_struc(n)%managact,&
               HYMAP3_routing_struc(n)%managtype,&
               HYMAP3_routing_struc(n)%managcoef)
       endif
    enddo

    !ag(7Sep2022) - read bifurcation pathway file
    do n=1, LIS_rc%nnest
       if(HYMAP3_routing_struc(n)%bifflag==1)then
          !ag(27Jul2025)
          !get number of bifurcations and elevations
          allocate(HYMAP3_routing_struc(n)%bifloc( &
               HYMAP3_routing_struc(n)%nbif,2))
          allocate(HYMAP3_routing_struc(n)%bifelv( &
               HYMAP3_routing_struc(n)%nbif))
          allocate(HYMAP3_routing_struc(n)%bifdelv( &
               HYMAP3_routing_struc(n)%nbifelv))
          allocate(HYMAP3_routing_struc(n)%bifman( &
               HYMAP3_routing_struc(n)%nbifelv))
          allocate(HYMAP3_routing_struc(n)%bifwth( &
               HYMAP3_routing_struc(n)%nbif, &
               HYMAP3_routing_struc(n)%nbifelv))
          allocate(HYMAP3_routing_struc(n)%biflen( &
               HYMAP3_routing_struc(n)%nbif))
          allocate(HYMAP3_routing_struc(n)%bifout( &
               HYMAP3_routing_struc(n)%nbif))
          allocate(HYMAP3_routing_struc(n)%bifout_pre( &
               HYMAP3_routing_struc(n)%nbif))
          allocate(HYMAP3_routing_struc(n)%bifsto( &
               HYMAP3_routing_struc(n)%nbif, &
               HYMAP3_routing_struc(n)%nbifelv))
          HYMAP3_routing_struc(n)%bifout_pre=0.

          call HYMAP3_get_bifurcation_pathways( &
               HYMAP3_routing_struc(n)%biffile,&
               LIS_rc%gnc(n),LIS_rc%gnr(n), &
               HYMAP3_routing_struc(n)%sindex,&
               HYMAP3_routing_struc(n)%nbif, &
               HYMAP3_routing_struc(n)%nbifelv,&
               HYMAP3_routing_struc(n)%bifloc,&
               HYMAP3_routing_struc(n)%bifdelv,&
               HYMAP3_routing_struc(n)%bifman,&
               HYMAP3_routing_struc(n)%bifelv,&
               HYMAP3_routing_struc(n)%biflen,&
               HYMAP3_routing_struc(n)%bifwth)

          allocate(elevtn_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(grarea_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivelv_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivlen_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivwth_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(rivstomax_glb(LIS_rc%glbnroutinggrid(n)))
          allocate(fldhgt_glb(LIS_rc%glbnroutinggrid(n), &
               HYMAP3_routing_struc(n)%nz))
          allocate(fldstomax_glb(LIS_rc%glbnroutinggrid(n), &
               HYMAP3_routing_struc(n)%nz))

          call HYMAP3_gather_tiles(n, &
               HYMAP3_routing_struc(n)%elevtn,elevtn_glb)
          call HYMAP3_gather_tiles(n,HYMAP3_routing_struc(n)%grarea, &
               grarea_glb)
          call HYMAP3_gather_tiles(n,HYMAP3_routing_struc(n)%rivelv, &
               rivelv_glb)
          call HYMAP3_gather_tiles(n,HYMAP3_routing_struc(n)%rivlen, &
               rivlen_glb)
          call HYMAP3_gather_tiles(n, &
               HYMAP3_routing_struc(n)%rivwth(:,1),rivwth_glb)
          call HYMAP3_gather_tiles(n, &
               HYMAP3_routing_struc(n)%rivstomax(:,1),rivstomax_glb)

          do iz=1,HYMAP3_routing_struc(n)%nz
             call HYMAP3_gather_tiles(n, &
                  HYMAP3_routing_struc(n)%fldhgt(:,iz),fldhgt_glb(:,iz))
             call HYMAP3_gather_tiles(n, &
                  HYMAP3_routing_struc(n)%fldstomax(:,iz,1), &
                  fldstomax_glb(:,iz))
          enddo

          do ibif=1,HYMAP3_routing_struc(n)%nbif
             icg=HYMAP3_routing_struc(n)%bifloc(ibif,1)
             do ielv=1,HYMAP3_routing_struc(n)%nbifelv
                bifelv1 = &
                     HYMAP3_routing_struc(n)%bifelv(ibif) + &
                     HYMAP3_routing_struc(n)%bifdelv(ielv)
                call HYMAP3_get_volume_profile( &
                     HYMAP3_routing_struc(n)%nz,dble(elevtn_glb(icg)), &
                     dble(fldhgt_glb(icg,:)),dble(fldstomax_glb(icg,:)),&
                     dble(grarea_glb(icg)),dble(rivstomax_glb(icg)), &
                     dble(rivelv_glb(icg)),dble(rivlen_glb(icg)), &
                     dble(rivwth_glb(icg)),bifelv1,bifsto1)
                HYMAP3_routing_struc(n)%bifsto(ibif,ielv)=bifsto1
             enddo
          enddo
          deallocate(elevtn_glb)
          deallocate(fldhgt_glb)
          deallocate(grarea_glb)
          deallocate(rivelv_glb)
          deallocate(rivlen_glb)
          deallocate(rivwth_glb)
          deallocate(rivstomax_glb)
          deallocate(fldstomax_glb)
       endif
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model start mode:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%startmode,rc=status)
       call LIS_verify(status,&
            "HYMAP3 routing model start mode: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model restart interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP3 routing model restart interval: not defined")
       call LIS_parseTimeString(time,HYMAP3_routing_struc(n)%rstInterval)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP3 routing model restart file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP3_routing_struc(n)%rstfile,rc=status)
       call LIS_verify(status,&
            "HYMAP3 routing model restart file: not defined")
    enddo

    if(LIS_rc%lsm.eq."none") then
       call initrunoffdata(trim(LIS_rc%runoffdatasource)//char(0))
    endif

    do n=1, LIS_rc%nnest
       call ESMF_ArraySpecSet(realarrspec,rank=1, &
            typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status)

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

       call ESMF_stateAdd(LIS_runoff_state(n),(/sf_runoff_field/), &
            rc=status)
       call LIS_verify(status, 'ESMF_StateAdd failed for surface runoff')

       call ESMF_stateAdd(LIS_runoff_state(n),(/baseflow_field/), &
            rc=status)
       call LIS_verify(status, 'ESMF_StateAdd failed for base flow')

       !hkb (4Mar2016)
       !create LSM interface objects to store evapotranspiration and
       !potential evaporation only if source is readin
       if ( HYMAP3_routing_struc(n)%evapflag .ne. 0 ) then
          evapotranspiration_field = &
               ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=LIS_vecTile(n), name="Total Evapotranspiration", &
               rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')

          call ESMF_FieldGet(evapotranspiration_field,localDE=0,&
               farrayPtr=evapotranspiration,&
               rc=status)
          call LIS_verify(status)
          evapotranspiration = 0.0

          call ESMF_stateAdd(LIS_runoff_state(n), &
               (/evapotranspiration_field/),rc=status)
          call LIS_verify(status, &
               'ESMF_StateAdd failed for Total Evapotranspiration')
       endif

       HYMAP3_routing_struc(n)%mo = -1

       !ag (12Sep2019)
       call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling",&
            0, rc=status)
       call LIS_verify(status)

       if (HYMAP3_routing_struc(n)%enable2waycpl==1) then
          ! River Storage
          rivsto_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=LIS_vecTile(n), name="River Storage",rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')

          call ESMF_FieldGet(rivsto_field,localDE=0,farrayPtr=rivstotmp,&
               rc=status)
          call LIS_verify(status)
          rivstotmp = 0.0

          call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling",&
               HYMAP3_routing_struc(n)%enable2waycpl, rc=status)
          call LIS_verify(status)

          call ESMF_stateAdd(LIS_runoff_state(n),(/rivsto_field/), &
               rc=status)
          call LIS_verify(status, &
               'ESMF_StateAdd failed for River Storage')

          ! Flood Storage
          fldsto_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=LIS_vecTile(n), name="Flood Storage",rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')

          call ESMF_FieldGet(fldsto_field,localDE=0,farrayPtr=fldstotmp,&
               rc=status)
          call LIS_verify(status)
          fldstotmp = 0.0

          call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling",&
               HYMAP3_routing_struc(n)%enable2waycpl, rc=status)
          call LIS_verify(status)

          call ESMF_stateAdd(LIS_runoff_state(n),(/fldsto_field/), &
               rc=status)
          call LIS_verify(status, &
               'ESMF_StateAdd failed for Flood Storage')

          ! Flooded fraction
          fldfrc_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=LIS_vecTile(n), name="Flooded Fraction",rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')

          call ESMF_FieldGet(fldfrc_field,localDE=0,farrayPtr=fldfrctmp,&
               rc=status)
          call LIS_verify(status)
          fldfrctmp = 0.0

          call ESMF_AttributeSet(LIS_runoff_state(n),"2 way coupling", &
               HYMAP3_routing_struc(n)%enable2waycpl, rc=status)
          call LIS_verify(status)

          call ESMF_stateAdd(LIS_runoff_state(n),(/fldfrc_field/), &
               rc=status)
          call LIS_verify(status, &
               'ESMF_StateAdd failed for Flooded Fraction')
       endif

    enddo

    do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(LIS_runoff_state(n), &
            "Routing model evaporation option",&
            HYMAP3_routing_struc(n)%evapflag, rc=status)
       call LIS_verify(status)

       call LIS_registerAlarm("HYMAP3 router model alarm",&
            HYMAP3_routing_struc(n)%dt,HYMAP3_routing_struc(n)%dt)

       call LIS_registerAlarm("HYMAP3 router output alarm",&
            HYMAP3_routing_struc(n)%dt, &
            HYMAP3_routing_struc(n)%outInterval)

       call LIS_registerAlarm("HYMAP3 router restart alarm",&
            HYMAP3_routing_struc(n)%dt, &
            HYMAP3_routing_struc(n)%rstInterval)
    enddo


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
            ESMF_GridCreate(name="HYMAP3 Patch Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = patchDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LIS_verify(status, &
            'ESMF_GridCreate failed in HYMAP3_routing_init')

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
            ESMF_GridCreate(name="HYMAP3 Tile Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = gridDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LIS_verify(status, &
            'ESMF_GridCreate failed in HYMAP3_routing_init')
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
          Routing_DAvalid = Routing_DAvalid.and. &
               LIS_rc%Routing_DAinst_valid(i)
       enddo

       allocate(LIS_Routing_State(LIS_rc%nnest, LIS_rc%nperts))
       allocate(LIS_Routing_Incr_State(LIS_rc%nnest, LIS_rc%nperts))

       do n=1,LIS_rc%nnest
          do k=1,LIS_rc%nperts
             write(LIS_logunit,*) &
                  '[INFO] Opening constraints for prognostic state variables ',&
                  LIS_rc%progattribFile(k)
             ftn = LIS_getNextUnitNumber()
             open(ftn, file = LIS_rc%progattribFile(k),status='old')
             read(ftn,*)
             read(ftn,*) LIS_rc%nstvars(k)
             read(ftn,*)

             allocate(vname(LIS_rc%nstvars(k)))
             allocate(stmin(LIS_rc%nstvars(k)))
             allocate(stmax(LIS_rc%nstvars(k)))

             call ESMF_ArraySpecSet(arrspec1,rank=1, &
                  typekind=ESMF_TYPEKIND_R4,&
                  rc=status)
             call LIS_verify(status, &
                  "ESMF_ArraySpecSet failed in LIS_routing_init")

             write(unit=temp,fmt='(i2.2)') n
             read(unit=temp,fmt='(2a1)') nestid

             write(unit=temp,fmt='(i3.3)') k
             read(unit=temp,fmt='(3a1)') caseid

             LIS_Routing_State(n,k) = &
                  ESMF_StateCreate(name="Routing State"//&
                  nestid(1)//nestid(2)&
                  //'_'//caseid(1)//caseid(2)//caseid(3), rc=status)
             call LIS_verify(status, &
                  "ESMF_StateCreate failed in LIS_routing_init")

             LIS_Routing_Incr_State(n,k) = &
                  ESMF_StateCreate(name="Routing Incr State"//&
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

                call ESMF_AttributeSet(varField,"Max Value",stmax(i), &
                     rc=status)
                call LIS_verify(status,&
                     "ESMF_AttribteSet failed in LIS_routing_init")

                call ESMF_AttributeSet(varField,"Min Value",stmin(i), &
                     rc=status)
                call LIS_verify(status,&
                     "ESMF_AttributeSet failed in LIS_routing_init")

                call ESMF_AttributeSet(VarIncrField,"Max Value", &
                     stmax(i),rc=status)
                call LIS_verify(status,&
                     "ESMF_AttributeSet failed in LIS_routing_init")

                call ESMF_AttributeSet(VarIncrField,"Min Value", &
                     stmin(i),rc=status)
                call LIS_verify(status,&
                     "ESMF_AttributeSet failed in LIS_routing_init")

                call ESMF_StateAdd(LIS_Routing_State(n,k),(/varField/), &
                     rc=status)
                call LIS_verify(status,&
                     "ESMF_StateAdd failed in LIS_routing_init")

                call ESMF_StateAdd(LIS_Routing_Incr_State(n,k), &
                     (/VarIncrField/), rc=status)
                call LIS_verify(status,&
                     "ESMF_StateAdd failed in LIS_routing_init")
                !--------------------------------------------------------
                ! Initially set the fresh increments available status to
                ! false.
                !--------------------------------------------------------
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
                allocate(routing_pert%ccorr(LIS_rc%nstvars(k), &
                     LIS_rc%nstvars(k)))

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
                        arrayspec=arrspec2, &
                        name=trim(routing_pert%vname(i)),&
                        rc=status)

                   call ESMF_StateAdd(LIS_Routing_Pert_State(n,k), &
                        (/pertField/),&
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
                         exit
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

                      call ESMF_AttributeSet(pertField, &
                           "Standard Deviation",&
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
             call perturbsetup(trim(LIS_rc%perturb_state(i))//char(0), &
                  4, i, &
                  LIS_Routing_State(:,i), LIS_Routing_Pert_State(:,i))
          endif
       enddo
    end if

  end subroutine HYMAP3_routingInit
  !=============================================
  subroutine HYMAP3_vector_read_dims(n)

    !USES:
    use LIS_coreMod
    use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none

! !ARGUMENTS:
    integer, intent(in)         :: n

    integer :: ios,nid,ncId
    logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=trim(HYMAP3_routing_struc(n)%vecfile), &
         exist=file_exists)
    if(file_exists) then

       write(LIS_logunit,*)'[INFO] Read HYMAP3 vector size ',&
            'from ',trim(HYMAP3_routing_struc(n)%vecfile)

       ios = nf90_open(path=trim(HYMAP3_routing_struc(n)%vecfile),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios, &
            'Error in nf90_open in HYMAP3_vector_read_dims')

       ios = nf90_inq_dimid(nid,"elv_layer",ncId)
       call LIS_verify(ios, &
            'Error in nf90_inq_dimid in HYMAP3_vector_read_dims')

       ios = nf90_inquire_dimension(nid,ncId, &
            len=HYMAP3_routing_struc(n)%nz)

       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in HYMAP3_vector_read_dims')

       ios = nf90_inq_dimid(nid,"domain_size",ncId)

       call LIS_verify(ios, &
            'Error in nf90_inq_dimid in HYMAP3_vector_read_dims')

       ios = nf90_inquire_dimension(nid,ncId, &
            len=HYMAP3_routing_struc(n)%nseqall)

       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in HYMAP3_vector_read_dims')

       ios = nf90_close(nid)

       call LIS_verify(ios, &
            'Error in nf90_close in HYMAP3_vector_read_dims')
    else
       write(LIS_logunit,*) '[ERR] file: ', &
            trim(HYMAP3_routing_struc(n)%vecfile), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine HYMAP3_vector_read_dims

  !=============================================
  subroutine HYMAP3_vector_read_param(ctitle, n, z, array)

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
    real,         intent(inout) :: &
         array(HYMAP3_routing_struc(n)%nseqall,z)

    integer                     :: ftn
    logical                     :: file_exists
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=HYMAP3_routing_struc(n)%vecfile,exist=file_exists)
    if(file_exists) then

       call LIS_verify(nf90_open(path=HYMAP3_routing_struc(n)%vecfile,&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in HYMAP3_vector_read_param in HYMAP3_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in HYMAP3_vector_read_param in HYMAP3_routingMod')

       call LIS_verify(nf90_get_var(ftn,varid, array), &
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in HYMAP3_vector_read_param in HYMAP3_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP3_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '// &
            trim(HYMAP3_routing_struc(n)%vecfile)
       write(LIS_logunit,*) &
            '[ERR] failed in HYMAP3_vector_read_param in HYMAP3_routingMod'
       call LIS_endrun()
    endif

#endif
  end subroutine HYMAP3_vector_read_param

  !======================================================================
  subroutine read_netcdf_1d_real(infile, ntiles, jt, it, nt, yvar, var)

    !USES:
    use LIS_coreMod
    use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none

    character(*), intent(in)  :: infile,yvar
    integer,      intent(in)  :: ntiles,jt,it,nt
    real*8,         intent(out) :: var(ntiles)

    integer                   :: varid
    integer                     :: ftn
    logical                     :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    inquire(file=trim(infile),exist=file_exists)
    if(file_exists) then
       call LIS_verify(nf90_open(path=trim(infile),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_netcdf_1d_real in HYMAP3_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(yvar),varid), &
            'nf90_inq_varid failed for '//trim(yvar)//&
            ' in read_netcdf_1d_real in HYMAP3_routingMod')

       call LIS_verify(nf90_get_var(ftn,varid, var,&
            start=(/it,jt/),&
            count=(/nt,ntiles/)),&
            'nf90_get_var failed for '//trim(yvar)//&
            ' in read_netcdf_1d_real in HYMAP3_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP3_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(infile)
       write(LIS_logunit,*) &
            '[ERR] failed in read_netcdf_1d_real in HYMAP3_routingMod'
       call LIS_endrun()
    endif
#endif

  end subroutine read_netcdf_1d_real
  !===============================================================================
  subroutine HYMAP3_read_header_resop_yassin(n, yheader, inst, &
       local_index, down_local_index, glb_index, down_glb_index, &
       nc_index)

    use LIS_logMod

    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(in)    :: n,inst
    integer,      intent(inout) :: local_index(inst), &
         down_local_index(inst)
    integer,      intent(inout) :: glb_index(inst),down_glb_index(inst)
    integer,      intent(inout) :: nc_index(inst)

    logical                     :: file_exists
    integer                     :: ist,ix,iy,ix_down,iy_down
    integer :: ftn

    inquire(file=yheader,exist=file_exists)

    if(file_exists) then
      !get name of station files
       write(LIS_logunit,*) &
            '[INFO] [read_header] get stations info: name and coordinates'
       write(LIS_logunit,*)'[INFO] [read_header] ',yheader,inst
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(yheader), status='old')
       read(ftn,*)
       do ist=1,inst
          read(ftn,*,end=10)nc_index(ist),ix,iy
          call HYMAP3_map_gxy2l_index(n,ix,iy,local_index(ist))
          glb_index(ist) = HYMAP3_routing_struc(n)%sindex(ix,iy)
          down_glb_index(ist) = &
               HYMAP3_routing_struc(n)%next_glb(glb_index(ist))
          ix_down=HYMAP3_routing_struc(n)%nextx(ix,iy)
          iy_down=HYMAP3_routing_struc(n)%nexty(ix,iy)
          call HYMAP3_map_gxy2l_index(n,ix_down,iy_down, &
               down_local_index(ist))
          write(LIS_logunit,'(a,10i8)')'[INFO] [read_header] ',ist, &
               nc_index(ist),ix,iy,ix_down,iy_down,local_index(ist), &
               down_local_index(ist),glb_index(ist),down_glb_index(ist)
       enddo
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) '[ERR] header file '//trim(yheader)
       write(LIS_logunit,*) '[ERR] failed in read_header in HYMAP3_routingMod'
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) '[ERR] header file '//trim(yheader)
    write(LIS_logunit,*) &
         '[ERR ]failed in read_header in HYMAP3_routingMod'
    call LIS_endrun()

  end subroutine HYMAP3_read_header_resop_yassin
  !======================================================================
  subroutine HYMAP3_read_param_real(ctitle, n, z, array)

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
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)

    if(file_exists) then

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_real in HYMAP3_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP3_routingMod')

       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
            LIS_nss_halo_ind(n,LIS_localPet+1)+1),z/)),&
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP3_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP3_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '// &
            trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) &
            '[ERR] failed in read_param_real in HYMAP3_routingMod'
       call LIS_endrun()
    endif

#endif
  end subroutine HYMAP3_read_param_real

  !=============================================
  subroutine HYMAP3_read_param_real_2d(ctitle, n, array)

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
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)

    if(file_exists) then

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_real in HYMAP3_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP3_routingMod')

       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1)/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
            LIS_nss_halo_ind(n,LIS_localPet+1)+1)/)),&
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_real in HYMAP3_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP3_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '// &
            trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) &
            '[ERR] failed in read_param_real in HYMAP3_routingMod'
       call LIS_endrun()
    endif

#endif
  end subroutine HYMAP3_read_param_real_2d

  subroutine HYMAP3_read_param_int(ctitle, z, n, array)

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
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)

    if(file_exists) then

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_int in HYMAP3_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP3_routingMod')

       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1),z/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP3_routingMod')

       call LIS_verify(nf90_close(ftn), &
            'nf90_close failed in HYMAP3_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '// &
            trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) &
            '[ERR] failed in read_param_int in HYMAP3_routingMod'
       call LIS_endrun()
    endif

#endif
  end subroutine HYMAP3_read_param_int

  subroutine HYMAP3_read_param_int_2d(ctitle, n, array)

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
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)

    if(file_exists) then

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_int in HYMAP3_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP3_routingMod')

       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1)/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
            LIS_nss_halo_ind(n,LIS_localPet+1)+1)/)),&
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP3_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP3_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '// &
            trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) &
            '[ERR] failed in read_param_int in HYMAP3_routingMod'
       call LIS_endrun()
    endif

#endif
  end subroutine HYMAP3_read_param_int_2d

  subroutine HYMAP3_read_param_int_2d_global(ctitle, n, array)

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
    integer                     :: varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)

    if(file_exists) then

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in read_param_int in HYMAP3_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP3_routingMod')

       call LIS_verify(nf90_get_var(ftn,varid, array), &
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int in HYMAP3_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP3_routingMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '// &
            trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) &
            '[ERR] failed in read_param_int in HYMAP3_routingMod'
       call LIS_endrun()
    endif

#endif
  end subroutine HYMAP3_read_param_int_2d_global

  !=============================================
  !=============================================
  subroutine HYMAP3_get_data_resop_alt(n, resopdir, resopheader, &
       nresop, &
       ntresop, resoploc, resoploc_dwn, resoploc_glb, resoploc_dwn_glb, &
       resopoutmin, tresop, resop, resoptype)

    use LIS_constantsMod, only: LIS_CONST_PATH_LEN

    implicit none

    character(*), intent(in)  :: resopdir,resopheader
    integer,      intent(in)  :: n,nresop,ntresop
    integer,      intent(out) :: resoploc(nresop), &
         resoploc_dwn(nresop),resoploc_glb(nresop), &
         resoploc_dwn_glb(nresop),resoptype(nresop)
    real*8,       intent(out) :: tresop(nresop,ntresop)
    real,         intent(out) :: resop(nresop,ntresop), &
         resopoutmin(nresop)

    integer                   :: res
    character(50)             :: resopname(nresop)
    character(LIS_CONST_PATH_LEN)            :: yfile

    call HYMAP3_read_header_resop(n,trim(resopheader),nresop,resopname, &
         resoploc,resoploc_dwn,resoploc_glb,resoploc_dwn_glb, &
         resopoutmin,resoptype)

    do res=1,nresop
       yfile=trim(resopdir)//trim(resopname(res))//'.txt'
       call HYMAP3_read_time_series(ntresop,trim(yfile),tresop(res,:), &
            resop(res,:))
    enddo
  end subroutine HYMAP3_get_data_resop_alt
  !=============================================
  !=============================================
  subroutine HYMAP3_get_sea_level_data(n, sealeveldir, sealevelheader, &
       outletlist, nsealevel, ntsealevel, noutlet, nseqall, outletid,  &
       tsealevel, sealevel)

    use LIS_constantsMod, only: LIS_CONST_PATH_LEN

    implicit none

    character(*), intent(in)  :: sealeveldir,sealevelheader,outletlist
    integer,      intent(in)  :: n
    integer,      intent(in)  :: nsealevel,ntsealevel,noutlet,nseqall
    integer,      intent(out) :: outletid(nseqall)
    real*8,       intent(out) :: tsealevel(nsealevel,ntsealevel)
    real,         intent(out) :: sealevel(nsealevel,ntsealevel)

    integer                   :: ins
    integer                   :: outlet_index(noutlet), &
         sealevelloc(noutlet)
    character(50)             :: sealevelname(nsealevel)
    character(LIS_CONST_PATH_LEN) :: yfile

    !read outlet list file
    call HYMAP3_read_outlet_list(n,trim(outletlist),noutlet, &
         sealevelloc,outlet_index)
    do ins=1,noutlet
       if(outlet_index(ins)>0)then
          outletid(outlet_index(ins))=sealevelloc(ins)
       endif
    enddo
    !read sea level header file
    call HYMAP3_read_sea_level_header(trim(sealevelheader),nsealevel, &
         sealevelname)
    do ins=1,nsealevel
       yfile=trim(sealeveldir)//trim(sealevelname(ins))//'.txt'
       call HYMAP3_read_time_series(ntsealevel,trim(yfile), &
            tsealevel(ins,:),sealevel(ins,:))
    enddo
  end subroutine HYMAP3_get_sea_level_data
  !=============================================
  !=============================================
  subroutine HYMAP3_read_header_size(yheader, isize)

    use LIS_logMod

    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(out)   :: isize

    logical                     :: file_exists
    integer :: ftn

    inquire(file=yheader,exist=file_exists)
    if(file_exists) then
       write(LIS_logunit,*)'[INFO] read header size '//trim(yheader)
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(yheader), status='old')
       read(ftn,*,end=10)isize
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) '[ERR] header '//trim(yheader)
       write(LIS_logunit,*) &
            '[ERR] failed in opening file in HYMAP3_routingMod'
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) '[ERR] check header size in file '// &
         trim(yheader)
    call LIS_endrun()

  end subroutine HYMAP3_read_header_size
  !=============================================
  !=============================================
  !ag(27Jul2025)
  subroutine HYMAP3_read_header_size1(yheader, isize, isize1)

    use LIS_logMod

    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(out)   :: isize,isize1

    logical                     :: file_exists
    integer :: ftn

    inquire(file=yheader,exist=file_exists)
    if(file_exists) then
       write(LIS_logunit,*)'[INFO] read header size '//trim(yheader)
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(yheader), status='old')
       read(ftn,*,end=10)isize,isize1
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) '[ERR] header '//trim(yheader)
       write(LIS_logunit,*) &
            '[ERR] failed in opening file in HYMAP3_routingMod'
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) '[ERR] check header size in file '// &
         trim(yheader)
    call LIS_endrun()

  end subroutine HYMAP3_read_header_size1
  !=============================================
  !=============================================
  subroutine HYMAP3_read_outlet_list(n, yheader, inst, id, local_index)

    use LIS_logMod

    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(in)    :: n,inst
    integer,      intent(inout) :: id(inst),local_index(inst)

    logical                     :: file_exists
    integer                     :: ist,ix,iy
    integer :: ftn

    inquire(file=yheader,exist=file_exists)
    if(file_exists) then
      !get name of station files
       write(LIS_logunit,*) &
            '[INFO] [read_outlet_list] get stations info: name and coordinates'
       write(LIS_logunit,*)'[INFO] [read_outlet_list] ',yheader,inst
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(yheader), status='old')
       !first row is the file length
       read(ftn,*)
       do ist=1,inst
          read(ftn,*,end=10)id(ist),ix,iy
          call HYMAP3_map_gxy2l_index(n,ix,iy,local_index(ist))
       enddo
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) '[ERR] outlet list file '//trim(yheader)
       write(LIS_logunit,*) &
            '[ERR] failed in read_outlet_list in HYMAP3_routingMod'
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) '[ERR] outlet list file '//trim(yheader)
    write(LIS_logunit,*) &
         '[ERR] failed in read_outlet_list in HYMAP3_routingMod'
    call LIS_endrun()

  end subroutine HYMAP3_read_outlet_list
  !=============================================
  !=============================================
  subroutine HYMAP3_read_sea_level_header(yheader, inst, yqname)

    use LIS_logMod

    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(in)    :: inst
    character(*), intent(inout) :: yqname(inst)

    logical                     :: file_exists
    integer                     :: ist
    integer :: ftn

    inquire(file=yheader,exist=file_exists)
    if(file_exists) then
       !get name of station files
       write(LIS_logunit,*) &
            '[INFO] [read_sea_level_header] get stations info: name and coordinates'
       write(LIS_logunit,*)'[INFO] [read_sea_level_header] ',yheader,inst
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(yheader), status='old')
       !first row is the file length
       read(ftn,*)
       do ist=1,inst
          read(ftn,*,end=10)yqname(ist)
          write(LIS_logunit,*)'[INFO] [read_sea_level_header] ', &
               ist,trim(yqname(ist))
       enddo
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) '[ERR] header file '//trim(yheader)
       write(LIS_logunit,*) &
            '[ERR] failed in read_sea_level_header in HYMAP3_routingMod'
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) '[ERR] header file '//trim(yheader)
    write(LIS_logunit,*) &
         '[ERR] failed in read_sea_level_header in HYMAP3_routingMod'
    call LIS_endrun()

  end subroutine HYMAP3_read_sea_level_header
  !=============================================
  !=============================================
  subroutine HYMAP3_get_discharge_data(insertdir, insertheader, nx, ny, &
       sindex, ninsert, ntinsert, insertloc, tinsert, insertdis)

    use LIS_constantsMod, only: LIS_CONST_PATH_LEN

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
    character(50)             :: insertname(ninsert)
    character(LIS_CONST_PATH_LEN)            :: yfile

    call HYMAP3_read_header(trim(insertheader),ninsert,insertname, &
         xinsert,yinsert)
    do ins=1,ninsert
       insertloc(ins)=sindex(xinsert(ins),yinsert(ins))
    enddo
    do ins=1,ninsert
       yfile=trim(insertdir)//trim(insertname(ins))//'.txt'
       call HYMAP3_read_time_series(ntinsert,trim(yfile), &
            tinsert(ins,:),insertdis(ins,:))
    enddo
  end subroutine HYMAP3_get_discharge_data
  !=============================================
  !=============================================
  subroutine HYMAP3_read_header_resop(n, yheader, inst, yqname,  &
       local_index, down_local_index, glb_index, down_glb_index, &
       outmin, resoptype)

    use LIS_logMod

    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(in)    :: n,inst
    character(*), intent(inout) :: yqname(inst)
    integer,      intent(inout) :: local_index(inst),&
         down_local_index(inst)
    integer,      intent(inout) :: glb_index(inst),down_glb_index(inst)
    integer,      intent(inout) :: resoptype(inst)
    real,         intent(inout) :: outmin(inst)

    logical                     :: file_exists
    integer                     :: ist,ix,iy,ix_down,iy_down
    integer :: ftn

    inquire(file=yheader,exist=file_exists)
    if(file_exists) then
       ftn = LIS_getNextUnitNumber()
       !get name of station files
       write(LIS_logunit,*) &
            '[INFO] [read_header] get stations info: name and coordinates'
       write(LIS_logunit,*)'[INFO] [read_header] ',yheader,inst
       open(ftn,file=trim(yheader), status='old')
       !ag(27Jul2025)
       read(ftn,*)
       do ist=1,inst
          read(ftn,*,end=10)yqname(ist),ix,iy,outmin(ist),resoptype(ist)
          call HYMAP3_map_gxy2l_index(n,ix,iy,local_index(ist))
          glb_index(ist) = HYMAP3_routing_struc(n)%sindex(ix,iy)
          down_glb_index(ist) = &
               HYMAP3_routing_struc(n)%next_glb(glb_index(ist))
          ix_down=HYMAP3_routing_struc(n)%nextx(ix,iy)
          iy_down=HYMAP3_routing_struc(n)%nexty(ix,iy)
          call HYMAP3_map_gxy2l_index(n,ix_down,iy_down, &
               down_local_index(ist))
          write(LIS_logunit,'(a,i5,a,4i5,f10.2,5i8)') &
               '[INFO] [read_header] ',ist,trim(yqname(ist)),ix,iy, &
               ix_down,iy_down,outmin(ist),resoptype(ist), &
               local_index(ist),down_local_index(ist),glb_index(ist), &
               down_glb_index(ist)
       enddo
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) '[ERR] header file '//trim(yheader)
       write(LIS_logunit,*) &
            '[ERR] failed in read_header in HYMAP3_routingMod'
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) &
         '[ERR] HYMAP3_read_header_resop: canno read header file ' &
         //trim(yheader)
    call LIS_endrun()

  end subroutine HYMAP3_read_header_resop
  !=============================================
  !=============================================
  subroutine HYMAP3_read_header(yheader, inst, yqname, ix, iy)

    use LIS_logMod

    implicit none

    character(*), intent(in)    :: yheader
    integer,      intent(in)    :: inst
    character(*), intent(inout) :: yqname(inst)
    integer,      intent(inout) :: ix(inst),iy(inst)

    logical                     :: file_exists
    integer                     :: ist
    integer :: ftn

    inquire(file=yheader,exist=file_exists)
    if(file_exists) then
       ftn = LIS_getNextUnitNumber()
       !get name of station files
       write(LIS_logunit,*) &
            '[INFO] [read_header] get stations info: name and coordinates'
       write(LIS_logunit,*)'[INFO] [read_header] ',yheader,inst
       open(ftn,file=trim(yheader), status='old')
       do ist=1,inst
          read(ftn,*,end=10)yqname(ist),ix(ist),iy(ist)
          write(LIS_logunit,*) &
               '[INFO] [read_header] ',ist,trim(yqname(ist)),ix(ist), &
               iy(ist)
       enddo
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) &
            '[ERR] HYMAP3_read_header:  Cannot read header file ' &
            //trim(yheader)
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) &
         '[ERR] HYMAP3_read_header:  Cannot read header file ' &
         //trim(yheader)
    call LIS_endrun()

  end subroutine HYMAP3_read_header
  !=============================================
  !=============================================
  subroutine HYMAP3_read_time_series(itmax, yfile, ztalt, zhalt)

    use LIS_logMod

    implicit none

    integer,      intent(in)    :: itmax
    character(*), intent(in)    :: yfile
    real*8,       intent(inout) :: ztalt(itmax)
    real,         intent(inout) :: zhalt(itmax)

    logical                     :: file_exists
    integer :: it
    integer :: ftn

    inquire(file=yfile,exist=file_exists)
    if(file_exists)then
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(yfile),status='old')
       do it=1,itmax
          read(ftn,*,end=10)ztalt(it),zhalt(it)
       enddo
10     continue
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) &
            '[ERR] HYMAP3_read_time_series: Cannot read '//trim(yfile)
       call LIS_endrun()
    endif

  end subroutine HYMAP3_read_time_series
  !=============================================
  !ag(27Apr2020)
  ! ================================================
  !BOP
  !
  ! !ROUTINE: HYMAP3_gather_tiles
  ! \label{HYMAP3_gather_tiles}
  !
  ! !INTERFACE:
  subroutine HYMAP3_gather_tiles_int(n, var, var_glb)

! !USES:
    use LIS_coreMod
    use LIS_mpiMod
    use LIS_routingMod
!
! !DESCRIPTION:
!  This subroutine gathers an individual variable
!  across different processors into a global array
!EOP

    implicit none

    integer, intent(in)     :: n
    integer, intent(in)     :: var(LIS_rc%nroutinggrid(n))
    integer, intent(out)    :: var_glb(LIS_rc%glbnroutinggrid(n))

    integer        :: tmpvar(LIS_rc%glbnroutinggrid(n))
    integer        :: i,l,ix,iy,ix1,iy1
    integer        :: status

#if (defined SPMD)
    external :: MPI_ALLGATHERV
#endif

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
          ix = HYMAP3_routing_struc(n)%seqx_glb(i+&
               LIS_routing_goffsets(n,l-1))
          iy = HYMAP3_routing_struc(n)%seqy_glb(i+&
               LIS_routing_goffsets(n,l-1))
          ix1 = ix + LIS_ews_halo_ind(n,l) - 1
          iy1 = iy + LIS_nss_halo_ind(n,l)-1
          var_glb(HYMAP3_routing_struc(n)%sindex(ix1,iy1)) = &
               tmpvar(i+LIS_routing_goffsets(n,l-1))
       enddo
    enddo

  end subroutine HYMAP3_gather_tiles_int
! ================================================
!BOP
!
! !ROUTINE: HYMAP3_gather_tiles
! \label{HYMAP3_gather_tiles}
!
! !INTERFACE:
  subroutine HYMAP3_gather_tiles(n,var,var_glb)
! !USES:
    use LIS_coreMod
    use LIS_mpiMod
    use LIS_routingMod
!
! !DESCRIPTION:
!  This subroutine gathers an individual variable
!  across different processors into a global array
!EOP

    implicit none

    integer, intent(in)        :: n
    real, intent(in)           :: var(LIS_rc%nroutinggrid(n))
    real, intent(out)          :: var_glb(LIS_rc%glbnroutinggrid(n))

    real           :: tmpvar(LIS_rc%glbnroutinggrid(n))
    integer        :: i,l,ix,iy,ix1,iy1
    integer        :: status

#if (defined SPMD)
    external :: MPI_ALLGATHERV
#endif

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
          ix = HYMAP3_routing_struc(n)%seqx_glb(i+&
               LIS_routing_goffsets(n,l-1))
          iy = HYMAP3_routing_struc(n)%seqy_glb(i+&
               LIS_routing_goffsets(n,l-1))
          ix1 = ix + LIS_ews_halo_ind(n,l) - 1
          iy1 = iy + LIS_nss_halo_ind(n,l)-1
          var_glb(HYMAP3_routing_struc(n)%sindex(ix1,iy1)) = &
               tmpvar(i+LIS_routing_goffsets(n,l-1))
       enddo
    enddo

  end subroutine HYMAP3_gather_tiles
! ================================================
!BOP
! !ROUTINE: HYMAP3_map_g2l
! \label{HYMAP3_map_g2l}
!
! !INTERFACE:
  subroutine HYMAP3_map_g2l(n, var_glb,var_local)
! !USES:
    use LIS_coreMod

! !DESCRIPTION:
! This subroutine maps a global array in the HYMAP3
! tile space to the local processor space.
!
!EOP
    implicit none

    integer, intent(in) :: n
    real, intent(in)    :: var_glb(LIS_rc%glbnroutinggrid(n))
    real, intent(out)   :: var_local(LIS_rc%nroutinggrid(n))

    integer             :: i, ix,iy,ix1,iy1

    do i=1,LIS_rc%nroutinggrid(n)
       ix = HYMAP3_routing_struc(n)%seqx(i)
       iy = HYMAP3_routing_struc(n)%seqy(i)
       ix1 = ix + LIS_ews_halo_ind(n,LIS_localPet+1) -1
       iy1 = iy + LIS_nss_halo_ind(n,LIS_localPet+1) -1
       var_local(i)  = var_glb(HYMAP3_routing_struc(n)%sindex(ix1,iy1))
    enddo

  end subroutine HYMAP3_map_g2l
! ================================================
!BOP
! !ROUTINE: HYMAP3_map_g2l_index
! \label{HYMAP3_map_g2l_index}
!
! !INTERFACE:
  subroutine HYMAP3_map_gxy2l_index(n, glb_x, glb_y, local_index)

! !USES:
    use LIS_coreMod

!
! !DESCRIPTION:
! This subroutine converts x,y from a global array in the HYMAP3
! tile space to the index in the local processor space.
!
!EOP
    implicit none

    integer, intent(in)  :: n
    integer, intent(out) :: local_index
    integer, intent(in)  :: glb_x, glb_y

    integer              :: ix,iy,iloc(1)

    ix = glb_x - LIS_ews_halo_ind(n,LIS_localPet+1) + 1
    iy = glb_y - LIS_nss_halo_ind(n,LIS_localPet+1) + 1

    if(minval(abs(ix-HYMAP3_routing_struc(n)%seqx) + &
         abs(iy-HYMAP3_routing_struc(n)%seqy))==0)then
       iloc = minloc(abs(ix-HYMAP3_routing_struc(n)%seqx) + &
            abs(iy-HYMAP3_routing_struc(n)%seqy))
       local_index=iloc(1)
    else
       local_index = HYMAP3_routing_struc(n)%imis
    endif

  end subroutine HYMAP3_map_gxy2l_index
! ================================================
!BOP
! !ROUTINE: HYMAP3_map_l2g_index
! \label{HYMAP3_map_l2g_index}
!
! !INTERFACE:
  subroutine HYMAP3_map_l2g_index(n, local_index, glb_index)
! !USES:
    use LIS_coreMod
!
! !DESCRIPTION:
!  This subroutine converts the local tile index into the
!  the global index.
!
!EOP
    implicit none

    integer, intent(in)             :: n
    integer, intent(in)             :: local_index
    integer, intent(out)            :: glb_index

    integer             :: ix,iy,ix1,iy1

    ix = HYMAP3_routing_struc(n)%seqx(local_index)
    iy = HYMAP3_routing_struc(n)%seqy(local_index)
    ix1 = ix + LIS_ews_halo_ind(n,LIS_localPet+1) -1
    iy1 = iy + LIS_nss_halo_ind(n,LIS_localPet+1) -1
    glb_index = HYMAP3_routing_struc(n)%sindex(ix1,iy1)

  end subroutine HYMAP3_map_l2g_index
! ================================================

end module HYMAP3_routingMod
