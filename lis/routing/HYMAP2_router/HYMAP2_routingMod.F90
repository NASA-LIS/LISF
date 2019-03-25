!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
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
!             large-scale river flow hydraulics in the Amazon Basin, Water!             Resour. Res., 49, doi:10.1002/wrcr.20212
!
! !REVISION HISTORY: 
! 08 Nov 2011: Augusto Getirana, Initial implementation in LIS based on the 
!                                HYMAP offline routing code. 
! 19 Jan 2016: Augusto Getirana, Inclusion of the local inertia formulation, 
!                                adaptive time step and reservoir operation. 
! 13 Apr 2016: Augusto Getirana, Inclusion of option for hybrid runs with a 
!                                river flow map. 
! 
! !USES: 
  use ESMF
  use LIS_topoMod
  
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
  !integer              :: nseqriv       !length of 1D sequnece for river
  integer              :: nseqall       !length of 1D sequnece for river and mouth  
  integer              :: nseqall_glb   !global length of 1D sequnece for river and mouth  
  integer, allocatable :: gdeltas(:)
  integer, allocatable :: goffsets(:)
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
  real,    allocatable :: fldgrd(:,:)    !floodplain gradient [-]
  real,    allocatable :: fldman(:)      !maning coefficient for floodplains [-]
  real,    allocatable :: fldhgt(:,:)    !floodplain height [m]
  real,    allocatable :: trnoff(:)      !runoff   concentration time parameter [day]
  real,    allocatable :: tbsflw(:)      !baseflow concentration time parameter [day]
  real,    allocatable :: cntime(:)      !concentration time [s]
  real,    allocatable :: fldstomax(:,:) !maximum floodplain storage [m3]
  real,    allocatable :: rivman(:)      !maning coefficient for rivers [-]
  real,    allocatable :: rivelv(:)      !elevation of river bed [m]
  real,    allocatable :: rivstomax(:)   !maximum river storage [m3]
  real,    allocatable :: rivlen(:)      !river length [m] 
  real,    allocatable :: rivwth(:)      !river width [m]
  real,    allocatable :: rivhgt(:)        !river heihgt [m]
  real                 :: rslpmin          !minimum slope [-]
  real,    allocatable :: rivare(:)        !river surface area [m2]
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
  character*100        :: rstfile
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

  character*50         :: LISdir     !if LIS output is being read
     
! === Local inertia variables =====================
  real,    allocatable  :: rivout_pre(:,:)
  real,    allocatable  :: rivdph_pre(:,:)
  real,    allocatable  :: fldout_pre(:,:)
  real,    allocatable  :: flddph_pre(:,:)
  real,    allocatable  :: fldelv1(:,:)
  real                  :: cadp
  real                  :: dstend
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
  character*50         :: resopdir
  character*50         :: resopheader
  !ag (29Jun2016)
  integer              :: floodflag
  character*100        :: HYMAP_dfile      

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
    use LIS_mpiMod
    
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
    character(100) :: diag_fname
    character(100) :: diag_dir
    integer, external  :: LIS_create_subdirs

    !ag (03Apr2017)
    character*20 :: evapflag0
    !character*20 :: pevap_comp_method
    integer       :: ic, c,r,ic_down
    integer       :: gdeltas
    
    allocate(HYMAP2_routing_struc(LIS_rc%nnest))
    
    HYMAP2_logunit=10001
    
    do n=1, LIS_rc%nnest
       HYMAP2_routing_struc(n)%rslpmin  = 1e-5  !minimum slope
       HYMAP2_routing_struc(n)%nz      = 10    !number of stages in the sub-grid discretization
       HYMAP2_routing_struc(n)%imis     = -9999 !undefined integer value
       !HYMAP2_routing_struc(n)%cadp     = 0.7   !alfa coefficient for adaptative time step as described in Bates et al., (2010) [-]
       HYMAP2_routing_struc(n)%dstend   = 25000 !river length to the ocean [m]
       HYMAP2_routing_struc(n)%grv      = 9.81  !gravity accerelation [m/s2]
       HYMAP2_routing_struc(n)%numout   = 0 
       HYMAP2_routing_struc(n)%fileopen = 0 
       HYMAP2_routing_struc(n)%dt_proc  = 0.
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
       else
         write(LIS_logunit,*)"[ERR] HYMAP2 routing method: unknown value ",trim(flowtype)
         stop
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
         stop          
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
         stop
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
     
    !ag (4Feb2016)
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP2 reservoir operation option:",rc=status)
    do n=1, LIS_rc%nnest
      call ESMF_ConfigGetAttribute(LIS_config,&
           HYMAP2_routing_struc(n)%resopflag,rc=status)
      call LIS_verify(status,&
           "HYMAP2 reservoir operation option: not defined")

      if(HYMAP2_routing_struc(n)%resopflag==1)then
        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 number of reservoirs:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%nresop,rc=status)
        call LIS_verify(status,&
             "HYMAP2 number of reservoirs: not defined")

        call ESMF_ConfigFindLabel(LIS_config,&
             "HYMAP2 reservoir operation input time series size:",rc=status)
        call ESMF_ConfigGetAttribute(LIS_config,&
             HYMAP2_routing_struc(n)%ntresop,rc=status)
        call LIS_verify(status,&
             "HYMAP2 reservoir operation input time series size: not defined")

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
        else
          write(LIS_logunit,*)"[ERR] HYMAP2 routing model time step method: unknown value"
          stop          
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
      
       
       write(LIS_logunit,*) '[INFO] Get number of HYMAP2 grid cells'
       call HYMAP2_get_vector_size(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            LIS_rc%gnc(n),LIS_rc%gnr(n),&
            LIS_ews_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1), &
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%nextx,&
            mask,HYMAP2_routing_struc(n)%nseqall)

       allocate(HYMAP2_routing_struc(n)%gdeltas(0:LIS_npes-1))
       allocate(HYMAP2_routing_struc(n)%goffsets(0:LIS_npes-1))

       gdeltas = HYMAP2_routing_struc(n)%nseqall      
       
#if (defined SPMD)
       call MPI_ALLREDUCE(HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%nseqall_glb,1,&
            MPI_INTEGER,MPI_SUM,&
            LIS_mpi_comm,status)

       call MPI_ALLGATHER(gdeltas,1,MPI_INTEGER,&
            HYMAP2_routing_struc(n)%gdeltas(:),1,MPI_INTEGER,&
            LIS_mpi_comm,status)

#else
       HYMAP2_routing_struc(n)%nseqall_glb = HYMAP2_routing_struc(n)%nseqall
#endif

       if(LIS_masterproc) then 
          HYMAP2_routing_struc(n)%goffsets = 0 
          do i=1,LIS_npes-1
             HYMAP2_routing_struc(n)%goffsets(i) = &
                  HYMAP2_routing_struc(n)%goffsets(i-1) +& 
                  HYMAP2_routing_struc(n)%gdeltas(i-1)
          enddo
       endif
#if (defined SPMD)
       call MPI_BCAST(HYMAP2_routing_struc(n)%goffsets, &
            LIS_npes, MPI_INTEGER,0, &
            LIS_mpi_comm, status)
#endif
       allocate(HYMAP2_routing_struc(n)%seqx(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%seqy(HYMAP2_routing_struc(n)%nseqall))

       allocate(HYMAP2_routing_struc(n)%seqx_glb(HYMAP2_routing_struc(n)%nseqall_glb))
       allocate(HYMAP2_routing_struc(n)%seqy_glb(HYMAP2_routing_struc(n)%nseqall_glb))

       allocate(HYMAP2_routing_struc(n)%sindex(LIS_rc%gnc(n),LIS_rc%gnr(n)))
       allocate(HYMAP2_routing_struc(n)%outlet(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%outlet_glb(HYMAP2_routing_struc(n)%nseqall_glb))
       allocate(HYMAP2_routing_struc(n)%next(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%next_glb(HYMAP2_routing_struc(n)%nseqall_glb))
       !allocate(HYMAP2_routing_struc(n)%nextx(HYMAP2_routing_struc(n)%nseqall))
       !allocate(HYMAP2_routing_struc(n)%nexty(HYMAP2_routing_struc(n)%nseqall)) 
       !allocate(HYMAP2_routing_struc(n)%mask(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%elevtn(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%nxtdst(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%grarea(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%fldgrd(HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%nz))
       allocate(HYMAP2_routing_struc(n)%fldman(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%fldhgt(HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%nz))
       allocate(HYMAP2_routing_struc(n)%cntime(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%trnoff(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%tbsflw(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%fldstomax(HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%nz))
       allocate(HYMAP2_routing_struc(n)%rivman(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rivelv(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rivstomax(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rivlen(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rivwth(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rivhgt(HYMAP2_routing_struc(n)%nseqall))
       allocate(HYMAP2_routing_struc(n)%rivare(HYMAP2_routing_struc(n)%nseqall))
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
       
    enddo
    
    do n=1,LIS_rc%nnest
       allocate(tmp_real(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(tmp_real_nz(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            HYMAP2_routing_struc(n)%nz))
    enddo

    write(LIS_logunit,*) '[INFO] Get HYMAP2 cell vector sequence'
    do n=1, LIS_rc%nnest
       call HYMAP2_get_sindex(LIS_rc%gnc(n),&
            LIS_rc%gnr(n),&
            HYMAP2_routing_struc(n)%nseqall_glb,&
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
            HYMAP2_routing_struc(n)%nseqall,MPI_INTEGER,&
            HYMAP2_routing_struc(n)%seqx_glb,&
            HYMAP2_routing_struc(n)%gdeltas(:),&
            HYMAP2_routing_struc(n)%goffsets(:),&
            MPI_INTEGER,&
            LIS_mpi_comm,status)

       call MPI_ALLGATHERV(HYMAP2_routing_struc(n)%seqy,&
            HYMAP2_routing_struc(n)%nseqall,MPI_INTEGER,&
            HYMAP2_routing_struc(n)%seqy_glb,&
            HYMAP2_routing_struc(n)%gdeltas(:),&
            HYMAP2_routing_struc(n)%goffsets(:),&
            MPI_INTEGER,&
            LIS_mpi_comm,status)
    enddo
#endif

    ctitle = 'HYMAP_river_width'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%rivwth)
    enddo
    

    ctitle = 'HYMAP_river_length'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%rivlen)
    enddo
    
    ctitle = 'HYMAP_river_height'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%rivhgt)  
    enddo
    
    ctitle = 'HYMAP_river_roughness'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%rivman)  
    enddo
    
    ctitle = 'HYMAP_floodplain_height'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real(ctitle,n, HYMAP2_routing_struc(n)%nz,&
            tmp_real_nz)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            HYMAP2_routing_struc(n)%nz,&
            HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real_nz,HYMAP2_routing_struc(n)%fldhgt)
    enddo

    ctitle = 'HYMAP_floodplain_roughness'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%fldman)
    enddo

    ctitle = 'HYMAP_grid_elevation'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       where(tmp_real<0)tmp_real=0
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%elevtn)
    enddo

    ctitle = 'HYMAP_grid_distance'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%nxtdst)
       where(HYMAP2_routing_struc(n)%outlet==1)HYMAP2_routing_struc(n)%nxtdst=HYMAP2_routing_struc(n)%dstend
    enddo

    ctitle = 'HYMAP_grid_area'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%grarea)
    enddo

    ctitle = 'HYMAP_runoff_time_delay'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%cntime)
       where(HYMAP2_routing_struc(n)%cntime==0)HYMAP2_routing_struc(n)%cntime=minval(HYMAP2_routing_struc(n)%cntime,HYMAP2_routing_struc(n)%cntime>0)
    enddo


    ctitle = 'HYMAP_runoff_time_delay_multiplier'
    do n=1, LIS_rc%nnest
       call HYMAP2_read_param_real_2d(ctitle,n,tmp_real)
       call HYMAP2_grid2vector(LIS_rc%lnc(n),LIS_rc%lnr(n),&
            1,HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%imis,&
            HYMAP2_routing_struc(n)%seqx,&
            HYMAP2_routing_struc(n)%seqy,&
            tmp_real,HYMAP2_routing_struc(n)%trnoff)	  
    enddo

    ctitle = 'HYMAP_baseflow_time_delay'
    do n=1, LIS_rc%nnest
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
    enddo

    !ag (23Nov2016)
    ctitle = 'HYMAP_runoff_dwi_ratio'
    do n=1, LIS_rc%nnest
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
    enddo
    !ag (23Nov2016)
    ctitle = 'HYMAP_baseflow_dwi_ratio'
    do n=1, LIS_rc%nnest
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
    enddo

    !ag (13Apr2016)
    ctitle = 'HYMAP_river_flow_type'
    do n=1, LIS_rc%nnest
       if(HYMAP2_routing_struc(n)%flowtype==0)then
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

    enddo

    !ag (20Sep2016) Correction for cases where parameter maps don't match
    do n=1, LIS_rc%nnest
       do i=1,HYMAP2_routing_struc(n)%nseqall
          if(HYMAP2_routing_struc(n)%rivwth(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%rivlen(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%rivhgt(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%rivman(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%fldhgt(i,1)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%elevtn(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%nxtdst(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%grarea(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%cntime(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%trnoff(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%tbsflw(i)==HYMAP2_routing_struc(n)%imis.or.&
               HYMAP2_routing_struc(n)%flowmap(i)==HYMAP2_routing_struc(n)%imis)then

             HYMAP2_routing_struc(n)%rivwth(i)=HYMAP2_routing_struc(n)%imis
             HYMAP2_routing_struc(n)%rivlen(i)=HYMAP2_routing_struc(n)%imis
             HYMAP2_routing_struc(n)%rivhgt(i)=HYMAP2_routing_struc(n)%imis
             HYMAP2_routing_struc(n)%rivman(i)=HYMAP2_routing_struc(n)%imis
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


    write(LIS_logunit,*) '[INFO] Processing data before running HYMAP'
    do n=1, LIS_rc%nnest
       write(LIS_logunit,*)'[INFO] Calculate maximum river storage'
       HYMAP2_routing_struc(n)%rivstomax = HYMAP2_routing_struc(n)%rivlen* &
            HYMAP2_routing_struc(n)%rivwth * HYMAP2_routing_struc(n)%rivhgt
       write(LIS_logunit,*)'[INFO] Calculate river bed elevation'
       HYMAP2_routing_struc(n)%rivelv = HYMAP2_routing_struc(n)%elevtn -&
            HYMAP2_routing_struc(n)%rivhgt
       write(LIS_logunit,*)'[INFO] Calculate river surface area'
       where(HYMAP2_routing_struc(n)%rivwth>0)HYMAP2_routing_struc(n)%rivare =&
            min(HYMAP2_routing_struc(n)%grarea, HYMAP2_routing_struc(n)%rivlen *&
            HYMAP2_routing_struc(n)%rivwth)
       write(LIS_logunit,*)'[INFO] Setting floodplain staging'
       call HYMAP2_set_fldstg(HYMAP2_routing_struc(n)%nz,&
            HYMAP2_routing_struc(n)%nseqall,&
            HYMAP2_routing_struc(n)%fldhgt,&
            HYMAP2_routing_struc(n)%grarea,&
            HYMAP2_routing_struc(n)%rivlen,&
            HYMAP2_routing_struc(n)%rivwth,&
            HYMAP2_routing_struc(n)%rivstomax,&
            HYMAP2_routing_struc(n)%fldstomax,&
            HYMAP2_routing_struc(n)%fldgrd,&
            HYMAP2_routing_struc(n)%rivare)		   
       !Start storages
       HYMAP2_routing_struc(n)%rivsto=0.0
       HYMAP2_routing_struc(n)%fldsto=0.0
       HYMAP2_routing_struc(n)%rnfsto=0.0
       HYMAP2_routing_struc(n)%bsfsto=0.0
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
            'nf90_open failed in read_param_real@HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_real@HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1),z/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_real@HYMAP2_routingMod')
       
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routindMod')
       
!       do l=1,z
!          array(:,:,l) = nint(array1(&
!               LIS_ews_halo_ind(n,LIS_localPet+1):&         
!               LIS_ewe_halo_ind(n,LIS_localPet+1), &
!               LIS_nss_halo_ind(n,LIS_localPet+1): &
!               LIS_nse_halo_ind(n,LIS_localPet+1),l))
!       enddo
    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_real@HYMAP2_routingMod'
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
            'nf90_open failed in read_param_real@HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_real@HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1)/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1)/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_real@HYMAP2_routingMod')
       
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routindMod')
       
    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_real@HYMAP2_routingMod'
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
            'nf90_open failed in read_param_int@HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int@HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),1/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1),z/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int@HYMAP2_routingMod')
       
       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routindMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_int@HYMAP2_routingMod'
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
            'nf90_open failed in read_param_int@HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int@HYMAP2_routingMod')
       

       call LIS_verify(nf90_get_var(ftn,varid, array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1)/), &
            count=(/(LIS_ewe_halo_ind(n,LIS_localPet+1)-&
            LIS_ews_halo_ind(n,LIS_localPet+1)+1),&
            (LIS_nse_halo_ind(n,LIS_localPet+1)-&
           LIS_nss_halo_ind(n,LIS_localPet+1)+1)/)),&
           'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int@HYMAP2_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routindMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_int@HYMAP2_routingMod'
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
            'nf90_open failed in read_param_int@HYMAP2_routingMod')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in read_param_int@HYMAP2_routingMod')
       
       call LIS_verify(nf90_get_var(ftn,varid, array), &
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in read_param_int@HYMAP2_routingMod')

       call LIS_verify(nf90_close(ftn),&
            'nf90_close failed in HYMAP2_routindMod')

    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in read_param_int@HYMAP2_routingMod'
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
    character(500)            :: yfile
    
    call HYMAP2_read_header_resop(trim(resopheader),nresop,resopname,xresop,yresop,resopoutmin)
    do res=1,nresop
      resoploc(res)=sindex(xresop(res),yresop(res))
    enddo
    do res=1,nresop
      yfile=trim(resopdir)//trim(resopname(res))//'.txt'
      call HYMAP2_read_resop_alt(ntresop,trim(yfile),tresop(res,:),resop(res,:))
    enddo
  end subroutine HYMAP2_get_data_resop_alt
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
      write(LIS_logunit,*) 'failed in read_header@HYMAP2_routingMod'
      call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) 'header file '//trim(yheader)
    write(LIS_logunit,*) 'failed in read_header@HYMAP2_routingMod'
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
      write(LIS_logunit,*) 'failed in read_header@HYMAP2_routingMod'
      call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) 'header file '//trim(yheader)
    write(LIS_logunit,*) 'failed in read_header@HYMAP2_routingMod'
    call LIS_endrun()

  end subroutine HYMAP2_read_header
  !=============================================
  !=============================================  
  subroutine HYMAP2_read_resop_alt(itmax,yfile,ztalt,zhalt)
  
    use LIS_logMod
    implicit none

    integer,      intent(in)    :: itmax
    character(*), intent(in)    :: yfile
    real*8,       intent(inout) :: ztalt(itmax)
    real,         intent(inout) :: zhalt(itmax)
    logical                     :: file_exists
    integer :: it
  
    inquire(file=yfile,exist=file_exists)
    
    if(file_exists)then 
      open(1,file=trim(yfile),status='old')
      do it=1,itmax
        read(1,*,end=10)ztalt(it),zhalt(it)
        write(LIS_logunit,*)ztalt(it),zhalt(it)
      enddo
10    continue
      close(1)
    else
      write(LIS_logunit,*) 'reservoir data time series file '//trim(yfile)
      write(LIS_logunit,*) 'failed in read_resop_alt@HYMAP2_routing_init'
      call LIS_endrun()
    endif  

  end subroutine HYMAP2_read_resop_alt
  !=============================================
  !=============================================  
end module HYMAP2_routingMod
