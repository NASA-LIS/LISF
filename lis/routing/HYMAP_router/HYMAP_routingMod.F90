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
module HYMAP_routingMod
!BOP
! 
! !MODULE: HYMAP_routingMod
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
!
! !REVISION HISTORY: 
! 03 Nov 2011: Augusto Getirana, Initial implementation in LIS based on the 
!                                HYMAP offline routing code. 
! 
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: HYMAP_routingInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: HYMAP_routing_struc
  
  type, public :: HYMAP_routing_dec
     
     real             :: dt
     integer          :: useens
     integer          :: inz            !! number of stages in the sub-grid discretization
! === Undefined Values ==================================
     integer          :: imis           !! real undefined value
! === River sequence ====================================
     integer, allocatable ::  seqx(:)       !! 1D sequence horizontal
     integer, allocatable ::  seqy(:)       !! 1D sequence vertical
     integer          ::  nseqriv       !! length of 1D sequnece for river
     integer          ::  nseqall       !! length of 1D sequnece for river and mouth  
! === Map ===============================================
     integer, allocatable :: nextx(:,:)       !! point downstream horizontal
     integer, allocatable :: nexty(:,:)       !! point downstream vertical
     integer, allocatable :: mask(:,:)        !! mask limiting modeled region (0: out; >=1: in)
     real,    allocatable :: elevtn(:,:)      !! river bed elevation [m]
     real,    allocatable :: nxtdst(:,:)      !! distance to the next grid [m]
     real,    allocatable :: grarea(:,:)      !! area of the grid [m^2] 
     real,    allocatable :: fldgrd(:,:,:)    !! floodplain gradient [-]
     real,    allocatable :: fldman(:,:)      !! maning coefficient for floodplains [-]
     real,    allocatable :: fldhgt(:,:,:)    !! floodplain height [m]
     real,    allocatable :: trnoff(:,:)      !! runoff   concentration time parameter [day]
     real,    allocatable :: tbsflw(:,:)      !! baseflow concentration time parameter [day]
     real,    allocatable :: cntime(:,:)      !! concentration time [s]
     real,    allocatable :: fldstomax(:,:,:) !! maximum floodplain storage [m3]
     real,    allocatable :: rivman(:,:)      !! maning coefficient for rivers [-]
     real,    allocatable :: rivelv(:,:)      !! elevation of river bed [m]
     real,    allocatable :: rivstomax(:,:)   !! maximum river storage [m3]
     real,    allocatable :: rivlen(:,:)      !! river length [m] 
     real,    allocatable :: rivwth(:,:)      !! river width [m]
     real,    allocatable :: rivhgt(:,:)      !! river heihgt [m]
     real             :: rslpmin          !! minimum slope [-]
     real,    allocatable :: rivare(:,:)      !! river surface area [m2]
     real,    allocatable :: runoff0(:,:)     !! input runoff [mm.dt-1]
     real,    allocatable :: basflw0(:,:)     !! input baseflow [mm.dt-1]
 ! === Outputs ==========================================
     real,    allocatable :: rivsto(:,:,:)      !! river storage [m3]
     real,    allocatable :: rivdph(:,:,:)      !! river depth [m]
     real,    allocatable :: rivvel(:,:,:)      !! river flow velocity [m/s]
     real,    allocatable :: rivout(:,:,:)      !! river outflow  [m3/s]
     real,    allocatable :: evpout(:,:,:)      !! evaporation from open waters [m3]
     real,    allocatable :: fldout(:,:,:)      !! floodplain outflow  [m3/s]
     real,    allocatable :: fldsto(:,:,:)      !! floodplain storage   [m3]
     real,    allocatable :: flddph(:,:,:)      !! floodplain depth [m]
     real,    allocatable :: fldvel(:,:,:)      !! floodplain flow velocity [m/s]
     real,    allocatable :: fldfrc(:,:,:)      !! area fraction [-]
     real,    allocatable :: fldare(:,:,:)      !! flooded area [m2]
     real,    allocatable :: sfcelv(:,:,:)      !! water surface elevation [m]
! === Linear reservoir components ================
     real,    allocatable :: rnfsto(:,:,:)      !! runoff   reservoir storage [m3]
     real,    allocatable :: bsfsto(:,:,:)      !! baseflow reservoir storage [m3]
     character(len=LIS_CONST_PATH_LEN) :: rstfile
     integer          :: numout
     integer          :: fileopen
     real             :: outInterval 
     real             :: rstInterval
     character*20     :: startMode
     integer          :: mo
     !ag (04Jun2012)
     integer          :: flowtype
     integer          :: linresflag
     integer          :: evapflag
  end type HYMAP_routing_dec

  type(HYMAP_routing_dec), allocatable :: HYMAP_routing_struc(:)

contains
 
!BOP
!
! !ROUTINE: HYMAP_routingInit
! \label{HYMAP_routingInit}
! 
  subroutine HYMAP_routingInit
    !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_routingMod    
    use HYMAP_initMod
    
    integer              :: n 
    integer              :: ftn 
    integer              :: status
    integer              :: i,j,k
    type(ESMF_Grid)      :: global_grid
    type(ESMF_DistGrid)  :: global_grid_dg
    type(ESMF_ArraySpec) :: realarrspec
    type(ESMF_Field)     :: sf_runoff_field
    type(ESMF_Field)     :: baseflow_field
    real, pointer        :: sfrunoff(:)
    real, pointer        :: baseflow(:)
    character*100        :: ctitle
    integer              :: iseq,ix,iy,jx,jy
    character*10         :: time
    character*20         :: flowtype

    allocate(HYMAP_routing_struc(LIS_rc%nnest))
    
    do n=1, LIS_rc%nnest
       HYMAP_routing_struc(n)%rslpmin  = 1e-5  !! minimum slope
       HYMAP_routing_struc(n)%inz      = 10    !! number of stages in the sub-grid discretization
       HYMAP_routing_struc(n)%imis     = -9999 !! undefined integer value
       HYMAP_routing_struc(n)%numout   = 0 
       HYMAP_routing_struc(n)%fileopen = 0 
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP routing model time step:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP routing model time step: not defined")

       call LIS_parseTimeString(time,HYMAP_routing_struc(n)%dt)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP routing model output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "HYMAP routing model output interval: not defined")

       call LIS_parseTimeString(time,HYMAP_routing_struc(n)%outInterval)
    enddo
    
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP run in ensemble mode:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP_routing_struc(n)%useens,rc=status)
       call LIS_verify(status,&
            "HYMAP run in ensemble mode: not defined")
    enddo


    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP routing method:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            flowtype,rc=status)
       call LIS_verify(status,&
            "HYMAP routing method: not defined")
       if(flowtype.eq."kinematic") then 
          HYMAP_routing_struc(n)%flowtype = 1
       elseif(flowtype.eq."diffusive") then 
          HYMAP_routing_struc(n)%flowtype = 2
       endif
    enddo
    
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP routing model linear reservoir flag:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP_routing_struc(n)%linresflag,rc=status)
       call LIS_verify(status,&
            "HYMAP routing model linear reservoir flag: not defined")
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "HYMAP routing model evaporation option:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            HYMAP_routing_struc(n)%evapflag,rc=status)
       call LIS_verify(status,&
            "HYMAP routing model evaporation option: not defined")
    enddo
    
     
    if(LIS_masterproc) then 
       write(LIS_logunit,*) '[INFO] Initializing HYMAP....'
       !allocate matrixes
       do n=1, LIS_rc%nnest

          allocate(HYMAP_routing_struc(n)%seqx(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%seqy(LIS_rc%gnc(n)*LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%nextx(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%nexty(LIS_rc%gnc(n),LIS_rc%gnr(n))) 
          allocate(HYMAP_routing_struc(n)%mask(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%elevtn(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%nxtdst(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%grarea(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%fldgrd(LIS_rc%gnc(n),LIS_rc%gnr(n),&
               HYMAP_routing_struc(n)%inz))
          allocate(HYMAP_routing_struc(n)%fldman(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%fldhgt(LIS_rc%gnc(n),LIS_rc%gnr(n),&
               HYMAP_routing_struc(n)%inz))
          allocate(HYMAP_routing_struc(n)%cntime(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%trnoff(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%tbsflw(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%fldstomax(LIS_rc%gnc(n),LIS_rc%gnr(n),&
               HYMAP_routing_struc(n)%inz))
          allocate(HYMAP_routing_struc(n)%rivman(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%rivelv(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%rivstomax(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%rivlen(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%rivwth(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%rivhgt(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%rivare(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%runoff0(LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(HYMAP_routing_struc(n)%basflw0(LIS_rc%gnc(n),LIS_rc%gnr(n)))

          if(HYMAP_routing_struc(n)%useens.eq.0) then 
             allocate(HYMAP_routing_struc(n)%rivsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%rivdph(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%rivvel(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%rivout(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%evpout(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%fldout(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%fldsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%flddph(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%fldvel(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%fldfrc(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%fldare(LIS_rc%gnc(n),LIS_rc%gnr(n),&
               1))
             allocate(HYMAP_routing_struc(n)%sfcelv(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%rnfsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
             allocate(HYMAP_routing_struc(n)%bsfsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  1))
          else
             allocate(HYMAP_routing_struc(n)%rivsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%rivdph(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%rivvel(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%rivout(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%evpout(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%fldout(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%fldsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%flddph(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%fldvel(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%fldfrc(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%fldare(LIS_rc%gnc(n),LIS_rc%gnr(n),&
               LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%sfcelv(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%rnfsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
             allocate(HYMAP_routing_struc(n)%bsfsto(LIS_rc%gnc(n),LIS_rc%gnr(n),&
                  LIS_rc%nensem(n)))
          endif

          HYMAP_routing_struc(n)%seqx=0.0
          HYMAP_routing_struc(n)%seqy=0.0
          HYMAP_routing_struc(n)%nextx=0.0
          HYMAP_routing_struc(n)%nexty=0.0 
          HYMAP_routing_struc(n)%mask=0.0
          HYMAP_routing_struc(n)%elevtn=0.0
          HYMAP_routing_struc(n)%nxtdst=0.0
          HYMAP_routing_struc(n)%grarea=0.0
          HYMAP_routing_struc(n)%fldgrd=0.0
          HYMAP_routing_struc(n)%fldman=0.0
          HYMAP_routing_struc(n)%fldhgt=0.0
          HYMAP_routing_struc(n)%cntime=0.0
          HYMAP_routing_struc(n)%trnoff=0.0
          HYMAP_routing_struc(n)%tbsflw=0.0
          HYMAP_routing_struc(n)%fldstomax=0.0
          HYMAP_routing_struc(n)%rivman=0.0
          HYMAP_routing_struc(n)%rivelv=0.0
          HYMAP_routing_struc(n)%rivstomax=0.0
          HYMAP_routing_struc(n)%rivlen=0.0
          HYMAP_routing_struc(n)%rivwth=0.0
          HYMAP_routing_struc(n)%rivhgt=0.0
          HYMAP_routing_struc(n)%rivare=0.0
          HYMAP_routing_struc(n)%runoff0=0.0
          HYMAP_routing_struc(n)%basflw0=0.0
          HYMAP_routing_struc(n)%rivsto=0.0
          HYMAP_routing_struc(n)%rivdph=0.0
          HYMAP_routing_struc(n)%rivvel=0.0
          HYMAP_routing_struc(n)%rivout=0.0
          HYMAP_routing_struc(n)%evpout=0.0
          HYMAP_routing_struc(n)%fldout=0.0
          HYMAP_routing_struc(n)%fldsto=0.0
          HYMAP_routing_struc(n)%flddph=0.0
          HYMAP_routing_struc(n)%fldvel=0.0
          HYMAP_routing_struc(n)%fldfrc=0.0
          HYMAP_routing_struc(n)%fldare=0.0
          HYMAP_routing_struc(n)%sfcelv=0.0
          HYMAP_routing_struc(n)%rnfsto=0.0
          HYMAP_routing_struc(n)%bsfsto=0.0
       enddo

       ctitle = 'HYMAP_river_width'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%rivwth)
       enddo
       
       ctitle = 'HYMAP_river_length'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%rivlen)
       enddo
       
       ctitle = 'HYMAP_river_height'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%rivhgt)  
       enddo
       
       ctitle = 'HYMAP_river_roughness'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%rivman)  
       enddo
       
       ctitle = 'HYMAP_floodplain_height'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,HYMAP_routing_struc(n)%inz,n,&
               HYMAP_routing_struc(n)%fldhgt)	  
       enddo
       
       ctitle = 'HYMAP_floodplain_roughness'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%fldman)
       enddo

       ctitle = 'HYMAP_flow_direction_x'
       do n=1, LIS_rc%nnest
          call read_param_int(ctitle,1,n,HYMAP_routing_struc(n)%nextx)	  
       enddo
       
       ctitle = 'HYMAP_flow_direction_y'
       do n=1, LIS_rc%nnest
          call read_param_int(ctitle,1,n,HYMAP_routing_struc(n)%nexty)	  
       enddo

       ctitle = 'HYMAP_grid_elevation'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%elevtn)
       enddo
       
       ctitle = 'HYMAP_grid_distance'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%nxtdst)
       enddo

       ctitle = 'HYMAP_grid_area'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%grarea)
       enddo
       
       ctitle = 'HYMAP_runoff_time_delay'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%cntime)
       enddo
       
       ctitle = 'HYMAP_runoff_time_delay_multiplier'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%trnoff)	  
       enddo
       
       ctitle = 'HYMAP_baseflow_time_delay'
       do n=1, LIS_rc%nnest
          call read_param_real(ctitle,1,n,HYMAP_routing_struc(n)%tbsflw)
       enddo

       ctitle = 'HYMAP_basin_mask'
       do n=1, LIS_rc%nnest
          call read_param_int(ctitle,1,n,HYMAP_routing_struc(n)%mask)	  
       enddo              
                                  
       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP routing model start mode:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP_routing_struc(n)%startmode,rc=status)
          call LIS_verify(status,&
               "HYMAP routing model start mode: not defined")
       enddo
       
       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP routing model restart interval:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
          call LIS_verify(status,&
               "HYMAP routing model restart interval: not defined")

          call LIS_parseTimeString(time,HYMAP_routing_struc(n)%rstInterval)
       enddo
       
       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP routing model restart file:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP_routing_struc(n)%rstfile,rc=status)
          call LIS_verify(status,&
               "HYMAP routing model restart file: not defined")
       enddo
       
       if(LIS_rc%lsm.eq."none") then 
          call initrunoffdata(trim(LIS_rc%runoffdatasource)//char(0))
       endif

       write(LIS_logunit,*) '[INFO] Processing data before running HYMAP'
       do n=1, LIS_rc%nnest
          
          call calc_1d_seq(LIS_rc%gnc(n),&
               LIS_rc%gnr(n),&
               HYMAP_routing_struc(n)%imis,& 
               HYMAP_routing_struc(n)%nextx,&
               HYMAP_routing_struc(n)%mask,&
               HYMAP_routing_struc(n)%seqx,&
               HYMAP_routing_struc(n)%seqy,&
               HYMAP_routing_struc(n)%nseqriv,&
               HYMAP_routing_struc(n)%nseqall)
          
          write(LIS_logunit,*)'[INFO] Calculate maximum river storage'
          HYMAP_routing_struc(n)%rivstomax = HYMAP_routing_struc(n)%rivlen* &
               HYMAP_routing_struc(n)%rivwth * HYMAP_routing_struc(n)%rivhgt
          
          do iseq=1, HYMAP_routing_struc(n)%nseqriv
             ix=HYMAP_routing_struc(n)%seqx(iseq)
             iy=HYMAP_routing_struc(n)%seqy(iseq)
             jx=HYMAP_routing_struc(n)%nextx(ix,iy)
             jy=HYMAP_routing_struc(n)%nexty(ix,iy)              
          enddo

!ag - 11 Dec 2012: remove zeros from groundwater time delay for simulations using the linear_reservoir_lis routine
          write(LIS_logunit,*)'[INFO] Remove zeros from groundwater time delay'
          where(HYMAP_routing_struc(n)%tbsflw<HYMAP_routing_struc(n)%cntime/86400.)&
               HYMAP_routing_struc(n)%tbsflw=HYMAP_routing_struc(n)%cntime/86400.
          write(LIS_logunit,*)'[INFO] Calculate river bed elevation'
          HYMAP_routing_struc(n)%rivelv = HYMAP_routing_struc(n)%elevtn -&
               HYMAP_routing_struc(n)%rivhgt
          write(LIS_logunit,*)'[INFO] Calculate river surface area'
          where(HYMAP_routing_struc(n)%rivwth>0)HYMAP_routing_struc(n)%rivare =&
               min(HYMAP_routing_struc(n)%grarea, HYMAP_routing_struc(n)%rivlen *&
               HYMAP_routing_struc(n)%rivwth) 
          write(LIS_logunit,*)'[INFO] Setting floodplain staging'
!write(777,'(100f20.7)')HYMAP_routing_struc(n)%fldhgt(17,5,:)
!stop
          call set_fldstg(LIS_rc%gnc(n),&
                          LIS_rc%gnr(n),&
                          HYMAP_routing_struc(n)%inz,&
                          HYMAP_routing_struc(n)%seqx,&
                          HYMAP_routing_struc(n)%seqy,&
                          HYMAP_routing_struc(n)%nseqall,&
                          HYMAP_routing_struc(n)%fldhgt,&
                          HYMAP_routing_struc(n)%grarea,&
                          HYMAP_routing_struc(n)%rivlen,&
                          HYMAP_routing_struc(n)%rivwth,&
                          HYMAP_routing_struc(n)%rivstomax,&
                          HYMAP_routing_struc(n)%fldstomax,&
                          HYMAP_routing_struc(n)%fldgrd,&
                          HYMAP_routing_struc(n)%rivare)
          !Start storages
          HYMAP_routing_struc(n)%rivsto=0.0
          HYMAP_routing_struc(n)%fldsto=0.0
          HYMAP_routing_struc(n)%rnfsto=0.0
          HYMAP_routing_struc(n)%bsfsto=0.0
       enddo

       do n=1, LIS_rc%nnest
          call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
               rc=status)
          call LIS_verify(status)
          
          global_grid_dg = ESMF_DistGridCreate(minIndex=(/1/),&
               maxIndex=(/LIS_rc%glbntiles(n)/),&
               regDecomp=(/1/),rc=status)
          call LIS_verify(status,&
               'ESMF_DistGridCreate failed in HYMAP_routingInit')
          
          global_grid = ESMF_GridCreate(name="Global Grid",&
               coordTypeKind=ESMF_TYPEKIND_R4,&
               distgrid = global_grid_dg, &
               gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
          call LIS_verify(status,&
               'ESMF_GridCreate failed in HYMAP_routingInit')
          
          !create LSM interface objects to store runoff and baseflow
          sf_runoff_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=global_grid, name="Surface Runoff",rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')
          
          baseflow_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=global_grid, name="Subsurface Runoff",rc=status)
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

          HYMAP_routing_struc(n)%mo = -1
       enddo

    else
       do n=1, LIS_rc%nnest
          allocate(HYMAP_routing_struc(n)%seqx(1))
          allocate(HYMAP_routing_struc(n)%seqy(1))
          allocate(HYMAP_routing_struc(n)%nextx(1,1))
          allocate(HYMAP_routing_struc(n)%nexty(1,1)) 
          allocate(HYMAP_routing_struc(n)%mask(1,1))
          allocate(HYMAP_routing_struc(n)%elevtn(1,1))
          allocate(HYMAP_routing_struc(n)%nxtdst(1,1))
          allocate(HYMAP_routing_struc(n)%grarea(1,1))
          allocate(HYMAP_routing_struc(n)%fldgrd(1,1,&
                         HYMAP_routing_struc(n)%inz))
          allocate(HYMAP_routing_struc(n)%fldman(1,1))
          allocate(HYMAP_routing_struc(n)%fldhgt(1,1,&
                         HYMAP_routing_struc(n)%inz))
          allocate(HYMAP_routing_struc(n)%cntime(1,1))
          allocate(HYMAP_routing_struc(n)%trnoff(1,1))
          allocate(HYMAP_routing_struc(n)%tbsflw(1,1))
          allocate(HYMAP_routing_struc(n)%rnfsto(1,1,1))
          allocate(HYMAP_routing_struc(n)%bsfsto(1,1,1))
          allocate(HYMAP_routing_struc(n)%fldstomax(1,1,&
                         HYMAP_routing_struc(n)%inz))
          allocate(HYMAP_routing_struc(n)%rivman(1,1))
          allocate(HYMAP_routing_struc(n)%rivelv(1,1))
          allocate(HYMAP_routing_struc(n)%rivstomax(1,1))
          allocate(HYMAP_routing_struc(n)%rivlen(1,1))
          allocate(HYMAP_routing_struc(n)%rivwth(1,1))
          allocate(HYMAP_routing_struc(n)%rivhgt(1,1))
          allocate(HYMAP_routing_struc(n)%rivare(1,1))
          allocate(HYMAP_routing_struc(n)%runoff0(1,1))
          allocate(HYMAP_routing_struc(n)%basflw0(1,1))
          allocate(HYMAP_routing_struc(n)%rivsto(1,1,1))
          allocate(HYMAP_routing_struc(n)%rivdph(1,1,1))
          allocate(HYMAP_routing_struc(n)%rivvel(1,1,1))
          allocate(HYMAP_routing_struc(n)%rivout(1,1,1))
          allocate(HYMAP_routing_struc(n)%evpout(1,1,1))
          allocate(HYMAP_routing_struc(n)%fldout(1,1,1))
          allocate(HYMAP_routing_struc(n)%fldsto(1,1,1))
          allocate(HYMAP_routing_struc(n)%flddph(1,1,1))
          allocate(HYMAP_routing_struc(n)%fldvel(1,1,1))
          allocate(HYMAP_routing_struc(n)%fldfrc(1,1,1))
          allocate(HYMAP_routing_struc(n)%fldare(1,1,1))
          allocate(HYMAP_routing_struc(n)%sfcelv(1,1,1))

       enddo
    endif

    do n=1,LIS_rc%nnest
       call LIS_registerAlarm("HYMAP router model alarm",&
            HYMAP_routing_struc(n)%dt,HYMAP_routing_struc(n)%dt)
       
       call LIS_registerAlarm("HYMAP router output alarm",&
            HYMAP_routing_struc(n)%dt,HYMAP_routing_struc(n)%outInterval)
       
       call LIS_registerAlarm("HYMAP router restart alarm",&
            HYMAP_routing_struc(n)%dt,HYMAP_routing_struc(n)%rstInterval)          
    enddo
  end subroutine HYMAP_routingInit
  !=============================================
  !=============================================
  subroutine read_param_real(ctitle,z,n,array)
    
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
    real*4,       intent(inout) :: array(LIS_rc%gnc(n),LIS_rc%gnr(n),z)
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
            'nf90_open failed in HYMAP_routing_init')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in HYMAP_routing_init')
       
       call LIS_verify(nf90_get_var(ftn,varid, array), &
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in HYMAP_routing_init')
       
    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in HYMAP_routing_init'
       call LIS_endrun()
    endif
    
#endif
  end subroutine read_param_real

  subroutine read_param_int(ctitle,z,n,array)
    
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
    integer,      intent(inout) :: array(LIS_rc%gnc(n),LIS_rc%gnr(n),z)
    integer                     :: ftn 
    logical                     :: file_exists
    integer                     :: status
    integer                     :: c,r,l
    character*100               :: cfile
    integer                     :: varid
    real                        :: array1(LIS_rc%gnc(n),LIS_rc%gnr(n),z)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    inquire(file=LIS_rc%paramfile(n),exist=file_exists)
    
    if(file_exists) then 

       call LIS_verify(nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE, ncid = ftn), &
            'nf90_open failed in HYMAP_routing_init')
       call LIS_verify(nf90_inq_varid(ftn,trim(ctitle),varid), &
            'nf90_inq_varid failed for '//trim(ctitle)//&
            ' in HYMAP_routing_init')
       
       call LIS_verify(nf90_get_var(ftn,varid, array1), &
            'nf90_get_var failed for '//trim(ctitle)//&
            ' in HYMAP_routing_init')
       do l=1,z
          do r=1,LIS_rc%gnr(n)
             do c=1,LIS_rc%gnc(n)
                array(c,r,l) = nint(array1(c,r,l))
             enddo
          enddo
       enddo
    else
       write(LIS_logunit,*) '[ERR] parameter input file '//trim(LIS_rc%paramfile(n))
       write(LIS_logunit,*) '[ERR] failed in HYMAP_routing_init'
       call LIS_endrun()
    endif
    
#endif
  end subroutine read_param_int
end module HYMAP_routingMod
