!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module NLDAS_routingMod
!BOP
! 
! !MODULE: NLDAS_routingMod
! 
! !DESCRIPTION: 
!  This module provides the definition of data structures used to control
!  the operation of the NLDAS routing scheme. It also provides the entry
!  method for the initialization of NLDAS routing variables. 
! 
!  Reference: Lohmann et al. (2004), "Streamflow and water balance
!  intercomparisons for four land surface models in the North American Land
!  Data Assimilation System project", Journal of Geophysical Research, 
!  109, D07S91, doi:10.1029/2003JD003517. 
!
! !REVISION HISTORY: 
! 6 May 2011: Sujay Kumar, Initial implementation in LIS based on the 
!                          offline routing code obtained from NCEP. 
! 
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: NLDAS_routingInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: NLDAS_routing_struc
  
  type, public :: NLDAS_routing_dec
  
     real       :: dt 
     integer    :: luh !length of the internal unit-hydrograph in dt
     integer    :: duh !length of the internal unit-hydrograph in days
     integer    :: max_t !maximum travel time of water through one grid box 
                         !in seconds
     integer    :: ltr !length of transport unit-hydrograph in dt
     real, allocatable :: uh_intern(:,:,:)
     real, allocatable :: uh_trans(:,:,:)
     real, allocatable :: runoff_intern(:,:,:)
     real, allocatable :: runoff_trans(:,:,:)
     integer, allocatable :: order(:,:)
     integer          :: order_n

     real,  allocatable   :: streamflow(:,:)

     character*100    :: initial_1
     character*100    :: initial_2
     character(len=LIS_CONST_PATH_LEN) :: rstfile
     real, allocatable    :: area(:,:)

     integer          :: numout
     integer          :: fileopen
     real             :: outInterval 
     real             :: rstInterval
     character*20     :: startMode

     integer          :: mo
  end type NLDAS_routing_dec

  type(NLDAS_routing_dec), allocatable :: NLDAS_routing_struc(:)

contains
 
!BOP
!
! !ROUTINE: NLDAS_routingInit
! \label{NLDAS_routingInit}
! 
! !INTERFACE: 
  subroutine NLDAS_routingInit
! !USES: 
    use LIS_coreMod,   only : LIS_rc, LIS_config, LIS_masterproc
    use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, LIS_seconds2time,&
         LIS_parseTimeString, LIS_registerAlarm
    use LIS_logMod,    only : LIS_getNextUnitNumber, LIS_releaseUnitNumber, &
         LIS_endrun, LIS_logunit, LIS_verify
    use LIS_routingMod, only : LIS_runoff_state

    integer              :: n 
    integer              :: ftn 
    integer              :: status
    character*100        :: uh_file_1, uh_file_2
    character(len=LIS_CONST_PATH_LEN) :: order_file
    integer              :: i,j,k
    integer              :: yr, mo, da, hr, mn, ss
    type(ESMF_Grid)      :: global_grid
    type(ESMF_DistGrid)  :: global_grid_dg
    type(ESMF_ArraySpec) :: realarrspec
    type(ESMF_Field)     :: sf_runoff_field
    type(ESMF_Field)     :: baseflow_field
    type(ESMF_Field)     :: sf_runoff_count_field
    type(ESMF_Time)      :: currTime
    character*10         :: time
    real, pointer        :: sfrunoff(:)
    real, pointer        :: baseflow(:)
    real, pointer        :: runoff_count(:)
    integer              :: ios
    
!currently this is setup to run on a single processor
    allocate(NLDAS_routing_struc(LIS_rc%nnest))

    do n=1, LIS_rc%nnest
       NLDAS_routing_struc(n)%dt = 3600
       NLDAS_routing_struc(n)%duh = 2
       NLDAS_routing_struc(n)%luh = NLDAS_routing_struc(n)%duh*86400/&
            NLDAS_routing_struc(n)%dt
       NLDAS_routing_struc(n)%max_t = 8*3600
       NLDAS_routing_struc(n)%ltr = NLDAS_routing_struc(n)%max_t/&
            NLDAS_routing_struc(n)%dt
       if(LIS_masterproc) then 
          allocate(NLDAS_routing_struc(n)%streamflow(&
               LIS_rc%gnc(n),LIS_rc%gnr(n)))
       else
          allocate(NLDAS_routing_struc(n)%streamflow(1,1))
       endif
       NLDAS_routing_struc(n)%numout = 0 
       NLDAS_routing_struc(n)%fileopen = 0 
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "NLDAS routing model output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "NLDAS routing model output interval: not defined")
       
       call LIS_parseTimeString(time,NLDAS_routing_struc(n)%outInterval)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "NLDAS routing model restart interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "NLDAS routing model restart interval: not defined")

       call LIS_parseTimeString(time,NLDAS_routing_struc(n)%rstInterval)
    enddo

    if(LIS_masterproc) then 
       write(LIS_logunit,*) 'Initializing the NLDAS routing scheme....'
       
       do n=1, LIS_rc%nnest
          allocate(NLDAS_routing_struc(n)%order(4,LIS_rc%gnc(n)*LIS_rc%gnr(n)))
          allocate(NLDAS_routing_struc(n)%uh_intern(&
               NLDAS_routing_struc(n)%luh,LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(NLDAS_routing_struc(n)%uh_trans(&
               NLDAS_routing_struc(n)%ltr,LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(NLDAS_routing_struc(n)%runoff_intern(&
               NLDAS_routing_struc(n)%luh,LIS_rc%gnc(n),LIS_rc%gnr(n)))
          allocate(NLDAS_routing_struc(n)%runoff_trans(&
               NLDAS_routing_struc(n)%ltr,LIS_rc%gnc(n),LIS_rc%gnr(n)))
          
          allocate(NLDAS_routing_struc(n)%area(&
               LIS_rc%gnc(n),LIS_rc%gnr(n)))
       enddo
       
       call ESMF_ConfigFindLabel(LIS_config,&
            "NLDAS routing internal unit hydrograph file:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, uh_file_1,rc=status)
          call LIS_verify(status,&
               "NLDAS routing internal unit hydrograph file: not defined")
          
          ftn = LIS_getNextUnitNumber()
          open(ftn,file=trim(uh_file_1),status='old',form='unformatted',&
               iostat = ios)
          if(ios.ne.0) then 
             write(LIS_logunit,*) 'File '//trim(uh_file_1)//& 
                  'does not exist '
             call LIS_endrun()
          endif
          read(ftn) (((NLDAS_routing_struc(n)%uh_intern(i,j,k), &
               i=1,NLDAS_routing_struc(n)%luh), j=1,LIS_rc%gnc(n)),&
               k=1,LIS_rc%gnr(n))
          
          call LIS_releaseUnitNumber(ftn)
       enddo
       
       call ESMF_ConfigFindLabel(LIS_config,&
            "NLDAS routing transport unit hydrograph file:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, uh_file_2,rc=status)
          call LIS_verify(status,&
               "NLDAS routing transport unit hydrograph file: not defined")
          
          ftn = LIS_getNextUnitNumber()
          open(ftn,file=trim(uh_file_2),status='old',form='unformatted',&
               iostat = ios)
          if(ios.ne.0) then 
             write(LIS_logunit,*) 'File '//trim(uh_file_2)//& 
                  'does not exist '
             call LIS_endrun()
          endif
          read(ftn) (((NLDAS_routing_struc(n)%uh_trans(i,j,k), &
               i=1,NLDAS_routing_struc(n)%ltr), j=1,LIS_rc%gnc(n)),&
               k=1,LIS_rc%gnr(n))
          
          call LIS_releaseUnitNumber(ftn)
          
       enddo
       
       call ESMF_ConfigFindLabel(LIS_config,&
            "NLDAS routing coordinates order file:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, order_file,rc=status)
          call LIS_verify(status,&
               "NLDAS routing coordinates order file: not defined")

          ftn = LIS_getNextUnitNumber()
          open(ftn,file=trim(order_file),status='old',form='unformatted',&
               iostat = ios)
          if(ios.ne.0) then 
             write(LIS_logunit,*) 'File '//trim(order_file)//& 
                  'does not exist '
             call LIS_endrun()
          endif
          read(ftn) ((NLDAS_routing_struc(n)%order(i,j), &
               i=1,4),j=1,LIS_rc%gnc(n)*LIS_rc%gnr(n))          
          call LIS_releaseUnitNumber(ftn)

          do i=1,LIS_rc%gnc(n)*LIS_rc%gnr(n)
             if(NLDAS_routing_struc(n)%order(1,i).eq.0) then 
                NLDAS_routing_struc(n)%order_n = i-1
                exit
             endif
          enddo

          call calc_area(LIS_rc%gnc(n),LIS_rc%gnr(n),&
               NLDAS_routing_struc(n)%area, &
               nint(NLDAS_routing_struc(n)%dt))
       enddo
       

       call ESMF_ConfigFindLabel(LIS_config,&
            "NLDAS routing initial condition for runoff:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               NLDAS_routing_struc(n)%initial_1,rc=status)
          call LIS_verify(status,&
               "NLDAS routing initial condition for runoff: not defined")
       enddo

       call ESMF_ConfigFindLabel(LIS_config,&
            "NLDAS routing initial condition for transport:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               NLDAS_routing_struc(n)%initial_2,rc=status)
          call LIS_verify(status,&
               "NLDAS routing initial condition for transport: not defined")
       enddo

       call ESMF_ConfigFindLabel(LIS_config,&
            "NLDAS routing model start mode:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,&
               NLDAS_routing_struc(n)%startmode,rc=status)
          call LIS_verify(status,&
               "NLDAS routing model start mode: not defined")
       enddo

       call ESMF_ConfigFindLabel(LIS_config,&
            "NLDAS routing model restart file:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,&
               NLDAS_routing_struc(n)%rstfile,rc=status)
          call LIS_verify(status,&
               "NLDAS routing model restart file: not defined")
       enddo

       do n=1, LIS_rc%nnest
          call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
               rc=status)
          call LIS_verify(status, &
               "ESMF_ArraySpecSet failed in NLDAS_routingMod")
       
          global_grid_dg = ESMF_DistGridCreate(minIndex=(/1/),&
               maxIndex=(/LIS_rc%glbntiles(n)/),&
               regDecomp=(/1/),rc=status)
          call LIS_verify(status,&
               'ESMF_DistGridCreate failed in NLDAS_routingInit')
          
          global_grid = ESMF_GridCreate(name="Global Grid",&
               coordTypeKind=ESMF_TYPEKIND_R4,&
               distgrid = global_grid_dg, &
               gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
          call LIS_verify(status,&
               'ESMF_GridCreate failed in NLDAS_routingInit')

!create LSM interface objects to store runoff and baseflow

          sf_runoff_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=global_grid, name="Surface Runoff",rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')

          baseflow_field =ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=global_grid, name="Subsurface Runoff",rc=status)
          call LIS_verify(status, 'ESMF_FieldCreate failed')

          sf_runoff_count_field=ESMF_FieldCreate(arrayspec=realarrspec,&
               grid=global_grid, name="Surface Runoff count",rc=status)
          call LIS_verify(status, "ESMF_FieldCreate failed")

          call ESMF_FieldGet(sf_runoff_field,localDE=0,farrayPtr=sfrunoff,&
               rc=status)
          call LIS_verify(status, &
               "ESMF_FieldGet failed in NLDAS_routingMod")
          sfrunoff = 0.0 

          call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,&
               rc=status)
          call LIS_verify(status, &
               "ESMF_FieldGet failed in NLDAS_routingMod")
          baseflow = 0.0

          call ESMF_FieldGet(sf_runoff_count_field,localDE=0, &
               farrayPtr=runoff_count, rc=status)
          call LIS_verify(status,&
               "ESMF_FieldGet failed in NLDAS_routingMod")
          runoff_count = 0.0

          call ESMF_stateAdd(LIS_runoff_state(n),(/sf_runoff_field/),rc=status)
          call LIS_verify(status, 'ESMF_StateAdd failed for surface runoff')

          call ESMF_stateAdd(LIS_runoff_state(n),(/baseflow_field/),rc=status)
          call LIS_verify(status, 'ESMF_StateAdd failed for base flow')

          call ESMF_stateAdd(LIS_runoff_state(n),(/sf_runoff_count_field/),rc=status)
          call LIS_verify(status, 'ESMF_StateAdd failed for runoff count')

          NLDAS_routing_struc(n)%mo = -1

       enddo
    endif

    do n=1,LIS_rc%nnest
       call LIS_registerAlarm("NLDAS router model alarm",&
            NLDAS_routing_struc(n)%dt,NLDAS_routing_struc(n)%dt)
       
       call LIS_registerAlarm("NLDAS router output alarm",&
            NLDAS_routing_struc(n)%dt,NLDAS_routing_struc(n)%outInterval)
       
       call LIS_registerAlarm("NLDAS router restart alarm",&
            NLDAS_routing_struc(n)%dt,NLDAS_routing_struc(n)%rstInterval)
    enddo

  end subroutine NLDAS_routingInit


  SUBROUTINE CALC_AREA(NX, NY, AREA, DT)
    
    IMPLICIT NONE
    
    INTEGER NX, NY
    REAL    AREA(NX,NY)
    REAL    PI, DPHI, RERD, E, SOUTH, LAT
    INTEGER DT
    
    PARAMETER (RERD  = 6378.136)
    PARAMETER (DPHI  = 0.125)
    PARAMETER (E     = 0.00669447)
    PARAMETER (SOUTH = 25.0625)
    INTEGER I, J
    
    PI = ATAN(1.0) * 4.0
    
    DO J = 1, NY
       DO I = 1, NX
          LAT = SOUTH + (J-1) * DPHI
          AREA(I,J) = (2.0*PI*RERD*DPHI/360.0)**2.0 * &
               COS((LAT)*PI/180.0) !in km2
          
       END DO
    END DO
    
  END SUBROUTINE CALC_AREA
  
end module NLDAS_routingMod
