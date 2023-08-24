!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

module RAPID_routingMod
!BOP
! 
! !MODULE: RAPID_routingMod
! 
! !DESCRIPTION: 
!  This module provides the definition of data structures used to control
!  the operation of the RAPID routing scheme. It also provides the entry
!  method for the initialization of RAPID routing variables. 
! 
!  Reference: 
!
! !REVISION HISTORY: 
! 17 Mar 2021: Yeosang Yoon: Initial implementation in LIS based on the 
!                            RAPID offline routing code. 
! 25 Oct 2022: Yeosang Yoon: Support to run with LSM ensemble mean runoff variables
! 
! !USES: 
  use ESMF
 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: RAPID_routingInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: RAPID_routing_struc
  
  type, public :: RAPID_routing_dec
     
     real             :: dt
     integer          :: imis                   !! real undefined value

     character*100    :: rstfile
     integer          :: numout
     integer          :: fileopen
     real             :: outInterval 
     real             :: rstInterval
     real             :: routingInterval
     character*20     :: startmode

     logical          :: bQinit       ! initial flow
     logical          :: bQfinal      ! write final flow
     logical          :: bV           ! compute volume
     logical          :: bhum         ! human-induced flows
     logical          :: bfor         ! forcing
     logical          :: bdam         ! dam model used
     logical          :: binfluence   ! output influence
     logical          :: buq          ! uncertainty quantif.

     integer          :: run_opt      ! run option
     integer          :: routing_opt  ! routing option        
     integer          :: phi_opt      ! phi option
     character*200    :: connectfile  ! river connectivity file
     integer          :: max_reach    ! max number of upstream reaches
     character*200    :: weightfile   ! river weight table
     character*200    :: basinIDfile  ! river basin ID file
     character*200    :: kfile        ! Muskingum parameter k file
     character*200    :: xfile        ! Muskingum parameter x file

     character*200    :: nmlfile
 
     integer          :: n_riv_tot    ! number of river connectivity 
     integer          :: n_riv_bas    ! number of river basins  
     integer          :: n_wei_table  ! number of reach in weight table file
  
     logical             :: initCheck     ! for rapid_init

     integer,allocatable :: riv_bas_id(:) ! for rapid output
     real,allocatable    :: riv_tot_lon(:)
     real,allocatable    :: riv_tot_lat(:)

     real,allocatable    :: Qout(:)  ! instantaneous flow 
     integer             :: useens   ! ensemble mode
end type RAPID_routing_dec

  type(RAPID_routing_dec), allocatable :: RAPID_routing_struc(:)

contains
 
!BOP
!
! !ROUTINE: RAPID_routingInit
! \label{RAPID_routingInit}
! 
  subroutine RAPID_routingInit
    !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_routingMod    
        
    integer              :: n 
    integer              :: ftn 
    integer              :: status
    integer              :: i,j,k
    type(ESMF_ArraySpec) :: realarrspec
    type(ESMF_Field)     :: sf_runoff_field
    type(ESMF_Field)     :: baseflow_field
    type(ESMF_Field)     :: sf_runoff_count_field
    real, pointer        :: sfrunoff(:)
    real, pointer        :: baseflow(:)
    character*10         :: time

    allocate(RAPID_routing_struc(LIS_rc%nnest))
 
!TODO: change setting   
    do n=1, LIS_rc%nnest
       RAPID_routing_struc(n)%imis     = -9999 !! undefined integer value
       RAPID_routing_struc(n)%numout   = 0
       RAPID_routing_struc(n)%fileopen = 0
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID routing model time step:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "RAPID routing model time step: not defined")

       call LIS_parseTimeString(time,RAPID_routing_struc(n)%dt)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID routing model output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "RAPID routing model output interval: not defined")

       call LIS_parseTimeString(time,RAPID_routing_struc(n)%outInterval)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID river routing time step:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "RAPID river routing time step: not defined")

       call LIS_parseTimeString(time,RAPID_routing_struc(n)%routingInterval)
    enddo
    
    write(LIS_logunit,*) '[INFO] Initializing the RAPID routing scheme....'

    do n=1, LIS_rc%nnest
       RAPID_routing_struc(n)%initCheck = .true.
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID routing model start mode:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%startmode,rc=status)
       call LIS_verify(status,&
            "RAPID routing model start mode: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID routing model restart interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "RAPID routing model restart interval: not defined")

       call LIS_parseTimeString(time,RAPID_routing_struc(n)%rstInterval)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID routing model restart file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%rstfile,rc=status)
       call LIS_verify(status,&
            "RAPID routing model restart file: not defined")
    enddo

    ! set high-level options that govern how the model is to run
    call ESMF_ConfigFindLabel(LIS_config,"RAPID initial flow:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%bQinit,rc=status)
       call LIS_verify(status,"RAPID initial flow: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID write final flow:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%bQfinal,rc=status)
       call LIS_verify(status,"RAPID write final flow: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID compute volume:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%bV,rc=status)
       call LIS_verify(status,"RAPID compute volume: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID human-induced flow:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%bhum,rc=status)
       call LIS_verify(status,"RAPID human-induced flow: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID upstream forcing:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%bfor,rc=status)
       call LIS_verify(status,"RAPID upstream forcing: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID dam model used:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%bdam,rc=status)
       call LIS_verify(status,"RAPID dam model used: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID output influence:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%binfluence,rc=status)
       call LIS_verify(status,"RAPID output influence: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID uncertainty quantification:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%buq,rc=status)
       call LIS_verify(status,"RAPID uncertainty quantification: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID run option:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%run_opt,rc=status)
       call LIS_verify(status,"RAPID run option: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID routing option:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%routing_opt,rc=status)
       call LIS_verify(status,"RAPID routing option: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID cost function phi option:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%phi_opt,rc=status)
       call LIS_verify(status,"RAPID cost function phi option: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID river connectivity file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%connectfile,rc=status)
       call LIS_verify(status,"RAPID river connectivity file: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID max number of upstream reaches:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%max_reach,rc=status)
       call LIS_verify(status,"RAPID max number of upstream reaches: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID river weight table:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%weightfile,rc=status)
       call LIS_verify(status,"RAPID river weight table: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID river basin ID file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%basinIDfile,rc=status)
       call LIS_verify(status,"RAPID river basin ID file: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID Muskingum parameter k file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%kfile,rc=status)
       call LIS_verify(status,"RAPID Muskingum parameter k file: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"RAPID Muskingum parameter x file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%xfile,rc=status)
       call LIS_verify(status,"RAPID Muskingum parameter x file: not defined")
    enddo

    ! namelist
    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID namelist file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%nmlfile,rc=status)
       call LIS_verify(status,&
            "RAPID namelist file: not defined")
    enddo

    ! ensemble mode
    call ESMF_ConfigFindLabel(LIS_config,&
         "RAPID run in ensemble mode:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            RAPID_routing_struc(n)%useens,rc=status)
       call LIS_verify(status,&
            "RAPID run in ensemble mode: not defined")
    enddo

    ! ensemble mode 0: open loop, ensemble mode 1: ensemble mean 
    ! EMK...Add do loop
    do n=1, LIS_rc%nnest
       if(RAPID_routing_struc(n)%useens>1) then
          write(LIS_logunit,*) "[ERR] Currently RAPID only supports ensemble modes 0 or 1" 
          call LIS_endrun()
       endif
    end do

    ! checks the size of static data for RAPID
    do n=1, LIS_rc%nnest
       call RAPID_check_domain_size(n)
    enddo

    ! for RAPID output
    do n=1, LIS_rc%nnest
       allocate(RAPID_routing_struc(n)%riv_bas_id(RAPID_routing_struc(n)%n_riv_bas))
       RAPID_routing_struc(n)%riv_bas_id=-9999

       allocate(RAPID_routing_struc(n)%riv_tot_lon(RAPID_routing_struc(n)%n_riv_bas))
       RAPID_routing_struc(n)%riv_tot_lon=-9999
       allocate(RAPID_routing_struc(n)%riv_tot_lat(RAPID_routing_struc(n)%n_riv_bas))
       RAPID_routing_struc(n)%riv_tot_lat=-9999
    enddo

    ! for RAPID restart
    do n=1, LIS_rc%nnest
       allocate(RAPID_routing_struc(n)%Qout(RAPID_routing_struc(n)%n_riv_bas))
       RAPID_routing_struc(n)%Qout=0
    enddo
     
    do n=1, LIS_rc%nnest
       call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status,&
            "ESMF_ArraySpecSet failed in RAPID_routingMod")

       !create LSM interface objects to store runoff and baseflow
       sf_runoff_field =ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_vecTile(n), name="Surface Runoff",rc=status)
       call LIS_verify(status, 'ESMF_FieldCreate failed')

       baseflow_field =ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_vecTile(n), name="Subsurface Runoff",rc=status)
       call LIS_verify(status, 'ESMF_FieldCreate failed')

       call ESMF_FieldGet(sf_runoff_field,localDE=0,farrayPtr=sfrunoff,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_FieldGet failed in RAPID_routingMod")
       sfrunoff = 0.0

       call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_FieldGet failed in RAPID_routingMod")
       baseflow = 0.0

       call ESMF_stateAdd(LIS_runoff_state(n),(/sf_runoff_field/),rc=status)
       call LIS_verify(status, 'ESMF_StateAdd failed for surface runoff')

       call ESMF_stateAdd(LIS_runoff_state(n),(/baseflow_field/),rc=status)
       call LIS_verify(status, 'ESMF_StateAdd failed for base flow')
    enddo
   
    do n=1,LIS_rc%nnest
       call LIS_registerAlarm("RAPID router model alarm",&
            RAPID_routing_struc(n)%dt,RAPID_routing_struc(n)%dt)
       
       call LIS_registerAlarm("RAPID router output alarm",&
            RAPID_routing_struc(n)%dt,RAPID_routing_struc(n)%outInterval)
       
       call LIS_registerAlarm("RAPID router restart alarm",&
            RAPID_routing_struc(n)%dt,RAPID_routing_struc(n)%rstInterval)          
    enddo
   
  end subroutine RAPID_routingInit

end module RAPID_routingMod
