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
     integer          :: useens
! === Undefined Values ==================================
     integer          :: imis           !! real undefined value
! === River sequence ====================================
! === Map ===============================================
     real,    allocatable :: runoff0(:,:)     !! input runoff [mm.dt-1]
     real,    allocatable :: basflw0(:,:)     !! input baseflow [mm.dt-1]
 ! === Outputs ==========================================
     real,    allocatable :: rivout(:,:,:)      !! river outflow  [m3/s]

     character*100    :: rstfile
     integer          :: numout
     integer          :: fileopen
     real             :: outInterval 
     real             :: rstInterval
     character*20     :: startMode
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
    real, pointer        :: sfrunoff(:)
    real, pointer        :: baseflow(:)
    character*10         :: time

    allocate(RAPID_routing_struc(LIS_rc%nnest))
 
!TODO: change setting   
    do n=1, LIS_rc%nnest
!       RAPID_routing_struc(n)%rslpmin  = 1e-5  !! minimum slope
!       RAPID_routing_struc(n)%inz      = 10    !! number of stages in the sub-grid discretization
       RAPID_routing_struc(n)%imis     = -9999 !! undefined integer value
!       RAPID_routing_struc(n)%numout   = 0 
!       RAPID_routing_struc(n)%fileopen = 0 
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
     
    if(LIS_masterproc) then 
       write(LIS_logunit,*) '[INFO] Initializing RAPID....'
       !allocate matrixes
       !do n=1, LIS_rc%nnest
          !TODO: setting variables
       !enddo
                                  
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
       
!       if(LIS_rc%lsm.eq."none") then 
!          call initrunoffdata(trim(LIS_rc%runoffdatasource)//char(0))
!       endif

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
       enddo
    endif

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
