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
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_lsmMod
!BOP
!
! !MODULE: LIS_lsmMod
! 
! !DESCRIPTION:
!  The code in this file provides interfaces to manage the operation
!  of different land surface models
!
!  \subsubsection{Overview}
!   This module defines the interface plugins for the incorporation of
!   different land surface models. These interfaces provide entry points
!   for introducing a new land surface scheme in LIS. The following 
!   implementations for each LSM are expected to specify 
!   methods to initialize and set the LSM specific variables and parameters, 
!   provide methods for model simulation, model output, and restart 
!   operations. A number of other optional interfaces need to specified
!   depending on the mode of operation of the LSM. For example, if the 
!   LSM is used for a coupled simulation with an atmospheric component, 
!   the {\tt LIS\_lsm\_setexport} interface need to be implemented. Similar
!   implementations need to be specified for the use of the LSM in 
!   data assimilation and optimization applications. 
! 
! !REVISION HISTORY: 
!  14 Nov 2002    Sujay Kumar  Initial Specification
! 
  use ESMF
  use LIS_coreMod
  use LIS_perturbMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use LIS_DAobs_pluginMod
  use LIS_routingMod
  use map_utils

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_lsm_init       ! initialize lsm variables, memory
  public :: LIS_setuplsm       ! set land surface parameters
  public :: LIS_lsm_readrestart    ! read the restart file
  public :: LIS_lsm_f2t     ! transfer forcing to model tiles
  public :: LIS_lsm_run       ! execute the land model 
  public :: LIS_setLSMDynParams    ! set the time dependent parameters
  public :: LIS_lsm_writerestart   ! write the restart file
  public :: LIS_lsm_setexport  ! set the variables to be exported to the
                               ! atmospheric component
  public :: LIS_lsm_perturb_states ! perturbs the prognostic variables
  public :: LIS_lsm_finalize   ! cleanup allocated structures
  public :: LIS_lsm_reset      ! reset structures
  public :: LIS_lsm_diagnoseVarsForDA  ! DA related updates to variables LSM
  public :: LIS_lsm_DAGetObsPred
  public :: LIS_lsm_DAGetStateVar
  public :: LIS_lsm_DASetStateVar
  public :: LIS_lsm_DAScaleStateVar
  public :: LIS_lsm_DADescaleStateVar
  public :: LIS_lsm_DAUpdateState
  public :: LIS_lsm_DAQCState
  public :: LIS_lsm_DAgetStateSpaceSize
  public :: LIS_lsm_DAextractStateVector
  public :: LIS_lsm_DAsetFreshIncrementsStatus
  public :: LIS_lsm_DAgetFreshIncrementsStatus
  public :: LIS_lsm_DAsetAnlysisUpdates
  public :: LIS_lsm_DAmapTileSpaceToObsSpace
  public :: LIS_lsm_DAgetStateVarNames
  public :: LIS_lsm_DAobsTransform
  public :: LIS_lsm_DAmapObsToLSM
  public :: LIS_lsm_DAqcObsState
  public :: LIS_lsm_getlatlons
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------  
  public :: LIS_LSM_State
  public :: LIS_LSM_Pert_State
  public :: LIS_LSM_Incr_State

  public :: LIS_LSM2SUBLSM_State
  public :: LIS_SUBLSM2LSM_State
!EOP

  type(ESMF_State), allocatable :: LIS_LSM_State(:,:) !ESMF state of prognostic
                                                 !variables
  type(ESMF_State), allocatable :: LIS_LSM_Pert_State(:,:)!ESMF state of prognostic
                                                 !variable perturbations
  
  type(ESMF_State), allocatable :: LIS_LSM_Incr_State(:,:)!ESMF State of LSM state
                                                 !increments


  type(ESMF_State), allocatable :: LIS_LSM2SUBLSM_State(:,:)
  type(ESMF_State), allocatable :: LIS_SUBLSM2LSM_State(:,:)

!BOP
! !ROUTINE: LIS_lsm_setexport
! \label{LIS_lsm_setexport}
!
! !INTERFACE: 
  interface LIS_lsm_setexport
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure lsm_setexport_noesmf
! !DESCRIPTION: 
! This interface provides the entry point for specifying an 
! export state (a list of model specific variables) from a land surface
! model. The routine is used in a coupled simulation to provide feedback
! to a different model component such as an atmospheric model. 
!EOP
  end interface

contains

!BOP
! 
! !ROUTINE: LIS_lsm_init
! \label{LIS_lsm_init}
! 
! !INTERFACE:
  subroutine LIS_lsm_init()
! !USES:
    use LIS_surfaceModelDataMod
!
! !DESCRIPTION:
! Interface for initializing the land model. The intialization includes
! the allocation of memory for LSM specific variables and datastructures
! and specification of any runtime specific options. 
!
! The calling sequence is: 
! \begin{description}
!  \item[lsminit] (\ref{lsminit}) \newline
!    invokes the generic method in the registry to initialize the 
!    land surface model
!  \item[perturbinit](\ref{perturbinit}) \newline
!    invokes the generic method in the registry to perturb the 
!    prognostic variables using the specified algorithm.
! \end{description}
!EOP
    integer       :: n 
    integer       :: status
    character*1   :: nestid(2)
    character*1   :: caseid(3)
    character*100 :: temp
    integer       :: ftn

    integer              :: i,j,k
    type(pert_dec_type)  :: lsm_pert
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
    integer              :: max_index
    logical              :: name_found
    character*20         :: alglist(10)
    logical              :: LSM_DAvalid
    integer              :: rc

    TRACE_ENTER("lsm_init")
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%lsm,&
         label="Land surface model:",rc=rc)
    call LIS_verify(rc,'Land surface model: option not specified in the config file')
    

    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nSubLSMs,&
         label="Number of subLSMs:",default=0,rc=rc)
    call LIS_verify(rc,'Number of subLSMs: option not specified in the config file')

    allocate(LIS_rc%subLSM(LIS_rc%nSubLSMs))

    if(LIS_rc%nSubLSMs.ge.1) then        
       
       call ESMF_ConfigFindLabel(LIS_config,"subLSM models:",rc=rc)
       do i=1,LIS_rc%nSubLSMs
          call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%subLSM(i),rc=rc)
          call LIS_verify(rc,"subLSM models: not defined")
       enddo
       
       allocate(LIS_LSM2SUBLSM_State(LIS_rc%nnest, LIS_rc%nsubLSMs))
       allocate(LIS_SUBLSM2LSM_State(LIS_rc%nnest, LIS_rc%nsubLSMs))

       do n=1,LIS_rc%nnest
          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid
          
          do i=1,LIS_rc%nSubLSMs
             
             write(unit=temp,fmt='(i3.3)') i
             read(unit=temp,fmt='(3a1)') caseid
             
             LIS_LSM2SUBLSM_State(n,i) = ESMF_StateCreate(&
                  name="LSM2SUBLSM State"//&
                  nestid(1)//nestid(2)&
                  //'_'//caseid(1)//caseid(2)//caseid(3), rc=status)
             call LIS_verify(status, &
                  "ESMF_StateCreate failed for LSM2SUBLSM in LIS_lsm_init")
             
             LIS_SUBLSM2LSM_State(n,i) = ESMF_StateCreate(&
                  name="SUBLSM2LSM State"//&
                  nestid(1)//nestid(2)&
                  //'_'//caseid(1)//caseid(2)//caseid(3), rc=status)
             call LIS_verify(status, &
                  "ESMF_StateCreate failed for SUBLSM2LSM in LIS_lsm_init")
             
          enddo
       enddo
    endif
    
    allocate(LIS_rc%LSM_DAinst_valid(LIS_rc%ndas))

    do i=1,LIS_rc%ndas
       call LIS_isDAinstanceValid(LIS_rc%daset(i),&
            "LSM",LIS_rc%LSM_DAinst_valid(i))
    enddo
    
!---------------------------------------------------------------------
!  create the LSM state if data assimilation is being done, 
!  create the LSM perturbation state only if perturbation
!  option is turned on. 
!---------------------------------------------------------------------
    LSM_DAvalid = .false. 

    if(LIS_rc%ndas.gt.0.or.LIS_rc%nperts.gt.0) then 

       LSM_DAvalid = .true. 

       do i=1,LIS_rc%ndas
          LSM_DAvalid = LSM_DAvalid.and.LIS_rc%LSM_DAinst_valid(i)
       enddo

       if(LSM_DAvalid) then 
          allocate(LIS_LSM_State(LIS_rc%nnest, LIS_rc%nperts))
          allocate(LIS_LSM_Incr_State(LIS_rc%nnest, LIS_rc%nperts))
          
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
                     "ESMF_ArraySpecSet failed in LIS_lsm_init")
                
                write(unit=temp,fmt='(i2.2)') n
                read(unit=temp,fmt='(2a1)') nestid
                
                write(unit=temp,fmt='(i3.3)') k
                read(unit=temp,fmt='(3a1)') caseid
                
                LIS_LSM_State(n,k) = ESMF_StateCreate(name="LSM State"//&
                     nestid(1)//nestid(2)&
                     //'_'//caseid(1)//caseid(2)//caseid(3), rc=status)
                call LIS_verify(status, &
                     "ESMF_StateCreate failed in LIS_lsm_init")
                
                LIS_LSM_Incr_State(n,k) = ESMF_StateCreate(name="LSM Incr State"//&
                     nestid(1)//nestid(2)// &
                     '_'//caseid(1)//caseid(2)//caseid(3), rc=status)
                call LIS_verify(status,&
                     "ESMF_StateCreate failed in LIS_lsm_init")
                
                do i=1,LIS_rc%nstvars(k)
                   read(ftn,fmt='(a40)') vname(i)
                   read(ftn,*) stmin(i),stmax(i)
                   write(LIS_logunit,*) '[INFO] ',vname(i),stmin(i),stmax(i)
                   
                   varField = ESMF_FieldCreate(grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
                        arrayspec=arrspec1,name=trim(vname(i)), rc=status)
                   call LIS_verify(status, &
                        "ESMF_FieldCreate failed in LIS_lsm_init")
                   
                   varIncrField = ESMF_FieldCreate(grid=&
                        LIS_vecPatch(n,LIS_rc%lsm_index),&
                        arrayspec=arrspec1,name=trim(vname(i)), rc=status)
                   call LIS_verify(status,&
                        "ESMF_FieldCreate failed in LIS_lsm_init")
                   
                   call ESMF_AttributeSet(varField,"Max Value",stmax(i),rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttribteSet failed in LIS_lsm_init")
                   
                   call ESMF_AttributeSet(varField,"Min Value",stmin(i),rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet failed in LIS_lsm_init")
                   
                   
                   call ESMF_AttributeSet(VarIncrField,"Max Value",stmax(i),rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet failed in LIS_lsm_init")
                   
                   call ESMF_AttributeSet(VarIncrField,"Min Value",stmin(i),rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet failed in LIS_lsm_init")
                   
                   call ESMF_StateAdd(LIS_LSM_State(n,k),(/varField/),rc=status)
                   call LIS_verify(status,&
                        "ESMF_StateAdd failed in LIS_lsm_init")
                   
                   call ESMF_StateAdd(LIS_LSM_Incr_State(n,k), &
                        (/VarIncrField/), rc=status)
                   call LIS_verify(status,&
                        "ESMF_StateAdd failed in LIS_lsm_init")
!----------------------------------------------------------------------------
! Initially set the fresh increments available status to false. 
!----------------------------------------------------------------------------
                   call ESMF_AttributeSet(LIS_LSM_Incr_State(n,k), &
                        name="Fresh Increments Status", value=.false., &
                        rc=status)
                   call LIS_verify(status,&
                        "ESMF_AttributeSet failed in LIS_lsm_init")
                enddo
                deallocate(vname)
                deallocate(stmin)
                deallocate(stmax)
                call LIS_releaseUnitNumber(ftn)
             enddo
             LIS_sfmodel_struc(n)%models_used = &
                  trim(LIS_sfmodel_struc(n)%models_used)//&
                  trim(LIS_rc%lsm)
          enddo
       
          if(LIS_rc%nperts.gt.0) then 
             allocate(LIS_LSM_Pert_State(LIS_rc%nnest, LIS_rc%nperts))
             
             call ESMF_ArraySpecSet(arrspec2,rank=1,typekind=ESMF_TYPEKIND_R4,&
                  rc=status)
             call LIS_verify(status,&
                  "ESMF_ArraySpecSet failed in LIS_lsm_init")
             
             do n=1,LIS_rc%nnest
                allocate(ssdev(LIS_rc%ngrid(n)))
                do k=1,LIS_rc%nperts
                   if(LIS_rc%perturb_state(k).ne."none") then 
                      allocate(lsm_pert%vname(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%perttype(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%ssdev(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%stdmax(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%zeromean(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%tcorr(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%xcorr(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%ycorr(LIS_rc%nstvars(k)))
                      allocate(lsm_pert%ccorr(LIS_rc%nstvars(k),LIS_rc%nstvars(k)))
                      
                      write(unit=temp,fmt='(i2.2)') n
                      read(unit=temp,fmt='(2a1)') nestid
                      
                      LIS_LSM_Pert_State(n,k) = ESMF_StateCreate(&
                           name="LSM_Pert_State"//&
                           nestid(1)//nestid(2),&
                           rc=status)
                      call LIS_verify(status,&
                           "ESMF_StateCreate: LSM_Pert_State failed in LIS_lsm_init")
                      
                      call LIS_readPertAttributes(LIS_rc%nstvars(k),&
                           LIS_rc%progpertAttribfile(k),&
                           lsm_pert)
                      
                      do i=1,LIS_rc%nstvars(k)
                         pertField = ESMF_FieldCreate(grid=&
                              LIS_vecPatch(n,LIS_rc%lsm_index),&
                              arrayspec=arrspec2,name=trim(lsm_pert%vname(i)),&
                              rc=status)
                         
                         call ESMF_StateAdd(LIS_LSM_Pert_State(n,k),(/pertField/),&
                              rc=status)
                         call LIS_verify(status,&
                              "ESMF_StateAdd failed in LIS_lsm_init")
                      enddo
                      
                      allocate(pertobjs(LIS_rc%nstvars(k)))
                      allocate(order(LIS_rc%nstvars(k)))
                      allocate(ccorr(LIS_rc%nstvars(k),LIS_rc%nstvars(k)))
                      order = -1
                      
                      call ESMF_StateGet(LIS_LSM_Pert_State(n,k),&
                           itemNameList=pertobjs,rc=status)
                      call LIS_verify(status,&
                           "ESMF_StateGet failed in LIS_lsm_init")
                      
                      do i=1,LIS_rc%nstvars(k)
                         do j=1,LIS_rc%nstvars(k)
                            if(lsm_pert%vname(j).eq.pertobjs(i)) then 
                               order(i) = j
                               exit;
                            endif
                         enddo
                      enddo
                      
                      do i=1,LIS_rc%nstvars(k)
                         do j=1,LIS_rc%nstvars(k)
                            ccorr(i,j) = lsm_pert%ccorr(order(i),order(j))
                         enddo
                      enddo
                      
                      do i=1,LIS_rc%nstvars(k)
                         call ESMF_StateGet(LIS_LSM_Pert_State(n,k),&
                              pertobjs(i),pertField,rc=status)
                         call LIS_verify(status,&
                              "ESMF_StateGet failed in LIS_lsm_init")
                         
                         call ESMF_AttributeSet(pertField,"Perturbation Type",&
                              lsm_pert%perttype(order(i)),&
                              rc=status)
                         call LIS_verify(status,&
                              "ESMF_AttributeSet: Perturbation Type failed in LIS_lsm_init")
                         
                         if(LIS_rc%ngrid(n).gt.0) then 
                            ssdev = lsm_pert%ssdev(order(i))
                            
                            call ESMF_AttributeSet(pertField,"Standard Deviation",&
                                 ssdev,itemCount=LIS_rc%ngrid(n),&
                                 rc=status)
                            call LIS_verify(status,&
                                 "ESMF_AttributeSet: Standard Deviation failed in LIS_lsm_init")
                         endif
                         call ESMF_AttributeSet(pertField,"Std Normal Max",&
                              lsm_pert%stdmax(order(i)),&
                              rc=status)
                         call LIS_verify(status,&
                              "ESMF_AttributeSet: Std Normal Max failed in LIS_lsm_init")
                         
                         call ESMF_AttributeSet(pertField,"Ensure Zero Mean",&
                              lsm_pert%zeromean(order(i)),&
                              rc=status)
                         call LIS_verify(status,&
                              "ESMF_AttributeSet: Ensure Zero Mean failed in LIS_lsm_init")
                         
                         call ESMF_AttributeSet(pertField,&
                              "Temporal Correlation Scale",&
                              lsm_pert%tcorr(order(i)), rc=status)
                         call LIS_verify(status,&
                              "ESMF_AttributeSet: Temporal Correlation Scale failed in LIS_lsm_init")
                         
                         call ESMF_AttributeSet(pertField,"X Correlation Scale",&
                              lsm_pert%xcorr(order(i)), rc=status)
                         call LIS_verify(status,&
                              "ESMF_AttributeSet: X Correlation Scale failed in LIS_lsm_init")
                         
                         call ESMF_AttributeSet(pertField,"Y Correlation Scale",&
                              lsm_pert%ycorr(order(i)), rc=status)
                         call LIS_verify(status,&
                              "ESMF_AttributeSet: Y Correlation Scale failed in LIS_lsm_init")
                         
                         call ESMF_AttributeSet(pertField,&
                              "Cross Correlation Strength",&
                              ccorr(i,:), itemCount=LIS_rc%nstvars(k),&
                              rc=status)
                         call LIS_verify(status,&
                              "ESMF_AttributeSet: Cross Correlation Strength failed in LIS_lsm_init")
                         
                      enddo
                      deallocate(pertobjs)
                      deallocate(order)
                      deallocate(ccorr)
                      deallocate(lsm_pert%vname)
                      deallocate(lsm_pert%perttype)
                      deallocate(lsm_pert%ssdev)
                      deallocate(lsm_pert%stdmax)
                      deallocate(lsm_pert%zeromean)
                      deallocate(lsm_pert%tcorr)
                      deallocate(lsm_pert%xcorr)
                      deallocate(lsm_pert%ycorr)
                      deallocate(lsm_pert%ccorr)
                   endif
                enddo
                deallocate(ssdev)          
             enddo
          endif
       endif
    end if

    call lsminit(trim(LIS_rc%lsm)//char(0))
    
    if(LSM_DAvalid) then 
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
             call perturbinit(trim(alglist(i))//char(0), 2)
          enddo
       endif
       
       do i=1,LIS_rc%nperts   
          if(LIS_rc%perturb_state(i).ne."none") then 
             call perturbsetup(trim(LIS_rc%perturb_state(i))//char(0), 2, i, &
                  LIS_LSM_State(:,i), LIS_LSM_Pert_State(:,i))
          endif
       enddo
    endif

    do i=1,LIS_rc%nSubLSMs
       call sublsminit(trim(LIS_rc%subLSM(i))//char(0),i)
    enddo

    TRACE_EXIT("lsm_init")

  end subroutine LIS_lsm_init

!BOP
! !ROUTINE: LIS_setuplsm
! \label{LIS_setuplsm}
!
! !INTERFACE:
  subroutine  LIS_setuplsm
! !USES:
!
! !DESCRIPTION:
! The setup interfaces are used for the LSM specification of
! LSM bparameters. If a parameter falls outside the LIS-specified
! generic parameter list, the LSM is expected to provide routines for the 
! handling of any input data, specific to those parameters.
! 
! The calling sequence is: 
! \begin{description}
!  \item[lsmsetup] (\ref{lsmsetup}) \newline
!    invokes the generic method in the registry to set up the 
!    land surface model
! \end{description}
!EOP
    
    integer              :: i 

    TRACE_ENTER("lsm_setup")
    call lsmsetup(trim(LIS_rc%lsm)//char(0))

    do i=1,LIS_rc%nSubLSMs
       call sublsmsetup(trim(LIS_rc%subLSM(i))//char(0))
    enddo

    TRACE_EXIT("lsm_setup")
  end subroutine LIS_setuplsm

!BOP
! !ROUTINE: LIS_lsm_run
! \label{LIS_lsm_run}
!
! !INTERFACE:            
  subroutine LIS_lsm_run(n)
! !USES:

! !ARGUMENTS: 
    integer, intent(in) :: n 
!
! !DESCRIPTION:
! This interface provides the entry point to the LSM routines that
! invokes the land surface model physics. 
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmrun] (\ref{lsmrun}) \newline
!    invokes the generic method in the registry to run the
!    land surface model
! \end{description}
!EOP
    integer              :: i 

    TRACE_ENTER("lsm_run")
    call lsmrun(trim(LIS_rc%lsm)//char(0), n)

    do i=1,LIS_rc%nSubLSMs    
       if(LIS_rc%lsm.ne."none") then 
          
          !Get the export states from the LSM
          call lsm2sublsmgetexport(trim(LIS_rc%lsm)//"+"&
               //trim(LIS_rc%subLSM(i))//char(0), n,&
                LIS_LSM2SUBLSM_State(n,i))
       
          !Assign the LSM export states within the subLSM
          call sublsmsetlsmimport(&
               trim(LIS_rc%subLSM(i))//char(0), n,&
               LIS_LSM2SUBLSM_State(n,i))
       endif

       !Run the subLSM
       call sublsmrun(trim(LIS_rc%subLSM(i))//char(0), n)

       if(LIS_rc%lsm.ne."none") then 
          !Get the export states from the subLSM
          call sublsm2lsmgetexport(trim(LIS_rc%lsm)//"+"&
               //trim(LIS_rc%subLSM(i))//char(0), n,&
               LIS_SUBLSM2LSM_State(n,i))
          
          !Assign the subLSM export states within the LSM
          call lsmsetsublsmimport(trim(LIS_rc%lsm)//char(0), n,&
               LIS_SUBLSM2LSM_State(n,i))
       endif
    enddo


    TRACE_EXIT("lsm_run")

  end subroutine LIS_lsm_run

!BOP
! 
! !ROUTINE: LIS_lsm_perturb_states
! \label{LIS_lsm_perturb_states}
! 
! !INTERFACE: 
  subroutine LIS_lsm_perturb_states(n)
! !USES: 

! !ARGUMENTS: 
    integer,    intent(IN)  :: n 
!
! !DESCRIPTION:
! This interface provides the entry point to the LSM routines that 
! computes the perturbations on the LSM prognostic state variables
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[perturbmethod] (\ref{perturbmethod}) \newline
!    invokes the abstract method to invoke the perturbation
!    algorithm to perturb LSM prognostic variables
!  \item[lsmdagetstatevar] (\ref{lsmdagetstatevar}) \newline
!    obtains the list of prognostic variables
!  \item[applyLSMPert] (\ref{applyLSMPert}) \newline
!    applies the specified perturbations to the LSM state
!  \item[lsmdaqcstate] (\ref{lsmdaqcstate}) \newline
!   performs the QC of the perturbed LSM state
!  \item[lsmdasetstatevar] (\ref{lsmdasetstatevar}) \newline
!   assigns the prognostic variables back to the model states
! \end{description}
!EOP

    real                    :: curr_time
    integer                 :: k 

    TRACE_ENTER("lsm_perturb")
    do k=1, LIS_rc%nperts
       if(LIS_rc%LSM_DAinst_valid(k)) then
          if(LIS_rc%perturb_state(k).ne."none") then 
             curr_time = float(LIS_rc%hr)*3600+60*float(LIS_rc%mn)+float(LIS_rc%ss)
             if(mod(curr_time,real(LIS_rc%pertstateInterval(k))).eq.0) then
!------------------------------------------------------------------------
!   Returns the perturbed state based on the chosen algorithm
!------------------------------------------------------------------------
                call perturbmethod(trim(LIS_rc%perturb_state(k))//char(0),2, n,k,&
                     LIS_LSM_State(n,k),LIS_LSM_Pert_State(n,k))
!------------------------------------------------------------------------
!   Propagate step or applying the perturbations to prognostic states
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!   apply the lsm perturbations to the the LSM state
!------------------------------------------------------------------------
                call lsmdagetstatevar(trim(LIS_rc%lsm)//"+"//&
                     trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k))
                call applyLSMPert(n, k, LIS_LSM_State(n,k),LIS_LSM_Pert_State(n,k))
!------------------------------------------------------------------------
!   Diagnose the perturbed state (updates the model prognostic states)
!------------------------------------------------------------------------
                call lsmdaqcstate(trim(LIS_rc%lsm)//"+"//&
                     trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k))
                call lsmdasetstatevar(trim(LIS_rc%lsm)//"+"//&
                     trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k)) 
             endif
          endif
       endif
    enddo
    TRACE_EXIT("lsm_perturb")
  end subroutine LIS_lsm_perturb_states
!BOP
! !ROUTINE: LIS_lsm_readrestart
! \label{LIS_lsm_readrestart}
!
! !INTERFACE:
  subroutine LIS_lsm_readrestart
! !USES:

!
! !DESCRIPTION:
!  This interface provides the entry point to read LSM specific 
!  model restart files. 
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmrestart] (\ref{lsmrestart}) \newline
!    invokes the generic method in the registry to read
!    restart files for the land surface model 
! \end{description}
!EOP
    integer              :: i 

    TRACE_ENTER("lsm_readrst")
    call lsmrestart(trim(LIS_rc%lsm)//char(0))

    do i=1,LIS_rc%nSubLSMs
       call sublsmrestart(trim(LIS_rc%subLSM(i))//char(0))
    enddo
    TRACE_EXIT("lsm_readrst")
  end subroutine LIS_lsm_readrestart


!BOP
! !ROUTINE: LIS_setLSMDynparams
! \label{LIS_setLSMDynparams}
!
! !INTERFACE:
  subroutine LIS_setLSMDynparams(n)
! !USES:

! !ARGUMENTS: 
    integer, intent(in) :: n

! 
! !DESCRIPTION:
! This interface provides an entry point to update any 
! time dependent land surface model parameters in an LSM. 
! 
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmdynsetup] (\ref{lsmdynsetup}) \newline
!    invokes the generic method in the registry to set up 
!    time dependent parameters for the land surface model 
! \end{description}
!EOP
    integer              :: i 

    TRACE_ENTER("lsm_dynsetup")
    call lsmdynsetup(trim(LIS_rc%lsm)//char(0),n)

    do i=1,LIS_rc%nSubLSMs
       call sublsmdynsetup(trim(LIS_rc%subLSM(i))//char(0),n)
    enddo
    TRACE_EXIT("lsm_dynsetup")
  end subroutine LIS_setLSMDynparams

!BOP
! !ROUTINE: LIS_lsm_f2t
! \label{LIS_lsm_f2t}
!
! !INTERFACE:
  subroutine LIS_lsm_f2t(n)
! !USES: 

! !ARGUMENTS: 
    integer, intent(in) :: n
! 
! !DESCRIPTION:
!  
!  This interface is used to transfer the forcing variables to the actual 
!  model tile space. Any forcing perturbations that need to be applied to 
!  the input forcing is applied at this stage as well. 
!
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmf2t] (\ref{lsmf2t}) \newline
!    invokes the generic method in the registry to transfer the
!    forcing to the land surface model tiles
! \end{description}
!EOP
    integer              :: i 

    TRACE_ENTER("lsm_f2t")
    call lsmf2t(trim(LIS_rc%lsm)//"+"//trim(LIS_rc%runmode)//char(0),&
         n)

    do i=1,LIS_rc%nSubLSMs
       call sublsmf2t(trim(LIS_rc%subLSM(i))//"+"//&
            trim(LIS_rc%runmode)//char(0), n)
    enddo
    TRACE_EXIT("lsm_f2t")

  end subroutine LIS_lsm_f2t


!BOP
! !ROUTINE: LIS_lsm_writerestart
! \label{LIS_lsm_writerestart}
! 
! !INTERFACE:
  subroutine LIS_lsm_writerestart(n)
! !USES:

! !ARGUMENTS: 
    integer, intent(in) :: n 
! !DESCRIPTION:
!  This interface provides the entry point to read the LSM specific 
!  model restart files. 
!
! The arguments are: 
! \begin{description}
! \item[n]
!  index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmwrst] (\ref{lsmwrst}) \newline
!    invokes the generic method in the registry to write
!    restart files for the land surface model 
! \end{description}
!EOP

    integer              :: i 

    TRACE_ENTER("lsm_writerst")
    call lsmwrst(trim(LIS_rc%lsm)//char(0),n)

    do i=1,LIS_rc%nSubLSMs
       call sublsmwrst(trim(LIS_rc%subLSM(i))//char(0),n)
    enddo
    TRACE_EXIT("lsm_writerst")
  end subroutine LIS_lsm_writerestart

!BOP
! !ROUTINE: LIS_lsm_finalize
! \label{LIS_lsm_finalize}
!
! !INTERFACE:
  subroutine LIS_lsm_finalize()
! !USES:

!
! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures in the specific instance of a 
!  land surface model 
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmfinalize] (\ref{lsmfinalize}) \newline
!    invokes the generic method in the registry to cleanup the 
!    LSM related datastructures    
! \end{description}
!EOP

    integer              :: i 

    call lsmfinalize(trim(LIS_rc%lsm)//char(0))

    do i=1,LIS_rc%nSubLSMs
       call sublsmfinalize(trim(LIS_rc%subLSM(i))//char(0))
    enddo
  end subroutine LIS_lsm_finalize

!BOP
! !ROUTINE: lsm_setexport_noesmf
! \label{lsm_setexport_noesmf}
! 
! !INTERFACE:
  subroutine lsm_setexport_noesmf(n)
! !USES:    

! !ARGUMENTS: 
    integer, intent(in) :: n 
! !DESCRIPTION:
! This interface provides the entry point for specifying an 
! export state (a list of model specific variables) from a land surface
! model. The routine is used in a coupled simulation to provide feedback
! to a different model component such as an atmospheric model. 
!
! The arguments are: 
! \begin{description}
! \item[n]
!   index of the nest or domain
! \end{description}
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmcplsetexport] (\ref{lsmcplsetexport}) \newline
!    invokes the generic method in the registry to set the 
!    export state from the land surface model
! \end{description}
!EOP

    TRACE_ENTER("lsm_setexp")
    call lsmcplsetexport(trim(LIS_rc%lsm)//"+"&
         //trim(LIS_rc%runmode)//char(0), n)
    TRACE_EXIT("lsm_setexp")

  end subroutine lsm_setexport_noesmf




!BOP
! !ROUTINE: LIS_lsm_reset
! \label{LIS_lsm_reset}
!
! !INTERFACE:
  subroutine LIS_lsm_reset()
! !USES:

!
! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures in the specific instance of a 
!  land surface model 
!
! The calling sequence is: 
! \begin{description}
!  \item[lsmreset] (\ref{lsmreset}) \newline
!    invokes the generic method in the registry to cleanup the 
!    LSM related datastructures    
! \end{description}
!EOP
    TRACE_ENTER("lsm_reset")
    call lsmreset(trim(LIS_rc%lsm)//char(0))
    TRACE_EXIT("lsm_reset")
  end subroutine LIS_lsm_reset

!BOP
! 
! !ROUTINE: applyLSMPert
! \label{applyLSMPert}
! 
! !INTERFACE: 
  subroutine applyLSMPert(n, k, LIS_LSM_State, LIS_LSM_Pert_State)
! !USES: 

! !ARGUMENTS:     
    integer, intent(IN) :: n 
    integer, intent(IN) :: k
    type(ESMF_State)    :: LIS_LSM_State
    type(ESMF_State)    :: LIS_LSM_Pert_State
!
! !DESCRIPTION:
! 
!  This routine applies the specified perturbations to the 
!  LSM prognostic state variables. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]        index of the nest
!   \item[LSM\_State]   ESMF State with prognostic variables
!   \item[LSM\_Pert\_State]    ESMF State with prognostic variable 
!                       perturbations
!  \end{description}
!EOP

    integer                   :: i,j,t,f,m,t_unpert
    integer                   :: lsm_state_count
    integer                   :: status
    character*100,    allocatable :: lsm_state_objs(:)
    type(ESMF_Field), allocatable :: lsm_field(:)
    real, pointer             :: lsm_temp(:)
    real                      :: delta(LIS_rc%nstvars(k))
    integer, allocatable          :: typ(:)
    real, allocatable         :: stvar(:,:)
    real, allocatable         :: stvar_up(:,:)
    real, allocatable         :: pert(:,:)
    integer                   :: nxx

    allocate(stvar(LIS_rc%nstvars(k),&
         LIS_rc%npatch(n,LIS_rc%lsm_index)))
    allocate(stvar_up(LIS_rc%nstvars(k),&
         LIS_rc%npatch(n,LIS_rc%lsm_index)))
    allocate(pert(LIS_rc%nstvars(k),&
         LIS_rc%npatch(n,LIS_rc%lsm_index)))

    call ESMF_StateGet(LIS_LSM_State,itemCount=lsm_state_count,rc=status)
    call LIS_verify(status, &
         "ESMF_StateGet failed in applyLSMPert")
    
    allocate(lsm_state_objs(lsm_state_count))
    allocate(lsm_field(lsm_state_count))
    
    call ESMF_StateGet(LIS_LSM_State,itemNameList=lsm_state_objs,rc=status)
    call LIS_verify(status,&
         "ESMF_StateGet failed in applyLSMPert")        
    
    do i=1,lsm_state_count
       call ESMF_StateGet(LIS_LSM_State,lsm_state_objs(i),lsm_field(i),&
            rc=status)
       call LIS_verify(status,&
            "ESMF_StateGet failed in applyLSMPert")
       call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
       call LIS_verify(status,&
            "ESMF_FieldGet failed in applyLSMPert")
       stvar(i,:) = lsm_temp(:)
       stvar_up(i,:) = lsm_temp(:)
    enddo
    deallocate(lsm_field)

    allocate(lsm_field(lsm_state_count))
    allocate(typ(lsm_state_count))
    
    do i=1,lsm_state_count
       call ESMF_StateGet(LIS_LSM_Pert_State,&
            trim(lsm_state_objs(i)),&
            lsm_field(i),rc=status)
       call LIS_verify(status,&
            "ESMF_StateGet failed in applyLSMPert")

       call ESMF_AttributeGet(lsm_field(i),"Perturbation Type",typ(i),&
            rc=status)
       call LIS_verify(status,&
            "ESMF_AttributeGet failed in applyLSMPert")

       call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
       call LIS_verify(status,&
            "ESMF_FieldGet failed in applyLSMPert")
       pert(i,:) = lsm_temp(:)
    enddo

    if(LIS_rc%perturb_state(k).eq."uniform") then 

       nxx = LIS_rc%nensem(n)/LIS_rc%nmetforc
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
          do f=1,LIS_rc%nmetforc
             t_unpert = (i-1)*LIS_rc%nensem(n) + f*nxx
             do j=1,LIS_rc%nstvars(k)
                do m=1,nxx
                   t=(i-1)*LIS_rc%nensem(n)+(f-1)*nxx + m
!                if(typ(j).eq.0) then 
                   if(LIS_rc%pert_bias_corr.eq.1) then 
                      if(m.ne.nxx) then 
                         if(j.eq.m) then 
                            stvar(j,t) = stvar(j,t) + pert(j,t)
                         endif
                      endif
                   else
                      stvar(j,t) = stvar(j,t) + pert(j,t)
                   endif
                enddo
             enddo
          enddo
       enddo
    else

!    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!      do j=1,LIS_rc%nstvars(k)
!         if(typ(j).eq.0) then 
!            stvar(j,t) = stvar(j,t) + pert(j,t)
!         elseif(typ(j).eq.1) then 
!            stvar(j,t) = stvar(j,t) * pert(j,t)
!         endif
!      enddo
!   enddo

! If perturbation bias correction scheme is turned on, then 
! apply the perturbations to (nensem-1) ensemble members. Keep 
! the last ensemble member for each tile unperturbed. 
! This will be used later to apply the Ryu et al. JHM (2009) 
! perturbation bias correction
!
       nxx = LIS_rc%nensem(n)/LIS_rc%nmetforc
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
          do f=1,LIS_rc%nmetforc
             t_unpert = (i-1)*LIS_rc%nensem(n) + f*nxx
             do j=1,LIS_rc%nstvars(k)
                do m=1,nxx
                   t=(i-1)*LIS_rc%nensem(n)+(f-1)*nxx + m
                   if(typ(j).eq.0) then 
                      if(LIS_rc%pert_bias_corr.eq.1) then 
                         if(m.ne.nxx) then 
                            stvar(j,t) = stvar(j,t) + pert(j,t)
                         endif
                      else
                         stvar(j,t) = stvar(j,t) + pert(j,t)
                      endif
                   elseif(typ(j).eq.1) then 
                      if(LIS_rc%pert_bias_corr.eq.1) then 
                         if(m.ne.nxx) then 
                            stvar(j,t) = stvar(j,t) * pert(j,t)
                         endif
                      else
                         stvar(j,t) = stvar(j,t) * pert(j,t)
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       
       if(LIS_rc%pert_bias_corr.eq.1) then 
          do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
             do f=1,LIS_rc%nmetforc
                t_unpert = (i-1)*LIS_rc%nensem(n) + f*nxx
                do j=1,LIS_rc%nstvars(k)
                   delta(j) = 0.0
                   do m=1,nxx-1
                      t = (i-1)*LIS_rc%nensem(n)+(f-1)*nxx+m
                      delta(j) = delta(j) + (stvar(j,t)-stvar(j,t_unpert))
                   enddo
                enddo
                
                do j=1,LIS_rc%nstvars(k)
                   delta(j) = delta(j)/(nxx-1)
                   do m=1,nxx-1
                      t = (i-1)*LIS_rc%nensem(n)+(f-1)*nxx+m
                      stvar(j,t) = stvar(j,t) - delta(j)
                   enddo
                enddo
             end do
          enddo
       endif
    endif

    deallocate(lsm_state_objs)
    deallocate(lsm_field)
    deallocate(typ)
    
!set the perturbed variables 
   
   allocate(lsm_state_objs(lsm_state_count))
   allocate(lsm_field(lsm_state_count))

   call ESMF_StateGet(LIS_LSM_State,itemNameList=lsm_state_objs,rc=status)
   call LIS_verify(status,&
        "ESMF_StateGet failed in applyLSMPert")       
 
   do i=1,lsm_state_count
      call ESMF_StateGet(LIS_LSM_State,lsm_state_objs(i),lsm_field(i),&
           rc=status)
      call LIS_verify(status,&
           "ESMF_StateGet failed in applyLSMpert")
      call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
      call LIS_verify(status,&
           "ESMF_FieldGet failed in applyLSMPert")
      lsm_temp(:) = stvar(i,:) 
   enddo

   do i=1,lsm_state_count
      call ESMF_StateGet(LIS_LSM_Pert_State,&
           trim(lsm_state_objs(i)),&
           lsm_field(i),rc=status)
      call LIS_verify(status,&
           "ESMF_StateGet failed in applyLSMPert")
      
      call ESMF_FieldGet(lsm_field(i), localDE=0,farrayPtr=lsm_temp,rc=status)
      call LIS_verify(status,&
           "ESMF_FieldGet failed in applyLSMPert")
      lsm_temp(:) = stvar(i,:) - stvar_up(i,:)
    enddo

   deallocate(lsm_state_objs)
   deallocate(lsm_field)

   deallocate(stvar)
   deallocate(stvar_up)
   deallocate(pert)

 end subroutine applyLSMPert

!BOP
!
!ROUTINE: LIS_lsm_diagnoseVarsForDA
! \label{LIS_lsm_diagnoseVarsForDA}
!
! !INTERFACE:
 subroutine LIS_lsm_diagnoseVarsForDA(n)
! !ARGUMENTS:
   integer, intent(in) :: n 

! !DESCRIPTION:
! This routine diagnoses variables needed for obspred
! calculations for data assimilation from the LSM
!
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP   
   integer     :: k 

   TRACE_ENTER("lsm_diagDA")
   if(LIS_rc%ndas.gt.0) then
      do k=1, LIS_rc%nperts
         call lsmdadiagnosevars(trim(LIS_rc%lsm)//"+"//&
              trim(LIS_rc%daset(k))//char(0),n)      
      enddo
   endif
   TRACE_EXIT("lsm_diagDA")

 end subroutine LIS_lsm_diagnoseVarsForDA



!BOP
!
!ROUTINE: LIS_lsm_DAGetObsPred
! \label{LIS_lsm_DAGetObsPred}
!
! !INTERFACE:
  subroutine LIS_lsm_DAGetObsPred(n,k,Obs_Pred)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
    real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
! 
!  This interface is used to compute the "obspred" from the land surface model
!  obspred represents the model's estimate of the observation. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!   \item[Obs_Pred]   Obs pred object from the model. 
!  \end{description}
!EOP


    integer                :: m

    if(LIS_rc%LSM_DAinst_valid(k)) then
       call lsmdagetobspred(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0),n, k, Obs_pred)
    endif
  end subroutine LIS_lsm_DAGetObsPred


!BOP
!
!ROUTINE: LIS_lsm_DAGetStateVar
! \label{LIS_lsm_DAGetStateVar}
!
! !INTERFACE:
  subroutine LIS_lsm_DAGetStateVar(n,k)

! !ARGUMENTS:
    integer                :: n
    integer                :: k

!
! !DESCRIPTION:
! 
!  This interface invokes the land model to return the state vector
!  used in data assimilation. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%LSM_DAinst_valid(k)) then
       call lsmdagetstatevar(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k))
    endif
  end subroutine LIS_lsm_DAGetStateVar

!BOP
!
!ROUTINE: LIS_lsm_DASetStateVar
! \label{LIS_lsm_DASetStateVar}
!
! !INTERFACE:
  subroutine LIS_lsm_DASetStateVar(n,k)

! !ARGUMENTS:
    integer                :: n
    integer                :: k

!
! !DESCRIPTION:
! 
!  This interface invokes the land model to set the state vector
!  used in data assimilation. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call lsmdasetstatevar(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k))
    endif

  end subroutine LIS_lsm_DASetStateVar

!BOP
!
!ROUTINE: LIS_lsm_DAScaleStateVar
! \label{LIS_lsm_DAScaleStateVar}
!
! !INTERFACE:
  subroutine LIS_lsm_DAScaleStateVar(n,k)

! !ARGUMENTS:
    integer                :: n
    integer                :: k
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to allow the scaling
!  of the state vector variables
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call lsmdascalestatevar(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k))
    endif
  end subroutine LIS_lsm_DAScaleStateVar

!BOP
!
!ROUTINE: LIS_lsm_DAScaleStateVar
! \label{LIS_lsm_DAScaleStateVar}
!
! !INTERFACE:
  subroutine LIS_lsm_DADescaleStateVar(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to allow the descaling
!  of the state vector variables
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call lsmdadescalestatevar(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k), &
            LIS_LSM_Incr_State(n,k))
    endif

  end subroutine LIS_lsm_DADescaleStateVar

!BOP
!
!ROUTINE: LIS_lsm_DAUpdateState
! \label{LIS_lsm_DAUpdateState}
!
! !INTERFACE:
  subroutine LIS_lsm_DAUpdateState(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to enable the update of 
!  of the state vector variables based on the increments from the
!  DA algorithm. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call lsmdaupdatestate(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_LSM_State(n,k), &
            LIS_LSM_Incr_State(n,k))
    endif
  
  end subroutine LIS_lsm_DAUpdateState

!BOP
!
!ROUTINE: LIS_lsm_DAQCState
! \label{LIS_lsm_DAQCState}
!
! !INTERFACE:
  subroutine LIS_lsm_DAQCState(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k    
!
! !DESCRIPTION:
! 
!  This interface invokes the land model to enable the QC of
!  state vector used in DA
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!   \item[k]    index of the data assimilation instance   
!  \end{description}
!EOP

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call lsmdaqcstate(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0),n,LIS_LSM_State(n,k))
    endif
  end subroutine LIS_lsm_DAQCState

!BOP
!
!ROUTINE: LIS_lsm_DAgetStateSpaceSize
! \label{LIS_lsm_DAgetStateSpaceSize}
!
! !INTERFACE:
  subroutine LIS_lsm_DAgetStateSpaceSize(n,k,size)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
    integer                :: size
!
! !DESCRIPTION:
! 
!  This routine returns the spatial size of the state vector
!  used in DA
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]    index of the nest
!  \end{description}
!EOP
    if(LIS_rc%LSM_DAinst_valid(k)) then 
       size = LIS_rc%npatch(n,LIS_rc%lsm_index)
    endif

  end subroutine LIS_lsm_DAgetStateSpaceSize

!BOP
!
!ROUTINE: LIS_lsm_DAextractStateVector
! \label{LIS_lsm_DAextractStateVector}
!
! !INTERFACE:
  subroutine LIS_lsm_DAextractStateVector(n,k,state_size,stvar)

! !ARGUMENTS:
    integer                :: n
    integer                :: k
    integer                :: state_size
    real                   :: stvar(LIS_rc%nstvars(k),&
         LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
! 
!  This interface extracts the state vector variables from the 
!  ESMF state object. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]     index of the nest
!   \item[k]     index of the data assimilation instance   
!   \item[stvar] state vector from the LSM state
!  \end{description}
!EOP
    integer                :: status
    integer                :: v,t
    character*100,    allocatable     :: lsm_state_objs(:)
    type(ESMF_Field)                  :: lsm_field(LIS_rc%nstvars(k))
    real,         pointer             :: stdata(:)
    real,         pointer             :: stincrdata(:)

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       allocate(lsm_state_objs(LIS_rc%nstvars(k)))
       
       call ESMF_StateGet(LIS_LSM_State(n,k),itemNameList=lsm_state_objs,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_StateGet failed in enkf_increments")
    
       do v=1,LIS_rc%nstvars(k)
          call ESMF_StateGet(LIS_LSM_State(n,k),trim(lsm_state_objs(v)),&
               lsm_field(v),rc=status)
          call LIS_verify(status, &
               "ESMF_StateGet failed in enkf_increments")
          
          call ESMF_FieldGet(lsm_field(v),localDE=0, farrayPtr=stdata,rc=status)
          call LIS_verify(status,&
               "ESMF_FieldGet failed in enkf_increments")
          
          
          do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
             
             stvar(v,t) = stdata(t)     
          enddo
       enddo
       deallocate(lsm_state_objs)
    endif
  end subroutine LIS_lsm_DAextractStateVector
   
!BOP
!
!ROUTINE: LIS_lsm_DAgetFreshIncrementsStatus
! \label{LIS_lsm_DAgetFreshIncrementsStatus}
!
! !INTERFACE: 
  subroutine  LIS_lsm_DAgetFreshIncrementsStatus(n,k,setStatus)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
    logical                :: setStatus
!
! !DESCRIPTION:
! 
!  This interface gets the 'fresh increments status' attribute of
!  the LSM increments object. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]     index of the nest
!   \item[k]     index of the data assimilation instance   
!   \item[setStatus] fresh increments status attribute
!  \end{description}
!EOP
    integer                :: status

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call ESMF_AttributeGet(LIS_LSM_Incr_State(n,k),&
            "Fresh Increments Status", setstatus, rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeSet: Fresh Increments Status failed in enkf_increments')   
    endif

  end subroutine LIS_lsm_DAgetFreshIncrementsStatus

!BOP
!
!ROUTINE: LIS_lsm_DAsetFreshIncrementsStatus
! \label{LIS_lsm_DAsetFreshIncrementsStatus}
!
! !INTERFACE:  
  subroutine  LIS_lsm_DAsetFreshIncrementsStatus(n,k,setStatus)
! !ARGUMENTS:
    integer                :: n
    integer                :: k
    logical                :: setStatus
!
! !DESCRIPTION:
! 
!  This interface sets the 'fresh increments status' attribute in
!  the LSM increments object. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]     index of the nest
!   \item[k]     index of the data assimilation instance   
!   \item[setStatus] fresh increments status attribute
!  \end{description}
!EOP
    integer                :: status

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call ESMF_AttributeSet(LIS_LSM_Incr_State(n,k),&
            "Fresh Increments Status", setstatus, rc=status)
       call LIS_verify(status,&
            'ESMF_AttributeSet: Fresh Increments Status failed in enkf_increments')   
    endif
  end subroutine LIS_lsm_DAsetFreshIncrementsStatus

!BOP
!
!ROUTINE: LIS_lsm_DAsetAnlysisUpdates
! \label{LIS_lsm_DAsetAnlysisUpdates}
!
! !INTERFACE: 
  subroutine LIS_lsm_DAsetAnlysisUpdates(n,k,state_size,stvar,stincr)

! !ARGUMENTS:
    integer                :: n
    integer                :: k
    integer                :: state_size
    real                   :: stvar(LIS_rc%nstvars(k),&
         LIS_rc%npatch(n,LIS_rc%lsm_index))
    real                   :: stincr(LIS_rc%nstvars(k),&
         LIS_rc%npatch(n,LIS_rc%lsm_index))

!
! !DESCRIPTION:
! 
!  This interface sets the variables the state vector
!  and state increments vector objects after assimilation
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]      index of the nest
!   \item[k]      index of the data assimilation instance   
!   \item[stvar]  state vector variables
!   \item[stincr] state increments vector variables
!  \end{description}
!EOP

    integer                :: status
    integer                :: v,t
    character*100,    allocatable     :: lsm_state_objs(:)
    type(ESMF_Field)                  :: lsm_field(LIS_rc%nstvars(k))
    type(ESMF_Field)                  :: lsm_incr_field(LIS_rc%nstvars(k))
    real,         pointer             :: stdata(:)
    real,         pointer             :: stincrdata(:)

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       allocate(lsm_state_objs(LIS_rc%nstvars(k)))
       
       call ESMF_StateGet(LIS_LSM_State(n,k),itemNameList=lsm_state_objs,&
            rc=status)
       call LIS_verify(status, &
            "ESMF_StateGet failed in enkf_increments")
       
       do v=1,LIS_rc%nstvars(k)
          call ESMF_StateGet(LIS_LSM_State(n,k),trim(lsm_state_objs(v)),&
               lsm_field(v),rc=status)
          call LIS_verify(status, &
               "ESMF_StateGet failed in enkf_increments")
       
          call ESMF_StateGet(LIS_LSM_Incr_State(n,k),trim(lsm_state_objs(v)),&
               lsm_incr_field(v),rc=status)
          call LIS_verify(status, &
               "ESMF_StateGet failed in enkf_increments")
          
          call ESMF_FieldGet(lsm_field(v),localDE=0, farrayPtr=stdata,rc=status)
          call LIS_verify(status,&
               "ESMF_FieldGet failed in enkf_increments")
          
          call ESMF_FieldGet(lsm_incr_field(v),localDE=0,farrayPtr=stincrdata,&
               rc=status)
          call LIS_verify(status, &
               'ESMF_FieldGet failed in enkf_increments')
          
          do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
             stdata(t) =  stvar(v,t)
             stincrdata(t) = stincr(v,t)
          enddo
          
       enddo

       deallocate(lsm_state_objs)
    endif
  end subroutine LIS_lsm_DAsetAnlysisUpdates

  subroutine LIS_lsm_DAobsTransform(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call Lsmdaobstransform(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, LIS_OBS_State(n,k))       
    endif
  end subroutine LIS_lsm_DAobsTransform


  subroutine LIS_lsm_DAmapObsToLSM(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call lsmdamapobstolsm(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%daset(k))//char(0), n, k, LIS_OBS_State(n,k),&
            LIS_LSM_Incr_State(n,k))
       
    endif
  end subroutine LIS_lsm_DAmapObsToLSM


  subroutine LIS_lsm_DAqcObsState(n,k)
! !ARGUMENTS:
    integer                :: n
    integer                :: k

    if(LIS_rc%LSM_DAinst_valid(k)) then 

       call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
            //trim(LIS_rc%daset(k))//char(0),n, k, LIS_OBS_State(n,k))
       
    endif
  end subroutine LIS_lsm_DAqcObsState

!BOP
!
!ROUTINE: LIS_lsm_DAmapTileSpaceToObsSpace
! \label{LIS_lsm_DAmapTileSpaceToObsSpace}
!
! !INTERFACE:
  subroutine LIS_lsm_DAmapTileSpaceToObsSpace(n,k,tileid,st_id,en_id)

! !DESCRIPTION: 
! This routine derives the observation space location that maps to the
! input tile space of a given patch space
!
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [tileid]
!     location in the tile space 
!   \item [st\_id]
!     starting index of the observation space location
!   \item [en\_id]
!     ending index of the observation space location
!  \end{description} 
!EOP

    integer                         :: n
    integer                         :: k
    integer                         :: tileid
    integer                         :: st_id
    integer                         :: en_id

    real                            :: lat, lon, col,row
    integer                         :: gid,c,r

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
            LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%col,&
            LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
            LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%col,&
            LIS_surface(n,LIS_rc%lsm_index)%tile(tileid)%row))%lon
       
       call latlon_to_ij(LIS_obs_domain(n,k)%lisproj,lat,lon,&
            col,row)
       c = nint(col)
       r = nint(row)
       
       gid = -1
       if(c.ge.1.and.c.le.LIS_rc%obs_lnc(k).and.&
            r.ge.1.and.r.le.LIS_rc%obs_lnr(k)) then 
          gid = LIS_obs_domain(n,k)%gindex(c,r)
          
       endif
       call ij_to_latlon(LIS_obs_domain(n,k)%lisproj,real(c),real(r),&
            lat,lon)
       st_id = gid
       en_id = gid
    endif

  end subroutine LIS_lsm_DAmapTileSpaceToObsSpace

!BOP
!
!ROUTINE: LIS_lsm_DAgetStateVarNames
! \label{LIS_lsm_DAgetStateVarNames}
!
! !INTERFACE:
  subroutine LIS_lsm_DAgetStateVarNames(n,k,stateNames)
! !ARGUMENTS:
    integer            :: n
    integer            :: k
    character(len=*)   :: stateNames(LIS_rc%nstVars(k))
! 
! !DESCRIPTION:
! 
!  This interface extracts the variable names from the state vector object
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]      index of the nest
!   \item[k]      index of the data assimilation instance   
!   \item[stateNames  state vector names
!  \end{description}
!EOP
    integer            :: status

    if(LIS_rc%LSM_DAinst_valid(k)) then 
       call ESMF_StateGet(LIS_LSM_State(n,k),itemNameList=stateNames,&
            rc=status)
    endif
  end subroutine LIS_lsm_DAgetStateVarNames

!BOP
!
!ROUTINE: LIS_lsm_getlatlons
! \label{LIS_lsm_getlatlons}
!
! !INTERFACE:
  subroutine LIS_lsm_getlatlons(n,k,state_size,lats,lons)
! !ARGUMENTS:
    integer             :: n    
    integer             :: k
    integer             :: state_size
    real                :: lats(LIS_rc%npatch(n,LIS_rc%lsm_index))
    real                :: lons(LIS_rc%npatch(n,LIS_rc%lsm_index))
! 
! !DESCRIPTION:
! 
!  This routine extracts the lat/lon values corresponding to the
!  state vector space used by the LSM. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]      index of the nest
!   \item[k]      index of the data assimilation instance   
!   \item[stateNames  state vector names
!  \end{description}
!EOP

    integer             :: i,gid
    if(LIS_rc%LSM_DAinst_valid(k)) then 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          
          gid = LIS_domain(n)%gindex(&
               LIS_surface(n, LIS_rc%lsm_index)%tile(i)%col,&
               LIS_surface(n, LIS_rc%lsm_index)%tile(i)%row)
          
          lats(i) = LIS_domain(n)%grid(gid)%lat
          lons(i) = LIS_domain(n)%grid(gid)%lon
          
       enddo
    endif
    
  end subroutine LIS_lsm_getlatlons

end module LIS_lsmMod
