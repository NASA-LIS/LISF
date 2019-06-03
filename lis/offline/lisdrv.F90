!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
!BOP
!
! !ROUTINE: lisdrv
!  \label{lisdrv}
!   Main program for LIS in an offline simulation
!
! !DESCRIPTION:
!  Main driver program for LIS. It creates all of the LIS components and
!  connects components if needed.
!
!  \begin{description}
!  \item[lisComp]
!      Main LIS component.
!  \item[auxComp]
!      Auxiliary LIS components.
!  \end{description}
!
! !REVISION HISTORY:
!  14Nov02    Sujay Kumar  Initial Specification
!  21Oct05    Sujay Kumar  Modified to include the runmodes. Switched
!                          to a init,run,finalize mode
!  19Jul18    Dan Rosen    Switched to ESMF components
program lisdrv
! !USES:
  use ESMF
  use lis_esmf_comp, only : lisSS => SetServices
  use LIS_cplMod, only    : cplSS => SetServices
!EOP
  implicit none

!BOC
  integer                                 :: rc, urc
  type(ESMF_GridComp)                     :: lisComp
  integer                                 :: my_gbl_id
  type(ESMF_Config)                       :: Comp_config
  logical                                 :: multiComp
  logical                                 :: writeRH
  logical                                 :: readRH
  integer                                 :: auxCount
  character(len=32), allocatable          :: auxList(:)
  type(ESMF_GridComp), allocatable        :: auxComp(:)
  type(ESMF_CplComp), allocatable         :: cplComp(:)
  integer                                 :: petListBounds(2)
  integer,allocatable                     :: petList(:)
  integer                                 :: a, p
  type(ESMF_State)                        :: lisState
  type(ESMF_State), allocatable           :: auxState(:)

  ! Initialize ESMF
!  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, &
!    logkindflag=ESMF_LOGKIND_MULTI,rc=rc)
  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, &
    logkindflag=ESMF_LOGKIND_MULTI_ON_ERROR,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_LogSet(flush=.true.,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  Comp_config = ESMF_ConfigCreate(rc=rc)
  ! Need to read command line aguments for file name
  call ESMF_ConfigLoadFile(Comp_config,"lis.config",rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_ConfigFindLabel(Comp_config,label="Auxiliary components:", &
    isPresent=multiComp, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (multiComp) then
    ! Create lisComp petList
    call ESMF_ConfigGetAttribute(Comp_config,petListBounds,&
      label="LIS petlist bounds:", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    allocate(petList(petListBounds(2)-petListBounds(1)+1),stat=rc)
    if (ESMF_LogFoundAllocError(statusToCheck=rc, &
      msg="Allocation of petList memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    do p=petListBounds(1), petListBounds(2)
      petList(p-petListBounds(1)+1) = p ! PETs are 0 based
    enddo
    ! Create lisComp
    lisComp = ESMF_GridCompCreate(name="LIS_MAIN",petList=petList,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    deallocate(petList)
    ! Set MULTICOMP Attribute
    call ESMF_AttributeSet(lisComp, name="MULTICOMP", value=multiComp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Count auxiliary components
    auxCount = ESMF_ConfigGetLen(Comp_config,label="Auxiliary components:",&
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    allocate(auxList(auxCount))
    allocate(auxComp(auxCount))
    allocate(cplComp(auxCount))
    allocate(auxState(auxCount))
    ! Get auxiliary component list
    call ESMF_ConfigGetAttribute(Comp_config,auxList,label="Auxiliary components:",&
    rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(lisComp, name="CMPLIST", valueList=auxList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! Get auxiliary component petlist and create auxComp
    call ESMF_ConfigFindLabel(Comp_config,"Aux petlist bounds:",rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    do a=1, auxCount
      call ESMF_ConfigGetAttribute(Comp_config,petListBounds(1), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_ConfigGetAttribute(Comp_config,petListBounds(2), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      allocate(petList(petListBounds(2)-petListBounds(1)+1),stat=rc)
      if (ESMF_LogFoundAllocError(statusToCheck=rc, &
        msg="Allocation of petList memory failed.", &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      do p=petListBounds(1), petListBounds(2)
        petList(p-petListBounds(1)+1) = p ! PETs are 0 based
      enddo
      ! Create auxComp
      auxComp(a) = ESMF_GridCompCreate(name="LIS_"//TRIM(auxList(a)),&
        petList=petList,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      deallocate(petList)
      ! Create cplComp
      cplComp(a) = ESMF_CplCompCreate(name="LIS_CPL_"//TRIM(auxList(a)),&
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      ! Determine if cplComp reads or writes route handle
      call ESMF_ConfigGetAttribute(Comp_config,value=writeRH,&
        label="Write RH:", default=.false., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_ConfigGetAttribute(Comp_config,readRH,&
        label="Read RH:", default=.false., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_AttributeSet(cplComp(a), name="write_rh", value=writeRH, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_AttributeSet(cplComp(a), name="read_rh", value=readRH, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    enddo
  else
    auxCount = 0
    ! Create lisComp on all pets
    lisComp = ESMF_GridCompCreate(name="LIS", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
  endif

  ! Destroy ESMF config structure after components are created
  call ESMF_ConfigDestroy(Comp_config, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! SetServices
  call ESMF_GridCompSetServices(lisComp, lisSS, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  do a=1, auxCount
    call ESMF_GridCompSetServices(auxComp(a), lisSS, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! Call setServices for the cplComp
    call ESMF_CplCompSetServices(cplComp(a), cplSS, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
  enddo

  ! Call initialize for the lisComp
  lisState = ESMF_StateCreate(rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_GridCompInitialize(lisComp, exportState=lisState, &
    userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  do a=1, auxCount
    ! Create wrtState
    auxState(a) = ESMF_StateCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_GridCompInitialize(auxComp(a), importState=auxState(a), &
      userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_CplCompInitialize(cplComp(a), importState=lisState &
      , exportState=auxState(a), userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
  enddo

  ! Call run for physics component
  if ( .NOT. writeRH ) then
    call ESMF_GridCompRun(lisComp, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    do a=1, auxCount
      call ESMF_GridCompRun(auxComp(a), userRc=urc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    enddo
  endif

  do a=1, auxCount
    ! Call finalize for coupler component
    call ESMF_CplCompFinalize(cplComp(a), userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! Call finalize for auxiliary component
    call ESMF_GridCompFinalize(auxComp(a), userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
  enddo
  ! Call finalize for physics component
  call ESMF_GridCompFinalize(lisComp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Destroy lisComp
  do a=1, auxCount
    call ESMF_CplCompDestroy(cplComp(a), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_GridCompDestroy(auxComp(a), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ! Destroy phyState
    call ESMF_StateDestroy(auxState(a), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
  enddo
  call ESMF_GridCompDestroy(lisComp, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  ! Destroy phyState
  call ESMF_StateDestroy(lisState, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (multiComp) then
    deallocate(auxState)
    deallocate(cplComp)
    deallocate(auxComp)
    deallocate(auxList)
  endif

  ! Finalize ESMF
  call ESMF_Finalize()
!EOC
end program lisdrv


