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
module LIS_cplMod
!BOP
!
! !MODULE: LIS_cplMod
!
! !REVISION HISTORY: 
! 06 Aug 2018 : Daniel Rosen: Initial Version
! 
! !USES:
  use ESMF
! !DESCRIPTION: 
!
! 
!
!EOP
  implicit none
  
  private

  public SetServices
  public LIS_cplRun
  public LIS_cplBcast

  type type_coupledFB
    type(ESMF_FieldBundle) :: srcFields
    type(ESMF_FieldBundle) :: dstFields
    type(ESMF_RouteHandle) :: redistRH
  end type

  type(ESMF_VM) :: LIS_cpl_vm
  type(type_coupledFB),allocatable :: coupled_histFB(:)
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(cpl_cmp, rc)
    type(ESMF_CplComp)   :: cpl_cmp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call ESMF_CplCompSetEntryPoint(cpl_cmp, ESMF_METHOD_INITIALIZE, &
      userRoutine=esmf_cpl_ini, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_CplCompSetEntryPoint(cpl_cmp, ESMF_METHOD_FINALIZE, &
      userRoutine=esmf_cpl_fin, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine SetServices

  subroutine esmf_cpl_ini(cpl_cmp,imp_state,exp_state,clock,rc)
    ! ARGUMENTS
    type(ESMF_CplComp)  :: cpl_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc
    ! LOCAL VARIABLES
    integer                                 :: itemCount
    character (len=ESMF_MAXSTR),allocatable :: itemNameList(:)
    type(ESMF_StateItem_Flag),allocatable   :: itemTypeList(:)
    integer                                 :: srcFBcount
    integer                                 :: dstFBcount
    integer                                 :: iIndex
    integer                                 :: bIndex
    logical                                 :: writeRH
    logical                                 :: readRH
    character (len=3)                       :: fbStr

    rc = ESMF_SUCCESS

    TRACE_ENTER("LIS_cpl_init")
    ! Need to reconcile import and export states
    call ESMF_CplCompGet(cpl_cmp, vm=LIS_cpl_vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(cpl_cmp, name="write_rh", value=writeRH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeGet(cpl_cmp, name="read_rh", value=readRH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    TRACE_ENTER("LIS_cpl_rcle")
    call ESMF_StateReconcile(imp_state, vm=LIS_cpl_vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_StateReconcile(exp_state, vm=LIS_cpl_vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    TRACE_EXIT("LIS_cpl_rcle")

    TRACE_ENTER("LIS_cpl_getFB")
    ! Get fieldbundle from imp_state
    call ESMF_StateGet(imp_state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    allocate(itemNameList(itemCount),stat=rc)
    if (ESMF_LogFoundAllocError(statusToCheck=rc, &
      msg="Allocation of itemNameList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
    allocate(itemTypeList(itemCount),stat=rc)
    if (ESMF_LogFoundAllocError(statusToCheck=rc, &
      msg="Allocation of itemTypeList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(imp_state, itemorderflag=ESMF_ITEMORDER_ADDORDER &
      , itemNameList=itemNameList, itemTypeList=itemTypeList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    srcFBcount=0
    do iIndex=1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELDBUNDLE ) then
        srcFBcount=srcFBcount + 1
      endif
    enddo
    allocate(coupled_histFB(srcFBcount),stat=rc)
    if (ESMF_LogFoundAllocError(statusToCheck=rc, &
      msg="Allocation of coupled_histFB memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out

    bIndex=1
    do iIndex=1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELDBUNDLE ) then
        call ESMF_StateGet(imp_state, itemName=itemNameList(iIndex) &
          , fieldbundle=coupled_histFB(bIndex)%srcFields, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        bIndex=bIndex + 1
      endif
    enddo


    deallocate(itemNameList,stat=rc)
    if (ESMF_LogFoundDeallocError(statusToCheck=rc, &
      msg="Deallocation of itemNameList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
    deallocate(itemTypeList,stat=rc)
    if (ESMF_LogFoundDeallocError(statusToCheck=rc, &
      msg="Deallocation of itemTypeList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out

    ! Get fieldbundle from exp_state
    call ESMF_StateGet(exp_state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(itemNameList(itemCount),stat=rc)
    if (ESMF_LogFoundAllocError(statusToCheck=rc, &
      msg="Allocation of itemNameList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
    allocate(itemTypeList(itemCount),stat=rc)
    if (ESMF_LogFoundAllocError(statusToCheck=rc, &
      msg="Allocation of itemTypeList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(exp_state, itemorderflag=ESMF_ITEMORDER_ADDORDER &
      , itemNameList=itemNameList, itemTypeList=itemTypeList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    dstFBcount=0
    do iIndex=1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELDBUNDLE ) then
        dstFBcount=dstFBcount + 1
      endif
    enddo
    if ( srcFBcount.ne.dstFBcount) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_INCONS, &
        msg="Source FB count does not match Destination FB count.", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return
    endif
    bIndex=1
    do iIndex=1, itemCount
      if ( itemTypeList(iIndex) == ESMF_STATEITEM_FIELDBUNDLE ) then
        call ESMF_StateGet(exp_state, itemName=itemNameList(iIndex), &
          fieldbundle=coupled_histFB(bIndex)%dstFields, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        bIndex=bIndex + 1
      endif
    enddo
    deallocate(itemNameList,stat=rc)
    if (ESMF_LogFoundDeallocError(statusToCheck=rc, &
      msg="Deallocation of itemNameList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
    deallocate(itemTypeList,stat=rc)
    if (ESMF_LogFoundDeallocError(statusToCheck=rc, &
      msg="Deallocation of itemTypeList memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
    TRACE_EXIT("LIS_cpl_getFB")

    if ( writeRH ) then
      TRACE_ENTER("LIS_cpl_rdwt")
      call ESMF_LogWrite("LIS_cplMod: Redist Store and RH Write")
      do bIndex=lbound(coupled_histFB,1), ubound(coupled_histFB,1)
        write (fbStr, "(I0)") bIndex
        TRACE_ENTER("LIS_cpl_genrh")
        call ESMF_FieldBundleRedistStore(coupled_histFB(bIndex)%srcFields, &
          coupled_histFB(bIndex)%dstFields, &
          routehandle=coupled_histFB(bIndex)%redistRH, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        TRACE_EXIT("LIS_cpl_genrh")
        TRACE_ENTER("LIS_cpl_wrtrh")
        call ESMF_RouteHandleWrite(coupled_histFB(bIndex)%redistRH, &
          fileName="LIS_CPL.RH", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        TRACE_EXIT("LIS_cpl_wrtrh")
      enddo
      TRACE_EXIT("LIS_cpl_rdwt")
    elseif ( readRH ) then
      call ESMF_LogWrite("LIS_cplMod: RH Read")
      TRACE_ENTER("LIS_cpl_rdrd")
      do bIndex=lbound(coupled_histFB,1), ubound(coupled_histFB,1)
        write (fbStr, "(I0)") bIndex
        coupled_histFB(bIndex)%redistRH = ESMF_RouteHandleCreate( &
          fileName="LIS_CPL.RH", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite("LIS_cplMod: Finished RH Read")
      enddo
      TRACE_EXIT("LIS_cpl_rdrd")
    else
      TRACE_ENTER("LIS_cpl_rdst")
      call ESMF_LogWrite("LIS_cplMod: Redist Store")
      do bIndex=lbound(coupled_histFB,1), ubound(coupled_histFB,1)
        call ESMF_FieldBundleRedistStore(coupled_histFB(bIndex)%srcFields, &
          coupled_histFB(bIndex)%dstFields, &
          routehandle=coupled_histFB(bIndex)%redistRH, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite("LIS_cplMod: Finished Redist Store")
      enddo
      TRACE_EXIT("LIS_cpl_rdst")
    endif

    TRACE_EXIT("LIS_cpl_init")
  end subroutine esmf_cpl_ini

  !-----------------------------------------------------------------------------

  subroutine LIS_cplRun(n,rc)
    ! ARGUMENTS
    integer,intent(in),optional :: n
    integer,intent(out)         :: rc
    ! LOCAL VARIABLES
    integer :: bIndex

    rc = ESMF_SUCCESS

    TRACE_ENTER("LIS_cpl_run")
    if ( present(n) ) then
      if ( n.lt.lbound(coupled_histFB,1) .OR. n.gt.ubound(coupled_histFB,1) ) then
        call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
          msg="Coupling nest number is out of range.", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return
      endif
      call ESMF_FieldBundleRedist(coupled_histFB(n)%srcFields, &
        coupled_histFB(n)%dstFields, &
        routehandle=coupled_histFB(n)%redistRH, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      do bIndex=lbound(coupled_histFB,1), ubound(coupled_histFB,1)
        call ESMF_FieldBundleRedist(coupled_histFB(bIndex)%srcFields, &
          coupled_histFB(bIndex)%dstFields, &
          routehandle=coupled_histFB(bIndex)%redistRH, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      enddo
    endif
    TRACE_EXIT("LIS_cpl_run")
  end subroutine LIS_cplRun

  !-----------------------------------------------------------------------------

  subroutine LIS_cplBcast(value,rc)
    ! ARGUMENTS
    integer,intent(inout)  :: value
    integer,intent(out)    :: rc
    ! LOCAL VARIABLES
    integer,dimension(1)   :: value1D(1)

    rc = ESMF_SUCCESS

    TRACE_ENTER("LIS_cpl_bcst")
    value1D(1) = value
    call ESMF_VMBroadcast(LIS_cpl_vm, bcstData=value1D, count=1, &
      rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    value = value1D(1)
    TRACE_EXIT("LIS_cpl_bcst")
  end subroutine LIS_cplBcast

  !-----------------------------------------------------------------------------


  subroutine esmf_cpl_fin(cpl_cmp,imp_state,exp_state,clock, rc)
    ! ARGUMENTS
    type(ESMF_CplComp)  :: cpl_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc
    ! LOCAL VARIABLES
    integer :: bIndex

    rc = ESMF_SUCCESS

    do bIndex=lbound(coupled_histFB,1), ubound(coupled_histFB,1)
      call ESMF_FieldBundleRedistRelease(coupled_histFB(bIndex)%redistRH, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo
    deallocate(coupled_histFB,stat=rc)
    if (ESMF_LogFoundDeallocError(statusToCheck=rc, &
      msg="Deallocation of coupled_histFB memory failed.", &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return  ! bail out
  end subroutine esmf_cpl_fin

  !-----------------------------------------------------------------------------

end module LIS_cplMod
