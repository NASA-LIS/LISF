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
!  Main driver program for LIS. It performs four main functions
!
!  \begin{description}
!  \item[LIS\_config\_init]  
!      calls the routines to read the runtime configurations
!  \item[LIS\_Init]  
!      calls the initializations based on the runmode
!  \item[LIS\_run]
!      calls the run routines based on the runmode
!  \item[LIS\_finalize] 
!     calls the cleanup routines 
!  \end{description}
!
!  The routines invoked are :
!  \begin{description}
!   \item[LIS\_config\_init](\ref{LIS_config_init}) \newline
!    call to initialize configuration tool and read model independent
!    options
!   \item[lisinit](\ref{lisinit}) \newline
!    call to initialize lis based on the runmode
!   \item[lisrun](\ref{lisrun}) \newline
!    call to run lis based on the runmode
!   \item[lisfinalize](\ref{lisfinalize}) \newline
!    call to cleanup lis structures based on the runmode
!  \end{description}
! !REVISION HISTORY: 
!  14Nov02    Sujay Kumar  Initial Specification
!  21Oct05    Sujay Kumar  Modified to include the runmodes. Switched
!                          to a init,run,finalize mode
module lis_esmf_comp
! !USES:       
  use LIS_coreMod, only : LIS_config_init, LIS_rc
  use LIS_histDataMod, only : LIS_histFB, LIS_histSend, LIS_histRecv
  use sfoutput_runMod
  use ESMF
!EOP
  implicit none

  private

  public SetServices

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

!BOC
  subroutine SetServices(lis_cmp, rc)
    type(ESMF_GridComp)  :: lis_cmp
    integer, intent(out) :: rc
    ! LOCAL VARIABLES
    character(32) :: cname

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(lis_cmp, name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (cname .eq. 'LIS_SFOUTPUT') then
      call ESMF_GridCompSetEntryPoint(lis_cmp, ESMF_METHOD_INITIALIZE, &
        userRoutine=esmf_lis_sfoutput_ini, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_GridCompSetEntryPoint(lis_cmp, ESMF_METHOD_RUN, &
        userRoutine=esmf_lis_sfoutput_run, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_GridCompSetEntryPoint(lis_cmp, ESMF_METHOD_FINALIZE, &
        userRoutine=esmf_lis_sfoutput_fin, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    else
      call ESMF_GridCompSetEntryPoint(lis_cmp, ESMF_METHOD_INITIALIZE, &
        userRoutine=esmf_lis_main_ini, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_GridCompSetEntryPoint(lis_cmp, ESMF_METHOD_RUN, &
        userRoutine=esmf_lis_main_run, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_GridCompSetEntryPoint(lis_cmp, ESMF_METHOD_FINALIZE, &
        userRoutine=esmf_lis_main_fin, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif

  end subroutine SetServices

  !-----------------------------------------------------------------------------
  ! LIS Main Component - Ini, Run, Fin
  !-----------------------------------------------------------------------------

  subroutine esmf_lis_main_ini(lis_cmp,imp_state,exp_state,clock,rc)
    ! ARGUMENTS
    type(ESMF_GridComp) :: lis_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc
    ! LOCAL VARIABLES
    type(ESMF_VM)                  :: vm
    character(len=32)              :: cname
    logical                        :: multiComp
    integer                        :: cmpCount
    character(len=32), allocatable :: cmpList(:)
    integer                        :: i

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(lis_cmp, vm=vm, name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Initialize Applicaiton
    call LIS_config_init(cmd_args=.true., vm=vm)
    TRACE_ENTER("LIS_init")
    call lisinit(trim(LIS_rc%runmode)//char(0))
    TRACE_EXIT("LIS_init")

    call ESMF_AttributeGet(lis_cmp, name="MULTICOMP", &
      value=multiComp, defaultvalue=.FALSE., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (multiComp) then
      call ESMF_AttributeGet(lis_cmp, name="CMPLIST", &
        itemCount=cmpCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      allocate(cmpList(cmpCount),stat=rc)
      if (ESMF_LogFoundAllocError(statusToCheck=rc, &
        msg="Allocation of cmpList memory failed.", &
        line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeGet(lis_cmp, name="CMPLIST", &
        valueList=cmpList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      do i=1, cmpCount
        if (cmpList(i) .eq. "SFOUTPUT") then
          call ESMF_LogWrite(trim(cname)//": Connecting to SFOUTPUT",ESMF_LOGMSG_INFO)
          LIS_histSend = .TRUE.
          ! Store field bundle list in exp_state
          call ESMF_StateAdd(exp_state, fieldbundleList=LIS_histFB, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
        endif
      enddo 

      deallocate(cmpList,stat=rc)
      if (ESMF_LogFoundDeallocError(statusToCheck=rc, &
        msg="Dellocation of cmpList memory failed.", &
        line=__LINE__, file=__FILE__)) return

    endif

  end subroutine esmf_lis_main_ini

  !-----------------------------------------------------------------------------

  subroutine esmf_lis_main_run(lis_cmp,imp_state,exp_state,clock,rc)
    ! ARGUMENTS
    type(ESMF_GridComp) :: lis_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Run Application
    TRACE_ENTER("LIS_run")
    call lisrun(trim(LIS_rc%runmode)//char(0))
    TRACE_EXIT("LIS_run")
  end subroutine esmf_lis_main_run

  !-----------------------------------------------------------------------------

  subroutine esmf_lis_main_fin(lis_cmp,imp_state,exp_state,clock, rc)
    ! ARGUMENTS
    type(ESMF_GridComp) :: lis_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Finalize Application
    call lisfinalize(trim(LIS_rc%runmode)//char(0))
  end subroutine esmf_lis_main_fin

  !-----------------------------------------------------------------------------
  ! LIS Surface Model Output Component - Ini, Run, Fin
  !-----------------------------------------------------------------------------

  subroutine esmf_lis_sfoutput_ini(lis_cmp,imp_state,exp_state,clock,rc)
    ! USES
    use sfoutput_runMod
    ! ARGUMENTS
    type(ESMF_GridComp) :: lis_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc
    ! LOCAL VARIABLES
    type(ESMF_VM)                  :: vm
    character(ESMF_MAXSTR)         :: cname
    logical                        :: sendComp
    logical                        :: recvComp

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(lis_cmp, vm=vm, name=cname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Initialize Applicaiton
    call LIS_config_init(cmd_args=.true., vm=vm)
    LIS_rc%decompose_by_processes = .true.
    LIS_rc%npesx = 0
    LIS_rc%npesy = 0
    TRACE_ENTER("LIS_init")
    call lis_init_sfoutput()
    TRACE_EXIT("LIS_init")

    ! Store field bundle list in exp_state
    call ESMF_StateAdd(imp_state, fieldbundleList=LIS_histFB, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    LIS_histRecv = .TRUE.

  end subroutine esmf_lis_sfoutput_ini

  !-----------------------------------------------------------------------------

  subroutine esmf_lis_sfoutput_run(lis_cmp,imp_state,exp_state,clock,rc)
    ! USES
    use sfoutput_runMod 
    ! ARGUMENTS
    type(ESMF_GridComp) :: lis_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Run Application
    TRACE_ENTER("LIS_run")
    call lis_run_sfoutput() 
    TRACE_EXIT("LIS_run")
  end subroutine esmf_lis_sfoutput_run

  !-----------------------------------------------------------------------------

  subroutine esmf_lis_sfoutput_fin(lis_cmp,imp_state,exp_state,clock, rc)
    ! USES
    use sfoutput_runMod
    ! ARGUMENTS
    type(ESMF_GridComp) :: lis_cmp
    type(ESMF_State)    :: imp_state
    type(ESMF_State)    :: exp_state
    type(ESMF_Clock)    :: clock
    integer,intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Finalize Application
    call lis_final_sfoutput()
  end subroutine esmf_lis_sfoutput_fin

  !-----------------------------------------------------------------------------

!EOC
end module lis_esmf_comp

