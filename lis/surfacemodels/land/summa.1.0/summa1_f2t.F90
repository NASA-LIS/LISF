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
!BOP
!
! !ROUTINE: summa1_f2t
! \label{summa1_f2t}
!
! !REVISION HISTORY: 
! 21 Jul 2004: Sujay Kumar   Initial Specification
! 23 Oct 2007: Kristi Arsenault, Implemented code for LISv5.0
! 
! !INTERFACE:
subroutine summa1_f2t(n)

! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_surface
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod 
  use LIS_logMod,        only : LIS_verify
  use LIS_timeMgrMod,    only : LIS_Calendar
  use nrtype
  use globalData
  use read_force_module
  use summa1_lsmMod
  use var_lookup,        only : iLookFORCE

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!
!  Forcing-only option (summa1) for calling the forcing transfer routines.
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

  integer            :: t,tid,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  
  integer            :: iHRU
  integer(i4b)                     :: ix_gru                     ! index of GRU that corresponds to the global HRU
  integer(i4b)                     :: ix_hru                     ! index of local HRU that corresponds to the global HRU
  character(len=1024)              :: message=''                 ! error message
  integer(i4b)                     :: err=0                      ! error code
  integer(i4b)                     :: forcNcid=-9999
  integer(i4b)                     :: forcingStep=-999           ! index of current time step in current forcing file
  integer(i4b)                     :: iFile=1                    ! index of current forcing file from forcing file list
  integer(i4b)                     :: modelTimeStep=0            ! index of model time step
  integer                          :: wgid, wtid, igru
  type(ESMF_Time)                  :: currentTime, referenceTime
  type(ESMF_TimeInterval)          :: diffTime
  integer                          :: iref_seconds
  real(dp)                         :: ref_seconds


  if(LIS_FORC_Tair%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Qair%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_SWdown%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
    call LIS_verify(status)
  endif

  if(LIS_FORC_LWdown%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
    call LIS_verify(status)
  endif

  if(LIS_FORC_Wind_E%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Wind_N%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Psurf%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Rainf%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LIS_verify(status)
  endif


  ! TODO: support accumulating forcing fields

  wgid = -9999 ! working gindex --- grid-cell being processed
  wtid = -9999 ! working tile id with respect to SUMMA's gru/hru

   summa1_struc(n)%timeStruct%var(1) = LIS_rc%yr
   summa1_struc(n)%timeStruct%var(2) = LIS_rc%mo
   summa1_struc(n)%timeStruct%var(3) = LIS_rc%da
   summa1_struc(n)%timeStruct%var(4) = LIS_rc%hr
   summa1_struc(n)%timeStruct%var(5) = LIS_rc%mn

  call ESMF_TimeSet(referenceTime, yy = 1990, &
     mm = 1, &
     dd = 1, &
     h  = 0, &
     m  = 0, &
     s  = 0, &
     calendar = LIS_calendar, &
     rc = status)
  call LIS_verify(status,'error in ESMF_TimeSet:referenceTime in summa1_f2t')

  call ESMF_TimeSet(currentTime, yy = LIS_rc%yr, &
     mm = LIS_rc%mo, &
     dd = LIS_rc%da, &
     h  = LIS_rc%hr, &
     m  = LIS_rc%mn, &
     s  = LIS_rc%ss, &
     calendar = LIS_calendar, &
     rc = status)
  call LIS_verify(status,'error in ESMF_TimeSet:currentTime in summa1_f2t')

  diffTime = currentTime - referenceTime
  call ESMF_TimeIntervalGet(diffTime, s=iref_seconds)
  !call ESMF_TimeIntervalGet(diffTime, s_i8=iref_seconds)
  ref_seconds = iref_seconds

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
     igru = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
     if ( igru /= wgid ) then
        wtid = 1
        wgid = igru
     else
        wtid = wtid + 1
     endif

     summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%time)=ref_seconds

     if(LIS_FORC_Tair%selectOpt.eq.1) then
       summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%airtemp)=tmp(tid)
     endif

     if(LIS_FORC_Qair%selectOpt.eq.1) then
       summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%spechum)=q2(tid)
     endif

     if(LIS_FORC_SWdown%selectOpt.eq.1) then
       summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%SWRadAtm)=swd(tid)
     endif

     if(LIS_FORC_LWdown%selectOpt.eq.1) then
       summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%LWRadAtm)=lwd(tid)
     endif

     if ( (LIS_FORC_Wind_E%selectOpt.eq.1) .and. &
          (LIS_FORC_Wind_N%selectOpt.eq.1) ) then
          if ( uwind(tid) == LIS_rc%udef .or. vwind(tid) == LIS_rc%udef ) then
             summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%windspd)=sqrt(uwind(tid)**2+vwind(tid)**2)
          else
             summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%windspd)=LIS_rc%udef
          endif
     endif

     if(LIS_FORC_Psurf%selectOpt.eq.1) then
       summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%airpres)=psurf(tid)
     endif

     if(pcp(tid).ne.LIS_rc%udef) then
       summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%pptrate)=pcp(tid)
     else
       summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%pptrate)=0.0
     endif

     if(LIS_FORC_CRainf%selectOpt.eq.1) then 
        if(cpcp(tid).ne.LIS_rc%udef) then 
           summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%pptrate)=summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%pptrate)+cpcp(tid)
        endif
     endif

     if(LIS_FORC_Snowf%selectOpt.eq.1) then 
        if(snowf(tid).ne.LIS_rc%udef) then
           summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%pptrate)=summa1_struc(n)%forcStruct%gru(wgid)%hru(wtid)%var(iLookFORCE%pptrate)+snowf(tid)
        endif
     endif
  enddo

!from read_force
  yearLength = 365
  if((mod(LIS_rc%yr,4) .eq. 0 .and. mod(LIS_rc%yr, 100).ne.0) &!leap year
       .or.(mod(LIS_rc%yr,400) .eq.0)) then 
     yearLength = 366
  else 
     yearLength = 365
  endif
#if 0 

  modelTimeStep=1 !for now - SVK
  ! read forcing data 
  do iHRU=1,summa1_struc(n)%nHRU  ! loop through global HRUs

  ! get mapping
     ix_gru = index_map(iHRU)%gru_ix
     ix_hru = index_map(iHRU)%ihru

  ! read forcing data
     call read_force(&
          ! input
          modelTimeStep,             & ! intent(in):    time step index
          ix_gru,                    & ! intent(in):    index of gru
          ix_hru,                    & ! intent(in):    index of LOCAL hru
          iHRU,                      & ! intent(in):    index of GLOBAL hru
          ! input-output
          iFile,                     & ! intent(inout): index of current forcing file in forcing file list
          forcingStep,               & ! intent(inout): index of read position in time dimension in current netcdf file
          forcNcid,                  & ! intent(inout): netcdf file identifier for the current forcing file
          ! output
          summa1_struc(n)%timeStruct%var,  & ! intent(out):   time data structure (integer)
          summa1_struc(n)%forcStruct%gru(ix_gru)%hru(ix_hru)%var, & ! intent(out):   forcing data structure (double precision)
          err, message)                             ! intent(out):   error control
     !  call handle_err(err,message)
  
  end do  ! (end looping through global HRUs)
#endif

end subroutine summa1_f2t
