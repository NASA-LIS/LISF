!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_RTMMod

!BOP
!
! !MODULE: LIS_RTMMod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
!  21 Mar 2009: Sujay Kumar; Initial implementation
!
!EOP
  use ESMF 
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_RTM_init   
  public :: LIS_RTM_run
  public :: LIS_RTM_output
  public :: LIS_RTM_finalize
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LIS_sfcState
  public :: LIS_forwardState
  public :: LIS_rtm_struc

!EOP

  type(ESMF_State), allocatable :: LIS_sfcState(:)
  type(ESMF_State), allocatable :: LIS_forwardState(:)

  type, private :: rtm_type_dec
     character(len=LIS_CONST_PATH_LEN) :: outputdirname
     character(len=LIS_CONST_PATH_LEN) :: outputfilename
     character*100           :: fstatsname
     logical                 :: rtmAlarmCheck
     real                    :: rtmoutInterval
     real                    :: rtmInterval
     character*100           :: models_used
     logical                 :: stats_file_open
  end type rtm_type_dec

  type(rtm_type_dec), allocatable :: LIS_rtm_struc(:)

contains
 
!BOP
! 
! !ROUTINE: LIS_RTM_init
! \label{LIS_RTM_init}
! 
! !INTERFACE:  
  subroutine LIS_RTM_init
! 
! !DESCRIPTION: 
! 
! !USES:
    use LIS_logMod,       only : LIS_verify
    use LIS_fileIOMod,    only : LIS_create_stats_filename
    use LIS_timeMgrMod,   only : LIS_registerAlarm, LIS_parseTimeString
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
!EOP    
    implicit none
    
    integer                 :: n    
    integer                 :: yr, mo, da, hr, mn, ss
    character*1             :: nestid(2)
    character*10            :: time
    character*100           :: temp
    integer                 :: status
    character(len=LIS_CONST_PATH_LEN) :: statsfilename    

    TRACE_ENTER("rtm_init")
    allocate(LIS_rtm_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest

       LIS_rtm_struc(n)%models_used = ""
       LIS_rtm_struc(n)%stats_file_open = .true.
       call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%rtm, &
            label="Radiative transfer model:", rc=status)
       call LIS_verify(status, "Radiative transfer model: not defined")

       LIS_rtm_struc(n)%models_used = &
            trim(LIS_rtm_struc(n)%models_used)//&
            trim(LIS_rc%rtm)

       if(LIS_rc%rtm.ne."none") then 
          call ESMF_ConfigGetAttribute(LIS_config, time, &
               label="RTM invocation frequency:", rc=status)
          call LIS_verify(status, "RTM invocation frequency: not defined")
!------------------------------------------------------------------------- 
!if the invocation frequency is set to none, then the RTM is not run
! at a regular interval. It is invoked whenever the call to the forward
! model is made (e.g. from within the Data assimilation system when 
! the forward model run is required)
!------------------------------------------------------------------------- 
          if(time.ne."none") then 
             call LIS_parseTimeString(time,LIS_rtm_struc(n)%rtmInterval)

             call LIS_registerAlarm("RTM model alarm",&
                  LIS_rtm_struc(n)%rtmInterval, &
                  LIS_rtm_struc(n)%rtmInterval)
          else
             LIS_rtm_struc(n)%rtmInterval = -1
             LIS_rtm_struc(n)%rtmAlarmCheck = .true. 
          endif

          call ESMF_ConfigGetAttribute(LIS_config, time, &
               label="RTM history output frequency:", rc=status)
          call LIS_verify(status, "RTM history output frequency: not defined")
          
          call LIS_parseTimeString(time,LIS_rtm_struc(n)%rtmOutInterval)
          
          call LIS_registerAlarm("RTM output alarm",&
               LIS_rtm_struc(n)%rtmInterval, &
               LIS_rtm_struc(n)%rtmOutInterval)
! Create empty surface state object
       endif
    enddo

    allocate(LIS_sfcState(LIS_rc%nnest))
    allocate(LIS_forwardState(LIS_rc%nnest))
    
    if(LIS_rc%rtm.ne."none") then 
       do n=1,LIS_rc%nnest
          
          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid
          
          LIS_sfcState(n) = ESMF_StateCreate(name=&
               "RTM Surface State"//nestid(1)//nestid(2)&
               ,rc=status)
          call LIS_verify(status, 'Error in StateCreate in LIS_RTM_Init')
          
          LIS_forwardState(n) = ESMF_StateCreate(name=&
               "RTM Forward State"//nestid(1)//nestid(2)&
               ,rc=status)
          call LIS_verify(status, 'Error in StateCreate in LIS_RTM_Init')
          
       enddo
       
       call rtminitialize(trim(LIS_rc%rtm)//char(0))
    endif
#if 0 
!!!!!!!!!!!!!TESTING OUTPUT FILE  INITIALIZATION STEP (BEGIN)
       write(unit=outputdirname,FMT='(I2.2)') trim(LIS_rc%rtm)
       outputdirname='RTM_'//trim(outputdirname)
       write(unit=outputfilename,FMT='(I2.2)') trim(LIS_rc%rtm)
       outputfilename='RTM_'//trim(outputfilename)
       write(unit=statsfilename,FMT='(I2.2)') trim(LIS_rc%rtm)
       statsfilename='RTMStats_'//trim(statsfilename)
       call LIS_create_stats_filename(n, fstatsname,statsfilename)
#endif       
    TRACE_EXIT("rtm_init")
       
  end subroutine LIS_RTM_init

!BOP
! 
! !ROUTINE: LIS_RTM_run
! \label{LIS_RTM_run}
! 
! !INTERFACE:  
  subroutine LIS_RTM_run(n)

    implicit none
    
    integer,     intent(in) :: n 
!EOP
    
    TRACE_ENTER("rtm_run")
    call force2rtm(n)
    call surface2rtm(n)
    call geom2rtm(n)
    call RTM_run(n)
    TRACE_EXIT("rtm_run")

  end subroutine LIS_RTM_run
  
  subroutine surface2rtm(n)

    implicit none

    integer, intent(in)  :: n 

    integer             :: m
    integer             :: iret

    if(LIS_rc%rtm.ne."none") then 

       do m=1,LIS_rc%nsf_model_types
          if(LIS_rc%sf_model_type_select(m).eq.1) then 
             
             call lsm2rtm(trim(LIS_rc%rtm)//"+"//&
                  trim(LIS_rc%lsm)//char(0), n, LIS_sfcState(n))
             
          elseif(LIS_rc%sf_model_type_select(m).eq.3) then 
             
          endif
       enddo
    endif

  end subroutine surface2rtm

  subroutine force2rtm(n)

    implicit none
    
    integer,     intent(in) :: n 

    logical  :: alarmCheck
    integer  :: iret
!transfer forcing only if alarm is ringing

    if(LIS_rc%rtm.ne."none") then 
       if(LIS_rtm_struc(n)%rtmInterval.ne.-1) then 
          LIS_rtm_struc(n)%rtmAlarmCheck = .false. 
          LIS_rtm_struc(n)%rtmAlarmCheck = &
            LIS_isAlarmRinging(LIS_rc,"RTM model alarm")
       else
          LIS_rtm_struc(n)%rtmAlarmCheck = .true. 
       endif

       if(LIS_rtm_struc(n)%rtmAlarmCheck) then 
          call rtmf2t(trim(LIS_rc%rtm)//char(0), n)
       endif
    endif

  end subroutine Force2rtm

  subroutine geom2rtm(n)

    implicit none

    integer, intent(in) :: n 

    integer  :: iret

    if(LIS_rc%rtm.ne."none") then 
       if(LIS_rtm_struc(n)%rtmAlarmCheck) then 
          
          call geometry2rtm(trim(LIS_rc%rtm)//char(0), n)
          
       endif
    endif

  end subroutine Geom2rtm
  
!BOP
! 
! !ROUTINE: RTM_run
! \label{RTM_run}
! 
! !INTERFACE:  
  subroutine RTM_run(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc
    use LIS_logMod, only : LIS_logunit
! 
! !DESCRIPTION: 
! 
! !USES: 
    
    implicit none
    integer,   intent(in) :: n 
!EOP
    integer  :: iret

    if(LIS_rc%rtm.ne."none") then 

       if(LIS_rtm_struc(n)%rtmAlarmCheck) then 
          write(LIS_logunit,*) 'Running FORWARD Model'          
          call rtmrun(trim(LIS_rc%rtm)//char(0),n)
       endif
    endif

  end subroutine RTM_run

!BOP
! 
! !ROUTINE: LIS_RTM_output
! \label{LIS_RTM_output}
! 
! !INTERFACE: 
  subroutine LIS_RTM_output(n)
! !USES: 
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_histDataMod, only : LIS_histData
    use LIS_historyMod,  only : LIS_writeModelOutput
    use LIS_fileIOMod,  only : LIS_create_output_directory, &
         LIS_create_output_filename,  &
         LIS_create_stats_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    integer, intent(in)   :: n 
!EOP

    character(len=LIS_CONST_PATH_LEN) :: outfile, statsfile
    logical             :: alarmCheck
    integer             :: mo, da
    logical             :: open_stats

    TRACE_ENTER("rtm_out")
    if(LIS_rc%rtm.ne."none") then 
       
       alarmCheck = .false. 
       if(LIS_rc%time .ge. LIS_histData(n)%time) then 
          if(LIS_rc%wopt.ne."none") then 
             if(LIS_rc%output_at_specifictime.eq.1) then 
                if(LIS_histData(n)%month.eq.-1) then 
                   mo = LIS_rc%mo
                else
                   mo = LIS_histData(n)%month
                endif
                if(LIS_histData(n)%day.eq.-1) then 
                   da = LIS_rc%da
                else
                   da = LIS_histData(n)%day
                endif
                
                if(LIS_rc%mo.eq.mo.and.&
                     LIS_rc%da.eq.da.and.&
                     LIS_rc%hr.eq.LIS_histData(n)%hour.and.&
                     LIS_rc%mn.eq.LIS_histData(n)%min.and.&
                     LIS_rc%ss.eq.LIS_histData(n)%sec) then 
                   alarmCheck = .true. 
                endif
             else
                alarmCheck = LIS_isAlarmRinging(LIS_rc,&
                     "RTM output alarm")
             endif
             
             if(alarmCheck) then 
                open_stats = .false.
                if(LIS_masterproc) then 
                   call LIS_create_output_directory('RTM')
                   if (LIS_rtm_struc(n)%stats_file_open) then
                      call LIS_create_stats_filename(n,statsfile,'RTM')
                      LIS_rtm_struc(n)%stats_file_open = .false.
                      open_stats = .true.
                   endif
                endif

                call LIS_create_output_filename(n,outfile,&
                     model_name = 'RTM',&
                     writeint=LIS_rtm_struc(n)%rtmoutInterval)

                call LIS_writeModelOutput(n,outfile,statsfile,open_stats, &
                     outInterval = LIS_rtm_struc(n)%rtmoutInterval,       &
                     nsoillayers = 1,                                     &
                     lyrthk = (/1.0/),                                    &
                     nsoillayers2 = 1,                                    &
                     model_name = LIS_rtm_struc(n)%models_used,           &
                     group = 3)
             endif
          end if
       endif
    endif
    TRACE_EXIT("rtm_out")
  end subroutine LIS_RTM_output
  

!BOP
! 
! !ROUTINE: LIS_RTM_finalize
! \label{LIS_RTM_finalize}
! 
! !INTERFACE:  
  subroutine LIS_RTM_finalize
! 
! !DESCRIPTION: 
! 
! !USES: 

!EOP
  end subroutine LIS_RTM_finalize
  
end module LIS_RTMMod
