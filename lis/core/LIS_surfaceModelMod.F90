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

module LIS_surfaceModelMod
  use ESMF
  use LIS_coreMod
  use LIS_lakemodelMod
  use LIS_glaciermodelMod
  use LIS_openwatermodelMod
  use LIS_surfaceModelDataMod
  use LIS_lsmMod
  use LIS_routingMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_RTMMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LIS_surfaceModel_init 
  public :: LIS_surfaceModel_setup   
  public :: LIS_surfaceModel_readrestart
  public :: LIS_surfaceModel_run
  public :: LIS_surfaceModel_f2t
  public :: LIS_surfaceModel_output
  public :: LIS_surfaceModel_writerestart
  public :: LIS_surfaceModel_perturb_states
  public :: LIS_surfaceModel_finalize
  public :: LIS_surfaceModel_reset
  public :: LIS_surfaceModel_diagnoseVarsforDA
  public :: LIS_surfaceModel_setexport
  public :: LIS_surfaceModel_DAGetObsPred
  public :: LIS_surfaceModel_DAGetStateVar
  public :: LIS_surfaceModel_DASetStateVar
  public :: LIS_surfaceModel_DAScaleStateVar
  public :: LIS_surfaceModel_DADescaleStateVar
  public :: LIS_surfaceModel_DAUpdateState
  public :: LIS_surfaceModel_DAQCState
  public :: LIS_surfaceModel_DAgetStateSpaceSize
  public :: LIS_surfaceModel_DAextractStateVector
  public :: LIS_surfaceModel_DASetFreshIncrementsStatus
  public :: LIS_surfaceModel_DAGetFreshIncrementsStatus
  public :: LIS_surfaceModel_DAsetAnlysisUpdates
  public :: LIS_surfaceModel_DAmapTileSpaceToObsSpace
  public :: LIS_surfaceModel_DAgetStateVarNames
  public :: LIS_surfaceModel_DAobsTransform
  public :: LIS_surfaceModel_DAmapObsToModel
  public :: LIS_surfaceModel_DAqcObsState
  public :: LIS_surfaceModel_getlatlons

!BOP
! !ROUTINE: LIS_surfaceModel_setexport
! \label{LIS_surfaceModel_setexport}
!
! !INTERFACE: 
  interface LIS_surfaceModel_setexport
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure surfaceModel_setexport_noesmf
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
! !ROUTINE: LIS_surfaceModel_init
! \label{LIS_surfaceModel_init}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_init

! !DESCRIPTION:
! This routine initialized the surface models.
!EOP
        
    integer             :: m
    
    character*10 :: time
    integer      :: n
    character*3  :: fnest
    integer      :: rc

    TRACE_ENTER("sf_init")
    allocate(LIS_sfmodel_struc(LIS_rc%nnest))

    call ESMF_ConfigFindLabel(LIS_config,"Surface model output interval:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
       call LIS_verify(rc,'Surface model output interval: not defined')
       
       write(fnest,'(i3.3)') n
       if ( time == "dekad" ) then
          LIS_sfmodel_struc(n)%outInterval = 864000.0
          LIS_sfmodel_struc(n)%outIntervalType = "dekad"
          call LIS_registerAlarm("LIS surface model output alarm "//trim(fnest),&
                                 LIS_sfmodel_struc(n)%ts,&
                                 LIS_sfmodel_struc(n)%outInterval,&
                                 intervalType="dekad", dek_offset=0, when="end")
       else
          call LIS_parseTimeString(time,LIS_sfmodel_struc(n)%outInterval)
          LIS_sfmodel_struc(n)%outIntervalType = ""
    
          call LIS_registerAlarm("LIS surface model output alarm "//trim(fnest),&
                                 LIS_sfmodel_struc(n)%ts,&
                                 LIS_sfmodel_struc(n)%outInterval)
       endif

       LIS_sfmodel_struc(n)%models_used = ""
       LIS_sfmodel_struc(n)%stats_file_open = .true.
    enddo

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_init()
!          LIS_sfmodel_struc(n)%models_used = &
!               trim(LIS_sfmodel_struc(n)%models_used)//&
!               trim(LIS_rc%lsm)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          call LIS_lakemodel_init()
!          LIS_sfmodel_struc(n)%models_used = &
!               trim(LIS_sfmodel_struc(n)%models_used)//&
!               "+"//trim(LIS_rc%lakemodel)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 
          call LIS_glaciermodel_init()
!          LIS_sfmodel_struc(n)%models_used = &
!               trim(LIS_sfmodel_struc(n)%models_used)//&
!               "+"//trim(LIS_rc%lakemodel)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%openwater_index) then 
          call LIS_openwatermodel_init()
!          LIS_sfmodel_struc(n)%models_used = &
!               trim(LIS_sfmodel_struc(n)%models_used)//&
!               "+"//trim(LIS_rc%openwatermodel)
       endif
    enddo
    TRACE_EXIT("sf_init")
  end subroutine LIS_surfaceModel_init


!BOP
! 
! !ROUTINE: LIS_surfaceModel_setup
! \label{LIS_surfaceModel_setup}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_setup

! !DESCRIPTION:
! This routine sets up the surface models.
!EOP
    
    integer             :: m
 
    TRACE_ENTER("sf_setup")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_setuplsm()
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          call LIS_setuplakemodel()
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 
          call LIS_setupglaciermodel()
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%openwater_index) then 
          call LIS_setupopenwatermodel()
       endif
    enddo
    TRACE_EXIT("sf_setup")
    
  end subroutine LIS_surfaceModel_setup

!BOP
! 
! !ROUTINE: LIS_surfaceModel_readrestart
! \label{LIS_surfaceModel_readrestart}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_readrestart
! !DESCRIPTION:
! This routine reads the restart data for each
! surface model.
!EOP

    integer             :: m

    TRACE_ENTER("sf_readrst")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_readrestart
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          call LIS_lakemodel_readrestart()
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 
          call LIS_glaciermodel_readrestart()
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%openwater_index) then 
          call LIS_openwatermodel_readrestart()
       endif
    enddo
    TRACE_EXIT("sf_readrst")

  end subroutine LIS_surfaceModel_readrestart

!BOP
! 
! !ROUTINE: LIS_surfaceModel_run
! \label{LIS_surfaceModel_run}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_run(n)
! !ARGUMENTS: 
    integer, intent(in)   :: n 

! !DESCRIPTION:
! This routine runs the surface models.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    integer             :: m

    TRACE_ENTER("sf_run")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_run(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          call LIS_lakemodel_run(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 
          call LIS_glaciermodel_run(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%openwater_index) then 
          call LIS_openwatermodel_run(n)
       endif
    enddo
    TRACE_EXIT("sf_run")

  end subroutine LIS_surfaceModel_run

!BOP
! 
! !ROUTINE: LIS_surfaceModel_f2t
! \label{LIS_surfaceModel_f2t}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_f2t(n)
! !ARGUMENTS: 
    integer, intent(in)   :: n 

! !DESCRIPTION:
! This routine transfers met forcing to the surface model tiles.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    integer             :: m

    TRACE_ENTER("sf_f2t")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_f2t(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          call LIS_lakemodel_f2t(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 
          call LIS_glaciermodel_f2t(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%openwater_index) then 
          call LIS_openwatermodel_f2t(n)
       endif
    enddo
    TRACE_EXIT("sf_f2t")

  end subroutine LIS_surfaceModel_f2t

!BOP
! 
! !ROUTINE: LIS_surfaceModel_output
! \label{LIS_surfaceModel_output}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_output(n)
! !USES: 

    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_histDataMod, only : LIS_histData
    use LIS_historyMod,  only : LIS_writeModelOutput
    use LIS_fileIOMod,  only : LIS_create_output_directory, &
         LIS_create_output_filename,  &
         LIS_create_stats_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

! !ARGUMENTS: 
    integer, intent(in)   :: n 

! !DESCRIPTION:
! This routine writes the surface model output.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: outfile, statsfile
    logical             :: alarmCheck
    integer             :: mo, da
    logical             :: open_stats
    character*3         :: fnest

    TRACE_ENTER("sf_output")
    alarmCheck = .false. 
    if(LIS_rc%time .ge. LIS_histData(n)%time) then 

       if(LIS_rc%wopt.ne."none") then 
          if(LIS_rc%output_at_specifictime.eq.LIS_rc%lsm_index) then 
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
             write(fnest,'(i3.3)') n
             if ( LIS_sfmodel_struc(n)%outIntervalType == "dekad" ) then
                alarmCheck = LIS_isAlarmRinging(LIS_rc,&
                             "LIS surface model output alarm "//trim(fnest), &
                             LIS_sfmodel_struc(n)%outIntervalType)
             else
                alarmCheck = LIS_isAlarmRinging(LIS_rc,&
                             "LIS surface model output alarm "//trim(fnest))
             endif
          endif
       
          if(alarmCheck) then 
             open_stats = .false.
             call LIS_create_output_directory('SURFACEMODEL')             
             if(LIS_masterproc) then 
                if (LIS_sfmodel_struc(n)%stats_file_open) then
                   call LIS_create_stats_filename(n,statsfile,'SURFACEMODEL')
                   LIS_sfmodel_struc(n)%stats_file_open = .false.
                   open_stats = .true.
                endif
             endif

             call LIS_create_output_filename(n,outfile,&
                  model_name = 'SURFACEMODEL',&
                  writeint=LIS_sfmodel_struc(n)%outInterval)

             ! hkb-- added second set of soil layer thickness for CLSM
             if ( LIS_sfmodel_struc(n)%nsm_layers .eq.  &
                  LIS_sfmodel_struc(n)%nst_layers ) then   
                call LIS_writeModelOutput(n,outfile,statsfile, open_stats, &
                     outInterval=LIS_sfmodel_struc(n)%outInterval,         &
                     nsoillayers = LIS_sfmodel_struc(n)%nsm_layers,        &
                     lyrthk = LIS_sfmodel_struc(n)%lyrthk,                 &
                     nsoillayers2 = LIS_sfmodel_struc(n)%nst_layers,       &
                     model_name=LIS_sfmodel_struc(n)%models_used)
             else
                call LIS_writeModelOutput(n,outfile,statsfile, open_stats, &
                     outInterval=LIS_sfmodel_struc(n)%outInterval,         &
                     nsoillayers = LIS_sfmodel_struc(n)%nsm_layers,        &
                     lyrthk = LIS_sfmodel_struc(n)%lyrthk,                 &
                     nsoillayers2 = LIS_sfmodel_struc(n)%nst_layers,       &
                     model_name=LIS_sfmodel_struc(n)%models_used,          &
                     lyrthk2 = LIS_sfmodel_struc(n)%lyrthk2)
             endif
          endif
       end if
    endif
    TRACE_EXIT("sf_output")
  end subroutine LIS_surfaceModel_output

!BOP
! 
! !ROUTINE: LIS_surfaceModel_writerestart
! \label{LIS_surfaceModel_writerestart}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_writerestart(n)
! !ARGUMENTS: 
    integer, intent(in)   :: n 

! !DESCRIPTION:
! This routine writes the surface model restart data.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    integer             :: m

    TRACE_ENTER("sf_writerst")
    if(LIS_rc%wopt_rst.ne.0) then 
       do m=1,LIS_rc%nsf_model_types
          if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
             call LIS_lsm_writerestart(n)
          elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
             call LIS_lakemodel_writerestart(n)
          elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 
             call LIS_glaciermodel_writerestart(n)
          elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%openwater_index) then 
             call LIS_openwatermodel_writerestart(n)
          endif
       enddo
    elseif( LIS_rc%wopt_rst == 0 .and. LIS_rc%tscount(n) == 1 ) then
       write(LIS_logunit,*) "MSG:  No restart files being written ... "
       write(LIS_logunit,*) "MSG: ... if you want to write restart files, set:"
       write(LIS_logunit,*) "MSG: 'Output model restart files' to '1'."
       write(LIS_logunit,*) " "
    endif
    TRACE_EXIT("sf_writerst")

  end subroutine LIS_surfaceModel_writerestart

!BOP
! 
! !ROUTINE: LIS_surfaceModel_perturb_states
! \label{LIS_surfaceModel_perturb_states}
! 
! !INTERFACE: 
  subroutine LIS_surfaceModel_perturb_states(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc
! !ARGUMENTS: 
    integer,    intent(IN)  :: n 
!
! !DESCRIPTION:
! This routine perturbs the surface model states.
!
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP
    
    integer             :: m

    TRACE_ENTER("sf_perturb")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_perturb_states(n)
          call LIS_routing_perturb_states(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%openwater_index) then 

       endif
    enddo
    TRACE_EXIT("sf_perturb")
  end subroutine LIS_surfaceModel_perturb_states

!BOP
! 
! !ROUTINE: LIS_surfaceModel_diagnoseVarsforDA
! \label{LIS_surfaceModel_diagnoseVarsforDA}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_diagnoseVarsforDA(n)
    
! !ARGUMENTS: 
    integer,    intent(IN)  :: n 

! !DESCRIPTION:
! This routine diagnoses variables needed for obspred
! calculations for data assimilation.
!
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    integer             :: m

    TRACE_ENTER("sf_diagDAout")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_diagnoseVarsforDA(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
    TRACE_EXIT("sf_diagDAout")

  end subroutine LIS_surfaceModel_diagnoseVarsforDA

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAGetObsPred
! \label{LIS_surfaceModel_DAGetObsPred}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAGetObsPred(n,k,Obs_Pred)
! 
! !DESCRIPTION:
!   This routine computes the obspred for a given data assimilation 
!   scheme. Obspred is the model's estimate of the observation 
!   used in data assimilation. 
!EOP
    integer                :: n
    integer                :: k
    real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 

          call LIS_lsm_DAGetObsPred(n,k,Obs_Pred)
          call LIS_routing_DAGetObsPred(n,k,Obs_pred)

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAGetObsPred

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAGetStateVar
! \label{LIS_surfaceModel_DAGetStateVar}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAGetStateVar(n,k)
! 
! !DESCRIPTION:
! 
!  This routine retrieves the state variables used in data assimilation
!  from the model
!
!EOP
    integer                :: n
    integer                :: k
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAGetStateVar(n,k)
          call LIS_routing_DAGetStateVar(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAGetStateVar

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DASetStateVar
! \label{LIS_surfaceModel_DASetStateVar}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DASetStateVar(n,k)
! 
! !DESCRIPTION:
! 
!  This routine sets the state variables used in data assimilation
!  from the model
!
!EOP
    integer                :: n
    integer                :: k
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DASetStateVar(n,k)
          call LIS_routing_DASetStateVar(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DASetStateVar

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAScaleStateVar
! \label{LIS_surfaceModel_DAScaleStateVar}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAScaleStateVar(n,k)
! 
! !DESCRIPTION:
! 
!  This routine scales the state vector used in data assimilation. 
!  Scaling and descaling is done to enhance the numerical stability of 
!  matrix calculations. 
!
!EOP
    integer                :: n
    integer                :: k
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAScaleStateVar(n,k)
          call LIS_routing_DAScaleStateVar(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAScaleStateVar

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DADescaleStateVar
! \label{LIS_surfaceModel_DADescaleStateVar}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DADescaleStateVar(n,k)
! 
! !DESCRIPTION:
! 
!  This routine descales the state vector used in data assimilation. 
!  Scaling and descaling is done to enhance the numerical stability of 
!  matrix calculations. 
!
!EOP
    integer                :: n
    integer                :: k
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DADescaleStateVar(n,k)
          call LIS_routing_DADescaleStateVar(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DADescaleStateVar

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAUpdateState
! \label{LIS_surfaceModel_DAUpdateState}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAUpdateState(n,k)
! 
! !DESCRIPTION:
! 
!  This routine updates the model state in response to
!  analysis increments from data assimilation
!
!EOP
    integer                :: n
    integer                :: k
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAUpdateState(n,k)
          call LIS_routing_DAUpdateState(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAUpdateState

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAQCState
! \label{LIS_surfaceModel_DAQCState}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAQCState(n,k)
! 
! !DESCRIPTION:
! 
!  This routine allows the screening and masking of the model state vector
!  used in data assimilation. 
!
!EOP
    integer                :: n
    integer                :: k
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAQCState(n,k)
          call LIS_routing_DAQCState(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAQCState

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAgetStateSpaceSize
! \label{LIS_surfaceModel_DAgetStateSpaceSize}
! 
! !INTERFACE:
  function LIS_surfaceModel_DAgetStateSpaceSize(n,k) result(size)
! 
! !DESCRIPTION:
! 
!  This routine returns the size of the state vector space used 
!  in data assimilation
!
!EOP

    integer                :: n
    integer                :: k
    integer                :: size
    integer                :: m

    size = 0
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAgetStateSpaceSize(n,k,size)
          call LIS_routing_DAgetStateSpaceSize(n,k,size)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 
          
       endif
    enddo

  end function LIS_surfaceModel_DAgetStateSpaceSize

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAextractStateVector
! \label{LIS_surfaceModel_DAextractStateVector}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAextractStateVector(n,k,Nstate,state_size,stvar)

! 
! !DESCRIPTION:
! 
!  This routine extracts the state vector variables from the 
!  ESMF state object.
!
!EOP

    integer                :: n
    integer                :: k
    integer                :: Nstate
    integer                :: state_size
    real                   :: stvar(Nstate,state_size)
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAextractStateVector(n,k,state_size,stvar)
          call LIS_routing_DAextractStateVector(n,k,state_size,stvar)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo

  end subroutine LIS_surfaceModel_DAextractStateVector

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DASetFreshIncrementsStatus
! \label{LIS_surfaceModel_DASetFreshIncrementsStatus}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DASetFreshIncrementsStatus(n,k,setStatus)
    
! !DESCRIPTION: 
! 
! This routine sets the fresh increments status flag (true or false)
! based on whether assimilation occurred or not
!
!EOP

    integer                :: n
    integer                :: k
    logical                :: setStatus

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAsetFreshIncrementsStatus(n,k,setStatus)
          call LIS_routing_DAsetFreshIncrementsStatus(n,k,setStatus)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DASetFreshIncrementsStatus

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAGetFreshIncrementsStatus
! \label{LIS_surfaceModel_DAGetFreshIncrementsStatus}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAGetFreshIncrementsStatus(n,k,setStatus)
    
! !DESCRIPTION: 
! 
! This routine returns the fresh increments status flag (true or false)
!
!EOP
    integer                :: n
    integer                :: k
    logical                :: setStatus

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAgetFreshIncrementsStatus(n,k,setStatus)
          call LIS_routing_DAgetFreshIncrementsStatus(n,k,setStatus)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAGetFreshIncrementsStatus

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAsetAnlysisUpdates
! \label{LIS_surfaceModel_DAsetAnlysisUpdates}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_DAsetAnlysisUpdates(n,k,Nstate,state_size,&
       stvar,stincr)
! 
! !DESCRIPTION: 
!  This routine sets the variables the state vector
!  and state increments vector objects after assimilation
!
!EOP
    integer                :: n
    integer                :: k
    integer                :: Nstate
    integer                :: state_size
    real                   :: stvar(Nstate,state_size)
    real                   :: stincr(Nstate,state_size)

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAsetAnlysisUpdates(n,k,state_size,stvar,stincr)
          call LIS_routing_DAsetAnlysisUpdates(n,k,state_size,stvar,stincr)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo

  end subroutine LIS_surfaceModel_DAsetAnlysisUpdates

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAmapTileSpaceToObsSpace
! \label{LIS_surfaceModel_DAmapTileSpaceToObsSpace}
!
! !INTERFACE: 
  subroutine LIS_surfaceModel_DAmapTileSpaceToObsSpace(&
       n,k,tileid,st_id,en_id)
!
! !DESCRIPTION: 
!  This routine computes the mapping of the input tile space
!  index in the observation space. 
! 
!EOP
    integer                :: n 
    integer                :: k
    integer                :: tileid
    integer                :: st_id
    integer                :: en_id
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAmapTileSpaceToObsSpace(n,k,tileid,st_id,en_id)
          call LIS_routing_DAmapTileSpaceToObsSpace(n,k,tileid,st_id,en_id)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo

  end subroutine LIS_surfaceModel_DAmapTileSpaceToObsSpace

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAgetStateVarNames
! \label{LIS_surfaceModel_DAgetStateVarNames}
!
! !INTERFACE: 
  subroutine LIS_surfaceModel_DAgetStateVarNames(n,k,stateNames)
!
! !DESCRIPTION: 
! 
! This routine returns the variable names used in the state vector
!EOP
    integer                :: n 
    integer                :: k
    character(len=*)       :: stateNames(LIS_rc%nstVars(k))

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAgetStateVarNames(n,k,stateNames)
          call LIS_routing_DAgetStateVarNames(n,k,stateNames)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo

  end subroutine LIS_surfaceModel_DAgetStateVarNames

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAobsTransform
! \label{LIS_surfaceModel_DAobsTransform}
!
! !INTERFACE: 
  subroutine LIS_surfaceModel_DAobsTransform(n,k)
!
! !DESCRIPTION: 
!  This routine transforms the observations to be consistent with
!  the model state space. 
!EOP 

    integer                :: n 
    integer                :: k

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAobsTransform(n,k)
          call LIS_routing_DAobsTransform(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo

  end subroutine LIS_surfaceModel_DAobsTransform

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAmapObsToModel
! \label{LIS_surfaceModel_DAmapObsToModel}
!
! !INTERFACE: 
  subroutine LIS_surfaceModel_DAmapObsToModel(n,k)
!
! !DESCRIPTION: 
!  This routine maps the observations into the model state
!  (Typically used in direct insertion). 
!EOP 

    integer                :: n 
    integer                :: k

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAmapObsToLSM(n,k)
          call LIS_routing_DAmapObsToRouting(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAmapObsToModel

!BOP
! 
! !ROUTINE: LIS_surfaceModel_DAqcObsState
! \label{LIS_surfaceModel_DAqcObsState}
!
! !INTERFACE: 
  subroutine LIS_surfaceModel_DAqcObsState(n,k)
!
! !DESCRIPTION: 
!  This routine allows the screening and masking of the observation state
!  used in data assimilation. 
!EOP 

    integer                :: n 
    integer                :: k

    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_DAqcObsState(n,k)
          call LIS_routing_DAqcObsState(n,k)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 
          
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
  end subroutine LIS_surfaceModel_DAqcObsState

!BOP
! 
! !ROUTINE: LIS_surfaceModel_getlatlons
! \label{LIS_surfaceModel_getlatlons}
!
! !INTERFACE:
  subroutine LIS_surfaceModel_getlatlons(n,k,state_size,lats,lons)    
!
! !DESCRIPTION: 
!  This routine returns the lat/lon values corresponding to the state
!  vector
!EOP 

    integer                :: n
    integer                :: k
    integer                :: state_size
    real                   :: lats(state_size)
    real                   :: lons(state_size)
    
    integer                :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_getlatlons(n,k,state_size,lats,lons)
          call LIS_routing_getlatlons(n,k,state_size,lats,lons)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo

  end subroutine LIS_surfaceModel_getlatlons

!BOP
! 
! !ROUTINE: LIS_surfaceModel_finalize
! \label{LIS_surfaceModel_finalize}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_finalize

! !DESCRIPTION:
! This routine cleans up the surface models.
!EOP
    
    integer             :: m

    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_finalize()
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       endif
    enddo

  end subroutine LIS_surfaceModel_finalize

!BOP
! 
! !ROUTINE: LIS_surfaceModel_reset
! \label{LIS_surfaceModel_reset}
! 
! !INTERFACE:
  subroutine LIS_surfaceModel_reset

! !DESCRIPTION:
! This routine resets the surface models.
!EOP
    
    integer             :: m

    TRACE_ENTER("sf_reset")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_reset()
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
    TRACE_EXIT("sf_reset")

  end subroutine LIS_surfaceModel_reset


!BOP
! !ROUTINE: surfaceModel_setexport_noesmf
! \label{surfaceModel_setexport_noesmf}
! 
! !INTERFACE:
  subroutine surfaceModel_setexport_noesmf(n)
! !USES:    
    use LIS_LMLCMod
    use LISWRFGridCompMod, only : LISWRF_export
    use LIS_historyMod, only : LIS_tile2grid
! !ARGUMENTS: 
    integer, intent(in) :: n 
! !DESCRIPTION:
!
!EOP
    integer             :: m
    integer             :: c,r

    TRACE_ENTER("sf_setexp")
    do m=1,LIS_rc%nsf_model_types
       if(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lsm_index) then 
          call LIS_lsm_setexport(n)
       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%lake_index) then 

       elseif(LIS_rc%sf_model_type_select(m).eq.LIS_rc%glacier_index) then 

       endif
    enddo
!building grid space from surface model tiles

    call LIS_tile2grid(n,LISWRF_export(n)%avgsurft,LISWRF_export(n)%avgsurft_t)
    call LIS_tile2grid(n,LISWRF_export(n)%qh,LISWRF_export(n)%qh_t)
    call LIS_tile2grid(n,LISWRF_export(n)%eta_kinematic,LISWRF_export(n)%eta_kinematic_t)
    call LIS_tile2grid(n,LISWRF_export(n)%qle,LISWRF_export(n)%qle_t)
    call LIS_tile2grid(n,LISWRF_export(n)%qg,LISWRF_export(n)%qg_t)
    call LIS_tile2grid(n,LISWRF_export(n)%albedo,LISWRF_export(n)%albedo_t)
    call LIS_tile2grid(n,LISWRF_export(n)%znt,LISWRF_export(n)%znt_t)
    call LIS_tile2grid(n,LISWRF_export(n)%q1,LISWRF_export(n)%q1_t)
    call LIS_tile2grid(n,LISWRF_export(n)%smc1,LISWRF_export(n)%smc1_t)
    call LIS_tile2grid(n,LISWRF_export(n)%smc2,LISWRF_export(n)%smc2_t)
    call LIS_tile2grid(n,LISWRF_export(n)%smc3,LISWRF_export(n)%smc3_t)
    call LIS_tile2grid(n,LISWRF_export(n)%smc4,LISWRF_export(n)%smc4_t)
    call LIS_tile2grid(n,LISWRF_export(n)%stc1,LISWRF_export(n)%stc1_t)
    call LIS_tile2grid(n,LISWRF_export(n)%stc2,LISWRF_export(n)%stc2_t)
    call LIS_tile2grid(n,LISWRF_export(n)%stc3,LISWRF_export(n)%stc3_t)
    call LIS_tile2grid(n,LISWRF_export(n)%stc4,LISWRF_export(n)%stc4_t)
    call LIS_tile2grid(n,LISWRF_export(n)%chs2,LISWRF_export(n)%chs2_t)
    call LIS_tile2grid(n,LISWRF_export(n)%cqs2,LISWRF_export(n)%cqs2_t)
    call LIS_tile2grid(n,LISWRF_export(n)%snocvr,LISWRF_export(n)%snocvr_t)
    call LIS_tile2grid(n,LISWRF_export(n)%snow,LISWRF_export(n)%snow_t)
    call LIS_tile2grid(n,LISWRF_export(n)%snowh, LISWRF_export(n)%snowh_t)
    call LIS_tile2grid(n,LISWRF_export(n)%lispor, LISWRF_export(n)%lispor_t)

    call LIS_tile2grid(n,LISWRF_export(n)%rootmoist,LISWRF_export(n)%rootmoist_t)
    call LIS_tile2grid(n,LISWRF_export(n)%soilm,LISWRF_export(n)%soilm_t)
    call LIS_tile2grid(n,LISWRF_export(n)%qs,LISWRF_export(n)%qs_t)
    call LIS_tile2grid(n,LISWRF_export(n)%qsb,LISWRF_export(n)%qsb_t)
    call LIS_tile2grid(n,LISWRF_export(n)%cmc,LISWRF_export(n)%cmc_t)
    call LIS_tile2grid(n,LISWRF_export(n)%qsm,LISWRF_export(n)%qsm_t)
    call LIS_tile2grid(n,LISWRF_export(n)%emiss,LISWRF_export(n)%emiss_t)
    call LIS_tile2grid(n,LISWRF_export(n)%xice,LISWRF_export(n)%xice_t)
    call LIS_tile2grid(n,LISWRF_export(n)%sh2o1,LISWRF_export(n)%sh2o1_t)
    call LIS_tile2grid(n,LISWRF_export(n)%sh2o2,LISWRF_export(n)%sh2o2_t)
    call LIS_tile2grid(n,LISWRF_export(n)%sh2o3,LISWRF_export(n)%sh2o3_t)
    call LIS_tile2grid(n,LISWRF_export(n)%sh2o4,LISWRF_export(n)%sh2o4_t)
    call LIS_tile2grid(n,LISWRF_export(n)%relsmc1,LISWRF_export(n)%relsmc1_t)
    call LIS_tile2grid(n,LISWRF_export(n)%relsmc2,LISWRF_export(n)%relsmc2_t)
    call LIS_tile2grid(n,LISWRF_export(n)%relsmc3,LISWRF_export(n)%relsmc3_t)
    call LIS_tile2grid(n,LISWRF_export(n)%relsmc4,LISWRF_export(n)%relsmc4_t)
    call LIS_tile2grid(n,LISWRF_export(n)%xland,LISWRF_export(n)%xland_t)
    TRACE_EXIT("sf_setexp")

  end subroutine surfaceModel_setexport_noesmf

end module LIS_surfaceModelMod
