!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_irrigationMod
!BOP
!
! !MODULE: LIS_irrigationMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  11 Nov 2012: Sujay Kumar; Initial implementation
!  29 May 2019; Jessica Erlingis; Incorporate Wanshu Nie's max/min GVF update
!
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod

  implicit none
  
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_irrigation_init
  public :: LIS_irrigation_run
  public :: LIS_irrigation_output
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_irrig_state      !data structure containing irrigation variables
!EOP  

  type, private :: irrig_type_dec
     real               :: outInterval
     character*100      :: models_used
     logical            :: stats_file_open
  end type irrig_type_dec

  type(irrig_type_dec),allocatable :: LIS_irrig_struc(:)

  type(ESMF_State),    allocatable :: LIS_irrig_state(:)

contains

!BOP
! 
! !ROUTINE: LIS_irrigation_init
! \label{LIS_irrigation_init}
! 
! !DESCRIPTION:
!
! Allocates memory for data structures used for reading 
! irrigation datasets. The irrigationdepth field is updated by the external
! files. The irrigation water equivalent fields are expected to be set
! by the model. 
! 
! !INTERFACE:
  subroutine LIS_irrigation_init

! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   use LIS_timeMgrMod,   only : LIS_registerAlarm, LIS_parseTimeString
!EOP
    integer       :: n
    integer       :: status
    integer       :: rc
    integer       :: ios, nid
    character*100 :: temp
    character*10  :: time
    character*1   :: nestid(2)
    logical       :: file_exists
! ___________________________________________________

 !- Read in Config file irrigation inputs:

  ! Read in type of irrigation scheme selected (spray,flood,drip):
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_type,&
         label="Irrigation scheme:",default="none",rc=rc)
    call LIS_verify(rc,&
         'Irrigation scheme: option not specified in the config file')

    if( LIS_rc%irrigation_type .ne. "none" ) then 

       write(LIS_logunit,*) "[INFO] Irrigation scheme selected:  ",&
                             trim(LIS_rc%irrigation_type)
 
       allocate(LIS_irrig_state(LIS_rc%nnest))
       allocate(LIS_irrig_struc(LIS_rc%nnest))

     ! Frequency with which irrigation field is written out:
       call ESMF_ConfigGetAttribute(LIS_config,time,&
            label="Irrigation output interval:",rc=rc)
       call LIS_verify(rc,"Irrigation output interval: not defined")
       write(LIS_logunit,*) "[INFO] Irrigation output interval:  ",time

     ! Threshold for which irrigation is triggered, like for flood irrigation:
       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_thresh,&
            label="Irrigation threshold:",rc=rc)
       call LIS_verify(rc,"Irrigation threshold: not defined")
       write(LIS_logunit,*) "[INFO] and irrigation threshold:  ",&
                             LIS_rc%irrigation_thresh

     ! Parameters to control the GVF threshold based on the range of GVF
     ! (shdmax-shdmin) for which sprinkler irrigation is triggered:(WN)
       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_GVFparam1,&
            label="Irrigation GVF parameter 1:",rc=rc)
       call LIS_verify(rc,"Irrigation GVF parameter 1: not defined")
       write(LIS_logunit,*) "and irrigation GVF parameter 1:  ",&
                             LIS_rc%irrigation_GVFparam1

       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_GVFparam2,&
            label="Irrigation GVF parameter 2:",rc=rc)
       call LIS_verify(rc,"Irrigation GVF parameter 2: not defined")
       write(LIS_logunit,*) "and irrigation GVF parameter 2:  ",&
                             LIS_rc%irrigation_GVFparam2

     ! Max. soil layer depth for irrigation to reach to (available for flood only):
       LIS_rc%irrigation_mxsoildpth = 1
       if( LIS_rc%irrigation_type == "Flood" ) then
          call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_mxsoildpth,&
               label="Irrigation max soil layer depth:", default=1, rc=rc)
          call LIS_verify(rc,"Irrigation max soil layer depth: not defined")
          write(LIS_logunit,*) "[INFO]and irrigation max soil depth:  ",&
                                LIS_rc%irrigation_mxsoildpth
       endif

     ! JE Remove irrigated water from groundwater
       LIS_rc%irrigation_GWabstraction = 0 ! Default is no
       ! Need to add model sanity check here to make sure model contains GW (?)
       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_GWabstraction,&
            label="Groundwater abstraction for irrigation:",default=0,rc=rc)
       call LIS_verify(rc,"Groundwater abstraction for irrigation: not defined")
       write(LIS_logunit,*) "[INFO]and irrigation withdrawn from GW:  ",&
                             LIS_rc%irrigation_GWabstraction

     ! Register irrigation output interval:
       do n=1,LIS_rc%nnest
          call LIS_parseTimeString(time,LIS_irrig_struc(n)%outInterval)
          call LIS_registerAlarm("LIS irrigation output interval",&
               real(LIS_irrig_struc(n)%outInterval), &
               LIS_irrig_struc(n)%outInterval)
          LIS_irrig_struc(n)%models_used = trim(LIS_rc%irrigation_type)
          LIS_irrig_struc(n)%stats_file_open = .true.
       enddo

       do n=1,LIS_rc%nnest
          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid
          
          LIS_irrig_state(n) = ESMF_StateCreate(name="LSM Irrigation State"//&
               nestid(1)//nestid(2), rc=status)
          call LIS_verify(status, &
               "ESMF_StateCreate failed in LIS_irrigation_init")
       enddo

    !- Initiate the irrigation scheme selected in lis.config file:
       call irrigationschemeinit(trim(LIS_rc%irrigation_type)//char(0),&
            LIS_irrig_state)

    endif

  end subroutine LIS_irrigation_init

!BOP
! 
! !ROUTINE: LIS_irrigation_run
! \label{LIS_irrigation_run}
! 
! !INTERFACE:
  subroutine LIS_irrigation_run(n)
! !USES: 
    implicit none

! !ARGUMENTS: 
    integer  :: n 

! !DESCRIPTION:
! This routine runs the specified irrigation model.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    if(LIS_rc%irrigation_type.ne."none") then     

       call getirrigationlsmstates(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%irrigation_type)//char(0), n,LIS_irrig_state(n))
       call applyirrigationupdates(trim(LIS_rc%irrigation_type)//char(0),&
            n,LIS_irrig_state(n))
       
    endif

  end subroutine LIS_irrigation_run

!BOP
! 
! !ROUTINE: LIS_irrigation_output
! \label{LIS_irrigation_output}
! 
! !INTERFACE:
  subroutine LIS_irrigation_output(n)
! !USES: 

    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_historyMod, only : LIS_writeModelOutput
    use LIS_fileIOMod,  only : LIS_create_output_directory, &
         LIS_create_output_filename,  &
         LIS_create_stats_filename

! !ARGUMENTS: 
    integer, intent(in)   :: n 

! !DESCRIPTION:
! This routine writes the irrigation model output.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP
    
    logical           :: alarmCheck,open_stats
    character*100     :: outfile, statsfile

    if(LIS_rc%irrigation_type.ne."none") then 
       alarmCheck = LIS_isAlarmRinging(LIS_rc,&
            "LIS irrigation output interval")
       if(alarmCheck) then 
          open_stats = .false. 
          if(LIS_rc%wopt.ne."none") then 
             if(LIS_masterproc) then 
                call LIS_create_output_directory('IRRIGATION')
                call LIS_create_output_filename(n,outfile,&
                     model_name ="IRRIGATION")
                if(LIS_irrig_struc(n)%stats_file_open) then 
                   call LIS_create_stats_filename(n,statsfile,"IRRIGATION")
                   LIS_irrig_struc(n)%stats_file_open = .false. 
                   open_stats = .true. 
                endif
             endif
             call LIS_writeModelOutput(n,outfile,statsfile,              &
                  open_stats,outInterval=LIS_irrig_struc(n)%outInterval, &
                  nsoillayers=1, lyrthk = (/1.0/),                       &
                  nsoillayers2=1,                                        &
                  model_name=LIS_irrig_struc(n)%models_used,group=4)
          endif
       endif
    endif
    
  end subroutine LIS_irrigation_output

end module LIS_irrigationMod
