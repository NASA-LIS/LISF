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
!  10 Dec 2020: Hiroko Beaudoing; Incorporate crop calendar and concurrent
!                                 irrigation schemes
!                                 Made irrig_type_dec public (was private)
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
  public :: LIS_irrig_state      !data structure containing irrigation states
  public :: LIS_irrig_struc      !data structure containing irrigation variables
!EOP  

  type, public :: irrig_type_dec
     real               :: outInterval
     character*100      :: models_used
     logical            :: stats_file_open
     character*50       :: cropcalendar
     integer            :: cropseasons
     real               :: sprinkler_start  !sprinkler start time
     real               :: sprinkler_duration   !sprinkler duration
     real               :: sprinkler_thresh  !sprinkler threshold
     real               :: sprinkler_efcor   !sprinkler efficiency
     real               :: drip_start  !drip start time
     real               :: drip_duration   !drip duration
     real               :: drip_thresh  !drip threshold
     real               :: drip_efcor   !drip efficiency
     real               :: flood_start  !flood start time
     real               :: flood_duration   !flood duration
     real               :: flood_thresh  !flood threshold
     real               :: flood_efcor   !flood efficency
     real,allocatable   :: plantDay(:,:)
     real,allocatable   :: harvestDay(:,:)
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
       if( LIS_rc%irrigation_type == "Flood" .or. &
           LIS_rc%irrigation_type == "Concurrent" ) then
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

<<<<<<< HEAD
     ! HKB--added irrigation type specific configurations 
     ! Set trigger check start time [local hour] and duration in lis.config
       call ESMF_ConfigFindLabel(LIS_config,"Irrigation Sprinkler start time:",&
            rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%sprinkler_start,rc=rc)
          call LIS_verify(rc,"Irrigation Sprinkler start time: not defined")
          write(LIS_logunit,*) "[INFO] Irrigation Sprinkler start at :  ",&
                             LIS_irrig_struc(n)%sprinkler_start
       enddo

       call ESMF_ConfigFindLabel(LIS_config,"Irrigation Sprinkler duration:",&
             rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%sprinkler_duration,rc=rc)
          call LIS_verify(rc,"Irrigation Sprinkler duration: not defined")
          write(LIS_logunit,*) "[INFO] for [hrs] :  ",LIS_irrig_struc(n)%sprinkler_duration
       enddo

       call ESMF_ConfigFindLabel(LIS_config,"Irrigation Drip start time:", &
            rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%drip_start,rc=rc)
          call LIS_verify(rc,"Irrigation Drip start time: not defined")
          write(LIS_logunit,*) "[INFO] Irrigation Drip start at :  ",&
                             LIS_irrig_struc(n)%drip_start
       enddo

       call ESMF_ConfigFindLabel(LIS_config,"Irrigation Drip duration:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%drip_duration,rc=rc)
          call LIS_verify(rc,"Irrigation Drip duration: not defined")
          write(LIS_logunit,*) "[INFO] for [hrs] :  ",LIS_irrig_struc(n)%drip_duration
       enddo
                             
       call ESMF_ConfigFindLabel(LIS_config,"Irrigation Flood start time:", &
            rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%flood_start,rc=rc)
          call LIS_verify(rc,"Irrigation Flood start time: not defined")
          write(LIS_logunit,*) "[INFO] Irrigation Flood start at :  ",&
                             LIS_irrig_struc(n)%flood_start
       enddo

       call ESMF_ConfigFindLabel(LIS_config,"Irrigation Flood duration:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%flood_duration,rc=rc)
          call LIS_verify(rc,"Irrigation Flood duration: not defined")
          write(LIS_logunit,*) "[INFO] for [hrs] :  ",LIS_irrig_struc(n)%flood_duration
       enddo

     ! Threshold for which irrigation is triggered per irrigation type
       call ESMF_ConfigFindLabel(LIS_config, &
            "Irrigation threshold for Sprinkler:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%sprinkler_thresh,rc=rc)
          call LIS_verify(rc,"Irrigation threshold for Sprinkler: not defined")
          write(LIS_logunit,*) "[INFO] and irrigation thresholds for Sprinkler: ",&
                             LIS_irrig_struc(n)%sprinkler_thresh
       enddo
       call ESMF_ConfigFindLabel(LIS_config, &
            "Irrigation threshold for Drip:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%drip_thresh,rc=rc)
          call LIS_verify(rc,"Irrigation threshold for Drip: not defined")
          write(LIS_logunit,*) "[INFO] for Drip:  ",LIS_irrig_struc(n)%drip_thresh
       enddo
       call ESMF_ConfigFindLabel(LIS_config, &
            "Irrigation threshold for Flood:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%flood_thresh,rc=rc)
          call LIS_verify(rc,"Irrigation threshold for Flood: not defined")
          write(LIS_logunit,*) "[INFO] for Flood:  ",LIS_irrig_struc(n)%flood_thresh
       enddo
                            
     ! Irrigation efficiency correction per irrigation type
       call ESMF_ConfigFindLabel(LIS_config, &
            "Irrigation Sprinkler efficiency:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%sprinkler_efcor,rc=rc)
          call LIS_verify(rc,"Irrigation Sprinkler efficiency: not defined")
       enddo
       call ESMF_ConfigFindLabel(LIS_config, &
            "Irrigation Drip efficiency:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%drip_efcor,rc=rc)
          call LIS_verify(rc,"Irrigation Drip efficiency: not defined")
       enddo
       call ESMF_ConfigFindLabel(LIS_config, &
            "Irrigation Flood efficiency:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%flood_efcor,rc=rc)
          call LIS_verify(rc,"Irrigation Flood efficiency: not defined")
       enddo

     ! Crop calendar options
       call ESMF_ConfigFindLabel(LIS_config, &
            "Crop Calendar use:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%cropcalendar,default="none",rc=rc)
          call LIS_verify(rc,"Crop Calender use: option not specified in the config file")
       enddo
       call ESMF_ConfigFindLabel(LIS_config,"Crop seasons:",rc=rc)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config, &
               LIS_irrig_struc(n)%cropseasons,default=1,rc=rc)
          call LIS_verify(rc,"Crop Calender use: option not specified in the config file")
       enddo
       do n=1,LIS_rc%nnest
        if ( LIS_irrig_struc(n)%cropcalendar .ne. "none" ) then
           allocate(LIS_irrig_struc(n)%plantDay( &
              LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_irrig_struc(n)%cropseasons))
           allocate(LIS_irrig_struc(n)%harvestDay( &
              LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_irrig_struc(n)%cropseasons))
           LIS_irrig_struc(n)%plantDay = 0.0
           LIS_irrig_struc(n)%harvestDay = 0.0
        endif
       enddo
         
=======
!------Wanshu----irrigation scheduling based on DVEG On--------
     LIS_rc%irrigation_dveg  = 0 ! Default is no
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_dveg,&
            label="Irrigation scheduling based on dynamic vegetation:",default=0,rc=rc)
       call LIS_verify(rc,"Irrigation scheduling based on dynamic vegetation: not defined")
       write(LIS_logunit,*) "[INFO] Irrigation scheduling based on dynamic vegetation:  ",&
                             LIS_rc%irrigation_dveg
!------------------------------------------------------

!------Wanshu---GW abstraction based on irrigation groundwater ratio data--------
     LIS_rc%irrigation_SourcePartition  = 0 ! Default is no
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_SourcePartition,&
            label="Irrigation source water partition:",default=0,rc=rc)
       call LIS_verify(rc,"Irrigation source water partition: not defined")
       write(LIS_logunit,*) "[INFO] Irrigation source water partition:  ",&
                             LIS_rc%irrigation_SourcePartition
!------------------------------------------------------

>>>>>>> master
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
                if (LIS_irrig_struc(n)%stats_file_open) then
                   call LIS_create_stats_filename(n,statsfile,"IRRIGATION")
                   LIS_irrig_struc(n)%stats_file_open = .false.
                   open_stats = .true.
                endif
             endif

             call LIS_create_output_filename(n,outfile,&
                  model_name ="IRRIGATION")

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
