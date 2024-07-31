!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISlsmSM_ANNdataMod
!BOP
! 
! !MODULE: LISlsmSM_ANNdataMod
! 
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation output as "observations". 
! 
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
!
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISlsmSM_ANNdatainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lsmsmANNdata
!
!EOP
  
  type, public :: lsmsmANNdatadec
     integer :: nvars
     integer :: nest    
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character(len=LDT_CONST_PATH_LEN) :: odir
  end type lsmsmANNdatadec

  type(lsmsmANNdatadec)  :: lsmsmANNdata

contains

!BOP
! !ROUTINE: LISlsmSM_ANNdatainit
! \label{LISlsmSM_ANNdatainit}
! 
! !INTERFACE: 
  subroutine LISlsmSM_ANNdatainit()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS:     
!    integer,     intent(IN) :: i   ! index of the observation type
! 
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model output (from a LIS simulation) as observations.  
! 
!EOP
    real                    :: run_dd(8)
    character*20            :: stime
    integer                 :: n    
    integer                 :: rc 
    real                    :: ts
    character*3             :: fnest

    n = 1

    call ESMF_ConfigGetAttribute(LDT_config,stime, &
         label="LIS soil moisture output timestep:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output timestep: not defined')
    call LDT_parseTimeString(stime, ts)
   
    call LDT_update_timestep(LDT_rc, n, ts)
    
    write(fnest,'(i3.3)') n
    call LDT_registerAlarm("LIS LSM SM alarm "//trim(fnest),&
         ts,ts)

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmANNdata%format, &
         label="LIS soil moisture output format:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmANNdata%wopt, &
         label="LIS soil moisture output methodology:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmANNdata%wstyle, &
         label="LIS soil moisture output naming style:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmANNdata%map_proj, &
         label="LIS soil moisture output map projection:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmANNdata%nest, &
         label="LIS soil moisture output nest index:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,lsmsmANNdata%odir, &
         label="LIS soil moisture output directory:",rc=rc)
    call LDT_verify(rc,'LIS soil moisture output directory: not defined')

    if(lsmsmANNdata%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(6),rc=rc)       
       
    elseif(lsmsmANNdata%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS soil moisture domain y-dimension size: not defined')

    elseif(lsmsmANNdata%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"LIS soil moisture domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,run_dd(8),rc=rc)

    endif

  end subroutine LISlsmSM_ANNdatainit
  
end module LISlsmSM_ANNdataMod
