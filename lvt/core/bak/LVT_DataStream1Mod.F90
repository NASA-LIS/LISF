!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: LVT_DataStream1Mod
! \label(LVT_DataStream1Mod)
!
! !INTERFACE:
module LVT_DataStream1Mod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for handling LIS model output data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
! !USES:       
  implicit none


  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_DataStream1Init
  public :: LVT_readDataStream1
!  public :: LVT_readDAobsData

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP  
contains

!BOP
! 
! !ROUTINE: LVT_DataStream1Init
! \label{LVT_DataStream1Init}
!
! !INTERFACE:   
  subroutine LVT_DataStream1Init
! 
! !USES: 
    use ESMF
    use LVT_coreMod, only : LVT_rc, LVT_config
    use LVT_histDataMod, only : LVT_histDataInit, LVT_histData
    use LVT_logMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer         :: ftn
    integer         :: c,r,rc

!options: "pixel-by-pixel" (default), "region-based"
    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%sp_avg_mode, &
         label="Spatial averaging mode:",&
         rc=rc)
    call LVT_verify(rc,'Spatial averaging mode: not defined')
    
    if(LVT_rc%sp_avg_mode.eq."region-based") then 
       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%reg_maskfile, &
            label="Regional mask file for spatial averaging:",&
            rc=rc)
       call LVT_verify(rc,'Regional mask file for spatial averaging: not defined')
       
       allocate(LVT_rc%regmask(LVT_rc%lnc,LVT_rc%lnr))
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=LVT_rc%reg_maskfile,form='unformatted',access='direct',&
            recl=LVT_rc%lnc*LVT_rc%lnr*4)
       read(ftn,rec=1) LVT_rc%regmask
       call LVT_releaseUnitNumber(ftn)
       
       LVT_rc%regmask_max = -1000.0
       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             if(LVT_rc%regmask(c,r).gt.LVT_rc%regmask_max) then 
                LVT_rc%regmask_max = LVT_rc%regmask(c,r)
             endif
          enddo
       enddo

    endif

    if(LVT_rc%computeEnsMetrics.eq.1) then 
       LVT_rc%npts = LVT_rc%lis_ntiles       
    else
       LVT_rc%npts = LVT_rc%ngrid
    endif

    ! if the LIS output is structured in tile space, then the arrays are 
    ! sized to be ntiles. Else the arrays are sized to be ngrid. If the ouput is 
    ! written as 2d grid, then then it will be stored in the 1-d land space. 

    call LVT_histDataInit(LVT_rc%npts)

  end subroutine LVT_DataStream1Init

!BOP
! 
! !ROUTINE: LVT_readDataStream1
! \label{LVT_readDataStream1}
!
! !INTERFACE: 
  subroutine LVT_readDataStream1
! 
! !USES: 
    use LVT_coreMod,   only : LVT_rc
    use LVT_fileIOMod, only : LVT_create_output_filename
    use LVT_logMod,    only : LVT_logunit
    use LVT_historyMod, only : LVT_readModelOutput
    use LVT_statsDataMod,only : LVT_stats

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine reads the LIS model data and selects the specified
!   set of variables to be analyzed. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    character*200 :: fname
    logical          :: file_exists
!hardcoded for now

    if(LVT_rc%chkTS) then 
       if(LVT_rc%output_src.eq."LSM") then 
          call LVT_create_output_filename(LVT_rc%nnest, fname, &
               'SURFACEMODEL',&
               writeint = LVT_rc%ts, wout = LVT_rc%lis_out_format,&
               style=LVT_rc%lis_out_style)
       elseif(LVT_rc%output_src.eq."Routing") then 
          call LVT_create_output_filename(LVT_rc%nnest, fname, &
               'ROUTING',&
               writeint = LVT_rc%ts, wout = LVT_rc%lis_out_format,&
               style=LVT_rc%lis_out_style)
       elseif(LVT_rc%output_src.eq."RTM") then 
          call LVT_create_output_filename(LVT_rc%nnest, fname, &
               'RTM',&
               writeint = LVT_rc%ts, wout = LVT_rc%lis_out_format,&
               style=LVT_rc%lis_out_style)
       elseif(LVT_rc%output_src.eq."Irrigation") then 
          call LVT_create_output_filename(LVT_rc%nnest, fname, &
               'IRRIGATION',&
               writeint = LVT_rc%ts, wout = LVT_rc%lis_out_format,&
               style=LVT_rc%lis_out_style)
       endif
       inquire(file=trim(fname),exist=file_exists)
       
       if(file_exists) then 
          ! The following line is modified by Shugong Wang: LVT should be LIS
          !write(LVT_logunit,*) 'Reading LVT output ',trim(fname)
          write(LVT_logunit,*) 'Reading LIS output ',trim(fname)
          call LVT_readModelOutput(trim(fname),LVT_rc%lis_out_format)
       else
          LVT_stats%datamask = 0.0  !initially set all to false. 
!toggle these print/stop statements on/off
!if LVT needs to work with missing LIS files
          write(LVT_logunit,*) 'LIS file ',trim(fname),' does not exist'
!          write(LVT_logunit,*) 'Program stopping.. '
!          stop
       endif
    endif

  end subroutine LVT_readDataStream1

#if 0 
!BOP
! 
! !ROUTINE: LVT_readDAobsData
! \label{LVT_readDAobsData}
!
! !INTERFACE: 
  subroutine LVT_readDAobsData
! 
! !USES: 
    use LVT_coreMod,   only : LVT_rc
    use LVT_logMod,    only : LVT_logunit
    use LVT_fileIOMod, only : LVT_create_daobs_filename
    use LVT_historyMod, only : LVT_readModelOutput

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads the processed (interpolated, qc'd and masked) data
!  assimilation observations from LIS. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    character*100 :: fname 
    logical          :: file_exists

    call LVT_create_daobs_filename(LVT_rc%nnest, fname)
    
    inquire(file=trim(fname),exist=file_exists)
    
    if(file_exists) then 
       write(LVT_logunit,*) 'Reading ',trim(fname)
       call LVT_readModelOutput(trim(fname),LVT_rc%lis_out_format)       
    endif
    
  end subroutine LVT_readDAobsData

#endif
end module LVT_DataStream1Mod
