!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: MODISlst_ANNdataMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module MODISlst_ANNdataMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MODISlst_ANNdatainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MODISlstobs
!EOP
  type, public :: MODISlstobsdec

     character(len=LDT_CONST_PATH_LEN) :: odir
     integer                    :: nc,nr
     integer, allocatable       :: n11(:)
     real                       :: gridDesc(50)
     
  end type MODISlstobsdec

  type(MODISlstobsdec), allocatable:: MODISlstobs(:)

contains
  
!BOP
! 
! !ROUTINE: MODISlst_ANNdatainit
! \label{MODISlst_ANNdatainit}
! 
! !INTERFACE: 
  subroutine MODISlst_ANNdatainit()
! !USES: 
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading SYNTHETIC soil moisture data. 
! 
!EOP
    integer                 :: npts
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 
    real                    :: ts
    character*3             :: fnest
    character*20            :: stime
    real                  :: cornerlat1, cornerlat2
    real                  :: cornerlon1, cornerlon2
    n = 1

    allocate(MODISlstobs(LDT_rc%nnest))

    ts = 86400.0
    call LDT_update_timestep(LDT_rc, n, ts)

    write(fnest,'(i3.3)') n
    call LDT_registerAlarm("MODIS LST data alarm "//trim(fnest),&
         ts,ts)

    call ESMF_ConfigFindLabel(LDT_config, &
         'MODIS LST data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, MODISlstobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'MODIS LST data directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       MODISlstobs(n)%gridDesc = 0
       
       cornerlat1 = max(-59.995,nint((LDT_rc%gridDesc(n,4)+59.995)/0.01)*0.01-59.995-2*0.01)
       cornerlon1 = max(-179.995,nint((LDT_rc%gridDesc(n,5)+179.995)/0.01)*0.01-179.995-2*0.01)
       cornerlat2 = min(89.995,nint((LDT_rc%gridDesc(n,7)+59.995)/0.01)*0.01-59.995+2*0.01)
       cornerlon2 = min(179.995,nint((LDT_rc%gridDesc(n,8)+179.995)/0.01)*0.01-179.995+2*0.01)
       
       MODISlstobs(n)%nr = nint((cornerlat2-cornerlat1)/0.01)+1
       MODISlstobs(n)%nc = nint((cornerlon2-cornerlon1)/0.01)+1
       
       allocate(MODISlstobs(n)%n11(MODISlstobs(n)%nc*MODISlstobs(n)%nr))
       
       !filling the items needed by the interpolation library
       MODISlstobs(n)%gridDesc(1) = 0  !input is EASE grid
       MODISlstobs(n)%gridDesc(2) = MODISlstobs(n)%nc
       MODISlstobs(n)%gridDesc(3) = MODISlstobs(n)%nr
       MODISlstobs(n)%gridDesc(4) = cornerlat1
       MODISlstobs(n)%gridDesc(5) = cornerlon1
       MODISlstobs(n)%gridDesc(7) = cornerlat2
       MODISlstobs(n)%gridDesc(8) = cornerlon2
       MODISlstobs(n)%gridDesc(6) = 128
       MODISlstobs(n)%gridDesc(9) = 0.01
       MODISlstobs(n)%gridDesc(10) = 0.01
       MODISlstobs(n)%gridDesc(20) = 64
       
       call upscaleByAveraging_input(MODISlstobs(n)%gridDesc,&
            LDT_rc%gridDesc(n,:),MODISlstobs(n)%nc*MODISlstobs(n)%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),MODISlstobs(n)%n11)
       
    enddo
  end subroutine MODISlst_ANNdatainit
     
end module MODISlst_ANNdataMod
