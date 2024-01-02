!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: MOD10A1_ANNdataMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module MOD10A1_ANNdataMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MOD10A1_ANNdatainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MOD10A1obs
!EOP
  type, public :: MOD10A1obsdec

     character(len=LDT_CONST_PATH_LEN) :: odir
     integer                    :: nc,nr
     integer, allocatable       :: n11(:)
     real                       :: gridDesc(50)
     type(proj_info)            :: mod_proj
  end type MOD10A1obsdec

  type(MOD10A1obsdec), allocatable:: MOD10A1obs(:)

contains
  
!BOP
! 
! !ROUTINE: MOD10A1_ANNdatainit
! \label{MOD10A1_ANNdatainit}
! 
! !INTERFACE: 
  subroutine MOD10A1_ANNdatainit()
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

    allocate(MOD10A1obs(LDT_rc%nnest))

    ts = 86400.0
    call LDT_update_timestep(LDT_rc, n, ts)

    write(fnest,'(i3.3)') n
    call LDT_registerAlarm("MOD10A1 data alarm "//trim(fnest),&
         ts,ts)

    call ESMF_ConfigFindLabel(LDT_config, &
         'MOD10A1 data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, MOD10A1obs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'MOD10A1 data directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       MOD10A1obs(n)%gridDesc = 0
       
       cornerlat1 = max(-59.995,nint((LDT_rc%gridDesc(n,4)+59.995)/0.01)*0.01-59.995-2*0.01)
       cornerlon1 = max(-179.995,nint((LDT_rc%gridDesc(n,5)+179.995)/0.01)*0.01-179.995-2*0.01)
       cornerlat2 = min(89.995,nint((LDT_rc%gridDesc(n,7)+59.995)/0.01)*0.01-59.995+2*0.01)
       cornerlon2 = min(179.995,nint((LDT_rc%gridDesc(n,8)+179.995)/0.01)*0.01-179.995+2*0.01)
       
       MOD10A1obs(n)%nr = nint((cornerlat2-cornerlat1)/0.01)+1
       MOD10A1obs(n)%nc = nint((cornerlon2-cornerlon1)/0.01)+1
       
       allocate(MOD10A1obs(n)%n11(MOD10A1obs(n)%nc*MOD10A1obs(n)%nr))
       
       !filling the items needed by the interpolation library
       MOD10A1obs(n)%gridDesc(1) = 0  !input is EASE grid
       MOD10A1obs(n)%gridDesc(2) = MOD10A1obs(n)%nc
       MOD10A1obs(n)%gridDesc(3) = MOD10A1obs(n)%nr
       MOD10A1obs(n)%gridDesc(4) = cornerlat1
       MOD10A1obs(n)%gridDesc(5) = cornerlon1
       MOD10A1obs(n)%gridDesc(7) = cornerlat2
       MOD10A1obs(n)%gridDesc(8) = cornerlon2
       MOD10A1obs(n)%gridDesc(6) = 128
       MOD10A1obs(n)%gridDesc(9) = 0.01
       MOD10A1obs(n)%gridDesc(10) = 0.01
       MOD10A1obs(n)%gridDesc(20) = 64
       
       call upscaleByAveraging_input(MOD10A1obs(n)%gridDesc,&
            LDT_rc%gridDesc(n,:),MOD10A1obs(n)%nc*MOD10A1obs(n)%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),MOD10A1obs(n)%n11)
       
       call map_set(PROJ_LATLON, MOD10A1obs(n)%gridDesc(4), &
            MOD10A1obs(n)%gridDesc(5), 0.0,&
            MOD10A1obs(n)%gridDesc(9), MOD10A1obs(n)%gridDesc(10),0.0,&
            MOD10A1obs(n)%nc,MOD10A1obs(n)%nr,MOD10A1obs(n)%mod_proj)

    enddo
  end subroutine MOD10A1_ANNdatainit
     
end module MOD10A1_ANNdataMod
