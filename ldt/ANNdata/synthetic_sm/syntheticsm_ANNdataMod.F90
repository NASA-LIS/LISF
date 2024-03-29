!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: syntheticsm_ANNdataMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module syntheticsm_ANNdataMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: syntheticsm_ANNdatainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: syntheticsmobs
!EOP
  type, public :: syntheticsmobsdec

     character(len=LDT_CONST_PATH_LEN) :: odir
     integer                :: mo
     real,    allocatable       :: smobs(:,:)
     logical                :: startmode 
  end type syntheticsmobsdec

  type(syntheticsmobsdec), allocatable:: syntheticsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: syntheticsm_ANNdatainit
! \label{syntheticsm_ANNdatainit}
! 
! !INTERFACE: 
  subroutine syntheticsm_ANNdatainit()
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

    n = 1

    allocate(syntheticsmobs(LDT_rc%nnest))

    call ESMF_ConfigGetAttribute(LDT_config,stime, &
         label="Synthetic soil moisture observation timestep:",rc=rc)
    call LDT_verify(rc,&
         'Synthetic soil moisture observation timestep: not defined')
    call LDT_parseTimeString(stime, ts)

    call LDT_update_timestep(LDT_rc, n, ts)

    write(fnest,'(i3.3)') n
    call LDT_registerAlarm("Synthetic sm alarm "//trim(fnest),&
         ts,ts)

    call ESMF_ConfigFindLabel(LDT_config, &
         'Synthetic soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, syntheticsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'Synthetic soil moisture observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       syntheticsmobs(n)%startmode = .true. 

       allocate(syntheticsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       syntheticsmobs(n)%smobs = -9999.0
    
    enddo
  end subroutine syntheticsm_ANNdatainit
     
end module syntheticsm_ANNdataMod
