!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: AquariusL2sm_obsMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module AquariusL2sm_obsMod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: AquariusL2sm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: AquariusL2smobs
!EOP
  type, public :: aquariussmobsdec

     character*100          :: odir
     integer                :: mo
     real,    allocatable       :: smobs(:,:)
     logical                :: startmode 
  end type aquariussmobsdec

  type(aquariussmobsdec), allocatable:: AquariusL2smobs(:)

contains
  
!BOP
! 
! !ROUTINE: AquariusL2sm_obsInit
! \label{AquariusL2sm_obsInit}
! 
! !INTERFACE: 
  subroutine AquariusL2sm_obsinit()
! !USES: 
    use ESMF
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading AquariusL2 soil moisture data. 
! 
!EOP
    integer                 :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(AquariusL2smobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'Aquarius L2 soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, AquariusL2smobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'Aquarius L2 soil moisture observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       AquariusL2smobs(n)%startmode = .true. 

       allocate(AquariusL2smobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       AquariusL2smobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
    enddo
  end subroutine AquariusL2sm_obsinit
     
end module AquariusL2sm_obsMod
