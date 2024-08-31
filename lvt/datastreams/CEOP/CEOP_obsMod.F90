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
module CEOP_obsMod

  use ESMF

  implicit none

  PRIVATE

  PUBLIC :: CEOP_obsInit
  PUBLIC :: CEOPobs
  
  type, public :: ceopstn
      real, allocatable    :: tskin(:)
      integer :: col
      integer :: row
  end type ceopstn

  type, public :: ceopobsdec
     character*50            :: odir
     integer                 :: nstns
     integer                 :: readsfc
     integer                 :: readflx
     integer                 :: readstm
     logical                 :: startFlag
     type(ESMF_TimeInterval) :: ts
     character*50,  allocatable  :: campaign(:)
     character*50,  allocatable  :: stnname(:)
     character*50,  allocatable  :: locname(:)
     type(ESMF_Time)         :: startTime
     type(ceopstn), allocatable  :: stn(:)
  end type ceopobsdec


  type(ceopobsdec), allocatable :: CEOPobs(:)

contains

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
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
  subroutine CEOP_obsInit(i)

    use LVT_coreMod,    only : LVT_rc, LVT_config
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod,     only : LVT_verify, LVT_logunit, &
         LVT_getNextUnitNumber, LVT_releaseUnitNumber
    
    implicit none
    
    integer,     intent(IN) :: i   ! index of the observation type
    integer                 :: rc,status
    character*100           :: stnlist
    type(ESMF_Time)         :: startTime, stopTime
    integer                 :: ftn
    integer                 :: nts
    integer                 :: k,j

    if(.not.allocated(ceopobs)) then 
       allocate(ceopobs(LVT_rc%nDataStreams))
    endif

    write(LVT_logunit,*) '[INFO] Initializing CEOP data reader....'
    ceopobs(i)%startFlag = .true. 
    call ESMF_ConfigGetAttribute(LVT_config, ceopobs(i)%odir, &
         label='CEOP data directory: ',rc=rc)
    call LVT_verify(rc, 'CEOP data directory: not defined')
    call ESMF_ConfigGetAttribute(LVT_config, stnlist, &
         label='CEOP station list file: ',rc=rc)
    call LVT_verify(rc, 'CEOP station list file: not defined')

    ftn=LVT_getNextUnitNumber()
    open(ftn,file=trim(stnlist),status='old')
    read(ftn,*) 
    read(ftn,*) ceopobs(i)%nstns
    read(ftn,*) 

    allocate(ceopobs(i)%campaign(ceopobs(i)%nstns))
    allocate(ceopobs(i)%stnname(ceopobs(i)%nstns))
    allocate(ceopobs(i)%locname(ceopobs(i)%nstns))

    do k=1,ceopobs(i)%nstns
       read(ftn,*) ceopobs(i)%campaign(k), ceopobs(i)%stnname(k),&
            ceopobs(i)%locname(k)
       write(LVT_logunit,*) '[INFO] ',ceopobs(i)%campaign(k), ceopobs(i)%stnname(k),&
            ceopobs(i)%locname(k)
    enddo
    
    call LVT_releaseUnitNumber(ftn)

    call ESMF_ConfigGetAttribute(LVT_config, ceopobs(i)%readsfc,&
         label='CEOP read surface meteorology data:',rc=rc)
    call LVT_verify(rc,'CEOP read surface meteorology data: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_config, ceopobs(i)%readflx,&
         label='CEOP read flux data:',rc=rc)
    call LVT_verify(rc,'CEOP read flux data: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, ceopobs(i)%readstm,&
         label='CEOP read soil moisture and temperature data:',rc=rc)
    call LVT_verify(rc,'CEOP read soil moisture and temperature data : not defined')

         
    call ESMF_TimeIntervalSet(ceopobs(i)%ts,&
         s = 1800,rc=status)
    call LVT_verify(status)

    call ESMF_TimeSet(startTime, yy = LVT_rc%syr, &
         mm = LVT_rc%smo, &
         dd = LVT_rc%sda, &
         h  = LVT_rc%shr, &
         m  = LVT_rc%smn, & 
         s  = LVT_rc%sss, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status)

    call ESMF_TimeSet(stopTime, yy = LVT_rc%eyr, &
         mm = LVT_rc%emo, &
         dd = LVT_rc%eda, &
         h  = LVT_rc%ehr, &
         m  = LVT_rc%emn, & 
         s  = LVT_rc%ess, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status)

    nts = nint((stopTime-startTime)/ceopobs(i)%ts)+1

!-------------------------------------------------------------------------
!  CEOP data contains soil moisture, soil temperature (6 layers) 
!  psurf, tair, qair, wind, swdown, lwdown, rainf, avgsurft, snowdepth
!  qle, qh, qg
!-------------------------------------------------------------------------
    if(ceopobs(i)%readsfc.eq.1) then 
       allocate(ceopobs(i)%stn(ceopobs(i)%nstns))
       do k=1,ceopobs(i)%nstns
          allocate(ceopobs(i)%stn(k)%tskin(nts))
          ceopobs(i)%stn(k)%tskin = LVT_rc%udef
!          allocate(ceopobs(i)%stn(k)%col)
!          allocate(ceopobs(i)%stn(k)%row)
          ceopobs(i)%stn(k)%col = -1
          ceopobs(i)%stn(k)%row = -1
       enddo
    endif

    write(LVT_logunit,*) '[INFO] Finished initializing CEOP data reader....'  
  end subroutine CEOP_obsInit


end module CEOP_obsMod
