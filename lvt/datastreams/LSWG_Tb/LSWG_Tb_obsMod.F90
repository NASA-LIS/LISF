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
! !MODULE: LSWG_Tb_obsMod
! \label(LSWG_Tb_obsMod)
!
! !INTERFACE:
module LSWG_Tb_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the LSWG_Tb SWE 
!  data. 
!
!  *NOTES* Currently only the 2007/12 to 2008/03 data has been tested. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Aug 2009   Sujay Kumar  Initial Specification
! 
!EOP

  PUBLIC :: LSWG_Tb_obsinit
  PUBLIC :: LSWG_Tbobs

  type, public :: lswgtbobsdec
     character*100        :: filename
     character*100        :: sname
     integer              :: nstns
     integer              :: nstates
     real                 :: udef
     integer              :: nts
     integer              :: numchannels
     integer              :: dataformat  !0 is AMSRE, etc., 1! AMSU
     type(ESMF_Clock)     :: clock
     type(ESMF_Time)      :: startTime, stopTime
     type(ESMF_TimeInterval) :: timestep
     logical                 :: start
     character*100           :: maskfile
     integer                 :: maskcol
     integer                 :: mask_option
     real                    :: cloud_pct
     real,  allocatable          :: Tb(:,:,:,:) !averaged, gridded Tb, in data channel space
     real,  allocatable          :: mask(:,:,:)
     integer, allocatable        :: data2rtm_channelmap(:)  !index into RTM channel
  end type lswgtbobsdec

  type(lswgtbobsdec), allocatable :: lswg_Tbobs(:)

contains
  
!BOP
! 
! !ROUTINE: LSWG_Tb_obsInit
! \label{LSWG_Tb_obsInit}
!
! !INTERFACE: 
  subroutine LSWG_Tb_obsinit(i)
! 
! !USES: 
#if 0 
    use LVT_coreMod, only : LVT_rc, LVT_config
    use LVT_RTMobsDataMod, only : LVT_RTMobsData, LVT_initializeRTMObsEntry
    use LVT_RTMhistDataMod, only : LVT_RTMhistData
    use LVT_timeMgrMod, only : LVT_calendar, LVT_clock
    use LVT_logMod,  only : LVT_verify, LVT_logunit,&
         LVT_getNextUnitNumber, LVT_releaseUnitNumber
#endif
    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading LSWG_Tb data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer            :: status, rc
    integer            :: ftn, k
    real*8             :: tdur
    integer            :: syr, smo, sda, shr, smn, sss
    integer            :: eyr, emo, eda, ehr, emn, ess
    integer            :: ts
    integer            :: siteid
    character*100      :: coordfile
    character*100      :: mdata

#if 0 
    call ESMF_ConfigGetAttribute(LVT_config, lswg_Tbobs%filename, &
         'LSWG Tb observation filename:',rc=status)
    call LVT_verify(status, 'LSWG Tb observation filename: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_config, lswg_Tbobs%sname, &
         'LSWG Tb satellite name:',rc=status)
    call LVT_verify(status, 'LSWG Tb satellite name: not defined')

!!$  Now moved to metadata file; vlevels is LIS-CRTM numchannels
!!$    call ESMF_ConfigGetAttribute(LVT_config, lswg_Tbobs%numchannels, &
!!$         'LSWG Tb number of channels:',rc=status)
!!$    call LVT_verify(status, 'LSWG Tb number of channels: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, lswg_Tbobs%dataformat, &
         'LSWG Tb data format:',rc=status)
    call LVT_verify(status, 'LSWG Tb data format: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, mdata, &
         'LSWG Tb metadata file:',rc=status)
    call LVT_verify(status, 'LSWG Tb metadata file: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,lswg_Tbobs%mask_option, &
         'LSWG Tb include cloud masking:',rc=status)
    call LVT_verify(status, 'LSWG Tb include cloud masking: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,lswg_Tbobs%maskfile, &
         'LSWG Tb cloud mask file:',rc=status)
    call LVT_verify(status, 'LSWG Tb cloud mask file: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,lswg_Tbobs%maskcol, &
         'LSWG Tb cloud mask column:',rc=status)
    call LVT_verify(status, 'LSWG Tb cloud mask column: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,lswg_Tbobs%cloud_pct, &
         'LSWG Tb cloud mask threshold(%):',rc=status)
    call LVT_verify(status, 'LSWG Tb cloud mask threshold(%): not defined')



    ftn=LVT_getNextUnitNumber()
    open(ftn,file=trim(mdata),status='old')
    write(LVT_logunit,*) '[INFO] Reading LSWG Tb metadata file ',trim(mdata)
    read(ftn,*)
    read(ftn,*) lswg_Tbobs%nstns, lswg_Tbobs%udef, syr, smo, &
         sda, shr, smn, eyr, emo, &
         eda, ehr, emn, ts   
    write(LVT_logunit,*) '[INFO] ',lswg_Tbobs%nstns, lswg_Tbobs%udef, syr, smo, &
         sda, shr, smn, eyr, emo, &
         eda, ehr, emn, ts    
    read(ftn,*)
    read(ftn,*) lswg_Tbobs%numchannels
print*, lswg_Tbobs%numchannels
    allocate(lswg_Tbobs%data2rtm_channelmap(lswg_Tbobs%numchannels))

    read(ftn,*)

!PUTTING BACK THE MAPPING
    write(LVT_logunit,*) '[INFO] Mapping of Tb data to LIS-CRTM Channels'
    do k=1,lswg_Tbobs%numchannels
       read(ftn,*) lswg_Tbobs%data2rtm_channelmap(k)
       write(LVT_logunit,*) '[INFO] LSWG channel index: ', k, &
            'LIS-RTM index: ',lswg_Tbobs%data2rtm_channelmap(k)
    enddo
       print *, lswg_Tbobs%data2rtm_channelmap

    call ESMF_TimeSet(lswg_Tbobs%starttime, yy=syr,&
         mm=smo, dd=sda,h=shr,m=smn,&
         calendar = LVT_calendar, rc=status)

    call ESMF_TimeSet(lswg_Tbobs%stoptime, yy=eyr,&
         mm=emo, dd=eda,h=ehr,m=emn,&
         calendar = LVT_calendar, rc=status)
    call LVT_releaseUnitNumber(ftn)

! WHAT DOES TS REFER TO?  OK, MATCHED TO NEAREST HOUR AS WILL BE MATCHED TO LIS THAT WAY
    call ESMF_TimeIntervalSet(lswg_Tbobs%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (lswg_Tbobs)')


    lswg_Tbobs%nts = nint((lswg_Tbobs%stoptime-lswg_Tbobs%starttime)/lswg_Tbobs%timestep)+1

!check last dimension
    call LVT_initializeRTMdataEntry(LVT_RTMobsData%Tb_obs,i,&
         "K",1,LVT_RTMhistData%Tb%vlevels)
    print *, LVT_RTMhistData%Tb%vlevels
!check last dimension
    allocate(lswg_Tbobs%Tb(LVT_rc%lnc, LVT_rc%lnr,lswg_Tbobs%nts,lswg_Tbobs%numchannels))

    allocate(lswg_Tbobs%mask(LVT_rc%lnc, LVT_rc%lnr,lswg_Tbobs%nts))
    lswg_Tbobs%mask = 1

    lswg_Tbobs%Tb = LVT_rc%udef

    lswg_Tbobs%start = .true. 
#endif   
  end subroutine LSWG_Tb_obsinit


end module LSWG_Tb_obsMod
