!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: GHCN_ANNdataMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module GHCN_ANNdataMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GHCN_ANNdatainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GHCNobs
!EOP
  type, public :: GHCNobsdec

     character(len=LDT_CONST_PATH_LEN) :: odir
     integer                    :: nstns
     real                       :: udef
     integer                    :: nts
     character*100, allocatable :: stnid(:)
     real,          allocatable :: stnlat(:)
     real,          allocatable :: stnlon(:)
     real,          allocatable :: stnelev(:)
     type(ESMF_Clock)           :: clock
     type(ESMF_Time)            :: startTime
     type(ESMF_TimeInterval)    :: timestep
     real,  allocatable         :: snod(:,:)
     logical              :: startflag
     integer              :: yr,syr
     integer              :: mo,smo
  end type GHCNobsdec

  type(GHCNobsdec), allocatable:: GHCNobs(:)

contains
  
!BOP
! 
! !ROUTINE: GHCN_ANNdatainit
! \label{GHCN_ANNdatainit}
! 
! !INTERFACE: 
  subroutine GHCN_ANNdatainit()
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
    integer            :: n 
    integer            :: status, rc
    integer            :: ftn, k
    real*8             :: tdur
    integer            :: syr, smo, sda, shr, smn, sss
    integer            :: eyr, emo, eda, ehr, emn, ess
    integer            :: ts
    character*50       :: stnid
    real               :: stnlat, stnlon,stnelev
    character(len=LDT_CONST_PATH_LEN) :: stationfile
    integer            :: c,r
    real               :: col,row

    n = 1

    allocate(GHCNobs(LDT_rc%nnest))

    call LDT_update_timestep(LDT_rc, n, 86400.0)

    call ESMF_ConfigFindLabel(LDT_config, &
         'GHCN data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, GHCNobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'GHCN data directory: not defined')
    enddo
    call ESMF_ConfigFindLabel(LDT_config, &
         'GHCN station file:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, stationfile, &
            rc=status)
       call LDT_verify(status, &
            'GHCN station file: not defined')
    enddo

    do n=1,LDT_rc%nnest
       call ESMF_TimeIntervalSet(GHCNobs(n)%timestep,&
            s=86400,rc=status)
       call LDT_verify(status, 'Error in setting timestep (GHCN_ANNdataMod)')
       ftn = LDT_getNextUnitNumber()
       open(ftn,file=trim(stationfile), status='old')
       write(LDT_logunit,*) '[INFO]  Reading GHCN station file ',trim(stationfile)
! first read to figure out the number of stations included in the 
! LVT/LIS domain
       status = 0 
       do while(status.eq.0) 
          read(ftn,*,iostat=status) stnid, stnlat, stnlon,stnelev
          if(status.eq.0) then 
             call latlon_to_ij(LDT_domain(n)%ldtproj, stnlat, stnlon, &
                  col, row)
             c = nint(col)
             r = nint(row)
             if(c.ge.1.and.c.le.LDT_rc%lnc(n).and.&
                  r.ge.1.and.r.le.LDT_rc%lnr(n)) then 
                GHCNobs(n)%nstns =  GHCNobs(n)%nstns +1
             endif
          else
             exit
          endif
       enddo
       call LDT_releaseUnitNumber(ftn)

       allocate(GHCNobs(n)%stnid(GHCNobs(n)%nstns))
       allocate(GHCNobs(n)%stnlat(GHCNobs(n)%nstns))
       allocate(GHCNobs(n)%stnlon(GHCNobs(n)%nstns))
       allocate(GHCNobs(n)%stnelev(GHCNobs(n)%nstns))
       
       ftn=LDT_getNextUnitNumber()
       open(ftn,file=trim(stationfile),status='old')
       
       status = 0 
       k = 0 
       do while(status.eq.0) 
          read(ftn,*,iostat=status) stnid, stnlat, stnlon,stnelev
          if(status.eq.0) then 
             call latlon_to_ij(LDT_domain(n)%ldtproj, stnlat, stnlon, &
                  col, row)
             c = nint(col)
             r = nint(row)
             if(c.ge.1.and.c.le.LDT_rc%lnc(n).and.&
                  r.ge.1.and.r.le.LDT_rc%lnr(n)) then 
                k = k + 1
                GHCNobs(n)%stnid(k) = stnid
                GHCNobs(n)%stnlat(k) = stnlat
                GHCNobs(n)%stnlon(k) = stnlon
                GHCNobs(n)%stnelev(k) = stnelev
             endif
          else
             exit
          endif
       enddo
       call LDT_releaseUnitNumber(ftn)
       
       GHCNobs(n)%nts = 366
       
       allocate(GHCNobs(n)%snod(GHCNobs(n)%nstns, GHCNobs(n)%nts))
       GHCNobs(n)%snod = LDT_rc%udef
       
       GHCNobs(n)%yr = -1
       GHCNobs(n)%mo = -1
       
    enddo
  end subroutine GHCN_ANNdatainit
     
end module GHCN_ANNdataMod
