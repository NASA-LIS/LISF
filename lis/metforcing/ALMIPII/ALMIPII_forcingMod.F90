!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ALMIPII_forcingMod
!BOP
! !MODULE: ALMIPII_forcingMod
! 
! !DESCRIPTION: 
!
!
! !USES: 
  use ESMF

  implicit none
  
  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_ALMIPII      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ALMIPII_struc
!EOP

  !Suffixes 1 are for bilinear (rlon1,rlat1)
  !Suffixes 2 are for conservative (rlon2,rlat2)
  !Suffixes 3 are for nearest neighbor (rlon3,rlat3)
  type, public ::  ALMIPII_type_dec 
     real               :: ts
     integer            :: nc, nr
     character*100      :: dir   !ALMIPII Forcing Directory
     character*100      :: filename_prefix   !ALMIPII filename prefix
     integer            :: mi
     integer            :: ncid
     
     type(ESMF_Time)    :: startTime
     type(ESMF_TimeInterval) :: timestep    
     type(ESMF_Time)    :: time1, time2

     logical            :: startFlag
     integer            :: syr
     real, allocatable  :: metdata1(:,:) 
     real, allocatable  :: metdata2(:,:) 
  end type ALMIPII_type_dec

  type(ALMIPII_type_dec), allocatable :: ALMIPII_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_ALMIPII
! \label{init_ALMIPII} 
!
! !REVISION HISTORY: 
! 14 Oct 2010: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_ALMIPII(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_verify, LIS_endrun
    use LIS_timeMgrMod, only : LIS_update_timestep
    use map_utils,      only : ij_to_latlon

    implicit none
    integer, intent(in) :: findex
    
! 
! !DESCRIPTION: 
!
!EOP
    integer :: n
    integer :: status

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the ALMIP II forcing'
       write(LIS_logunit,*) '[ERR]  reader is not set up to run in forecast'
       write(LIS_logunit,*) '[ERR]  mode.  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(ALMIPII_struc(LIS_rc%nnest))

    call readcrd_ALMIPII()
    
    do n=1, LIS_rc%nnest
       ALMIPII_struc(n)%ts = 1800
       call LIS_update_timestep(LIS_rc, n, ALMIPII_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9 

    do n=1,LIS_rc%nnest
       allocate(ALMIPII_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(ALMIPII_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       ALMIPII_struc(n)%metdata1 = 0
       ALMIPII_struc(n)%metdata2 = 0

#if 0 
       ALMIPII_struc(n)%nc = 28 
       ALMIPII_struc(n)%nr = 25

       gridDesci(n,1) = 0
       gridDesci(n,2) = ALMIPII_struc(n)%nc
       gridDesci(n,3) = ALMIPII_struc(n)%nr
       gridDesci(n,4) = 8.975
       gridDesci(n,5) = 1.475
       gridDesci(n,6) = 128
       gridDesci(n,7) = 10.175
       gridDesci(n,8) = 2.825
       gridDesci(n,9) = 0.05
       gridDesci(n,10) = 0.05
       gridDesci(n,20) = 64
       
       ALMIPII_struc(n)%mi = ALMIPII_struc(n)%nc*ALMIPII_struc(n)%nr
#endif

       ALMIPII_struc(n)%syr = -1
       ALMIPII_struc(n)%startFlag = .true. 

       call ESMF_TimeIntervalSet(ALMIPII_struc(n)%timestep, s=1800,rc=status)
       call LIS_verify(status, 'Error in timeintervalset in ALMIPII')

    enddo
  end subroutine init_ALMIPII
end module ALMIPII_forcingMod
