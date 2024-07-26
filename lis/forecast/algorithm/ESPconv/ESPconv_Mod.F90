!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module ESPconv_Mod
!BOP
!
! !MODULE: ESPconv_Mod
!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
! 25 Oct 2016: Sujay Kumar, initial specification
!
! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_forecastMod
  use LIS_numerRecipesMod
  use LIS_ran2_gasdev

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ESPconv_initialize
  public :: ESPconv_sampledate
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  type, public :: espconv_type_dec
     type(ESMF_TimeInterval), allocatable :: twInterval(:)
     integer                              :: iterId
     integer, allocatable                 :: stime(:,:)
     integer, allocatable                 :: etime(:,:)
     integer, allocatable                 :: sampleCounter(:)
     integer, allocatable                 :: include_fcstyr(:)
     integer                              :: ref_da
     integer                              :: ref_yr
     type(ESMF_TimeInterval)              :: offSet
     logical                              :: sampleFlag
  end type espconv_type_dec
  
  type(espconv_type_dec), allocatable :: ESPconv_struc(:)

!EOP
  SAVE

contains
!BOP
! 
! !ROUTINE: ESPconv_initialize
! \label{ESPconv_initialize}
! 
! !INTERFACE:
  subroutine ESPconv_initialize()
! !USES:

! !DESCRIPTION:        
!
!EOP
    implicit none
    
    integer :: n,i,k,t,rc
    character*10             :: time_val
    type(ESMF_Time)          :: startTime, endTime
    type(ESMF_TimeInterval)  :: deltaT
    real                     :: twInterval
    integer, allocatable     :: stime(:,:),etime(:,:)
    integer, allocatable     :: nIterations(:)
! ______________

    write(LIS_logunit,*) "[INFO] ESP conventional forecast method selected "

    allocate(ESPconv_struc(LIS_rc%nnest))

    allocate(LIS_forecast_struc(1)%seed(LIS_rc%nmetforc, NRANDSEED))
    allocate(ESPconv_struc(1)%twInterval(LIS_rc%nmetforc))
    allocate(ESPconv_struc(1)%sampleCounter(LIS_rc%nmetforc))
    allocate(ESPconv_struc(1)%stime(LIS_rc%nmetforc,3))
    allocate(ESPconv_struc(1)%etime(LIS_rc%nmetforc,3))
!    allocate(ESPconv_struc(1)%include_fcstyr(LIS_rc%nmetforc))
    do n=1,LIS_rc%nnest
      allocate(ESPconv_struc(n)%include_fcstyr(LIS_rc%nmetforc))
    enddo

    allocate(stime(LIS_rc%nmetforc, 3))
    allocate(etime(LIS_rc%nmetforc, 3))
    allocate(nIterations(LIS_rc%nmetforc))

    ESPconv_struc%iterId = -1
    ESPconv_struc%sampleFlag =.false.

    ! Starting from a forecast coldstart file
    if(LIS_forecast_struc(1)%startMode.eq."coldstart") then 
       LIS_forecast_struc(1)%seed = -100
       call init_randseed( LIS_forecast_struc(1)%seed )
    else   ! Starting from a forecast restart file
       call LIS_forecast_readrestart()
    endif
    
    ! Read in lis.config file entries: 
    ! Start forecast ESP year:
    call ESMF_ConfigFindLabel(LIS_config,"ESP conventional start time of the forcing archive:",rc=rc)
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc
          do i=1,3
             call ESMF_ConfigGetAttribute(LIS_config,stime(k,i),rc=rc)
             call LIS_verify(rc,"ESP conventional start time of the forcing archive: not defined")
          enddo
       enddo
    enddo
    ! Final forecast ESP year:
    call ESMF_ConfigFindLabel(LIS_config,"ESP conventional end time of the forcing archive:",rc=rc)
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc
          do i=1,3
             call ESMF_ConfigGetAttribute(LIS_config,etime(k,i),rc=rc)
             call LIS_verify(rc,"ESP conventional end time of the forcing archive: not defined")
          enddo
       enddo
    enddo
    ! Include targeted forecast year:
    call ESMF_ConfigFindLabel(LIS_config,"ESP conventional include targeted forecast year:",rc=rc)
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc
          call ESMF_ConfigGetAttribute(LIS_config,ESPconv_struc(n)%include_fcstyr(k),rc=rc)
          call LIS_verify(rc,"ESP conventional include targeted forecast year: not defined")
       enddo
    enddo
    
    ! Set ESMF-based Start and End times (as set in lis.config file):
    do n=1,LIS_rc%nnest
       do k=1,LIS_rc%nmetforc   ! number for forcing datasets
          call ESMF_TimeSet(startTime, yy = stime(k,1), &
               mm = stime(k,2), &
               dd = stime(k,3), &
               h  = 0, &
               m  = 0, & 
               s  = 0, &
               calendar = LIS_calendar, & 
               rc = rc)
          call LIS_verify(rc,'error in ESMF_TimeSet:startTime in ESPconv_initialize')
          
          call ESMF_TimeSet(endTime, yy = etime(k,1), &
               mm = etime(k,2), &
               dd = etime(k,3), &
               h  = 0, &
               m  = 0, &
               s  = 0, &
               calendar = LIS_calendar, & 
               rc = rc)
          call LIS_verify(rc,'error in ESMF_TimeSet:endTime in ESPconv_initialize')
          
          call ESMF_TimeIntervalSet(deltaT,s=86400*365,rc=rc)
          call LIS_verify(rc,'error in ESMF_TimeIntervalSet in ESPconv_initialize')

          
          ! Modify so the nIterations matches the number of years, designated by user:
          if( stime(k,2) == etime(k,2) ) then  ! If start == end month
             nIterations(k) = etime(k,1) - stime(k,1) + 1
          else    !  Months don't match
             nIterations(k) = etime(k,1) - stime(k,1) 
          endif

          ! Account for targeted forecast year, based on user input:
          if(LIS_rc%yr.le.etime(k,1).and.&
               LIS_rc%yr.ge.stime(k,1)) then ! current year is in the sample

             ! Do NOT include targeted forecast year:
             if( ESPconv_struc(n)%include_fcstyr(k) == 0 ) then
                niterations(k) = nIterations(k)-1
             endif
          endif
          if((LIS_rc%yr < stime(k,1) .or.  &
              LIS_rc%yr > etime(k,1)).and. &
              ESPconv_struc(n)%include_fcstyr(k) == 1 ) then

             write(LIS_logunit,*) "[WARN] Since the target forecast year,",LIS_rc%yr
             write(LIS_logunit,*) "[WARN] is outside of the archive data period selected,"
             write(LIS_logunit,*) "[WARN] the 'include forecast year' option is not activated."
          endif 

          write(LIS_logunit,*) "[INFO] ESP conventional forecast selected years: "
          write(LIS_logunit,*) "[INFO] ",stime(k,1),",",etime(k,1)

          write(LIS_logunit,*) "[INFO] ESPconv forcing dataset, number of members:",&
                k,nIterations(k)

          ESPconv_struc(n)%stime(k,:) = stime(k,:)
          ESPconv_struc(n)%etime(k,:) = etime(k,:)
          ESPconv_struc(n)%sampleCounter(k) = 0          
       enddo
    enddo

    write(LIS_logunit,*) "[WARN] Current ESP conventional setup is only based"
    write(LIS_logunit,*) "[WARN]  on start and end years, and it does not take"
    write(LIS_logunit,*) "[WARN]  into account month / day, at this time " 
    write(LIS_logunit,*) "[WARN]  (unless starting / ending months are different)."
    write(LIS_logunit,*) " "
    
    ! Initialize total number of forecast members for all forcing datasets:
    do n=1,LIS_rc%nnest
       LIS_forecast_struc(n)%nIterations = 0 
       do k=1,LIS_rc%nmetforc
          if(nIterations(k).gt.LIS_forecast_struc(n)%nIterations) then 
             LIS_forecast_struc(n)%nIterations = nIterations(k)
          endif
       enddo
    enddo

  end subroutine ESPconv_initialize


  subroutine ESPconv_sampledate(n, kk, k, yr, mo, da)
    
    integer     :: n     ! nest
    integer     :: kk    ! iteration id
    integer     :: k     ! forcing dataset (findex)
    integer     :: yr    ! target fcst runtime year
    integer     :: mo    ! target fcst runtime month
    integer     :: da    ! target fcst runtime day

    integer                 :: da1
    integer                 :: twInterval
    type(ESMF_Time)         :: stTime, currTime, currTime1
    type(ESMF_TimeInterval) :: deltaT
    integer                 :: rc
! _______________________________________

    ! Advance forecast year:
    yr = ESPconv_struc(n)%stime(k,1) + & 
         (kk-1)

    ! Do NOT include targeted forecast year:
    if( ESPconv_struc(n)%include_fcstyr(k) == 0 ) then
      if( yr >= LIS_rc%yr ) then  ! Withholding target forecast year
         if( ESPconv_struc(n)%stime(k,1) <= LIS_rc%yr ) then
           yr = yr + 1
         endif
      endif
    endif

    write(LIS_logunit,*) "[INFO] ESP conventional forecast member, year:",&
                          kk, yr

    ! If current date is 29 and sampled year is not a leap year:
    if(da.eq.29.and.&
         .not.((mod(yr,4) .eq. 0 .and. mod(yr, 100).ne.0) &
         .or.(mod(yr,400) .eq.0))) then 
       da = 28
    endif

!    if( mo==2 .and. da==29 )then
!      if((mod(yr,4) .eq. 0 .and. mod(yr, 100).ne.0) &
!         .or.(mod(yr,400) .eq.0)) then 
!        da=29
!      else
!        da = 28
!      endif
!    endif

         
  end subroutine ESPconv_sampledate

end module ESPconv_Mod

