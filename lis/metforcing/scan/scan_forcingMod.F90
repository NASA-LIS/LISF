!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module scan_forcingMod
!BOP
! !MODULE: scan_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various SCAN stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt scan\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 13Apr2007: Bailing Li:  Initial Specification
! 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_SCAN      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: scan_struc
!EOP

  type, public ::  scan_type_dec
     real          :: ts
     character(len=LIS_CONST_PATH_LEN) :: scandir 
     character*40  :: metadata 
     real          :: undef
     real*8        :: starttime,scantime1,scantime2
     integer       :: findtime1,findtime2,nstns
     logical       :: startRead
     integer, allocatable :: stnid(:)
     real, allocatable :: stnwt(:,:)
     real, allocatable :: stnlat(:),stnlon(:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type scan_type_dec

  type(scan_type_dec), allocatable :: scan_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_SCAN
! \label{init_SCAN}
! 
! !INTERFACE:
  subroutine init_SCAN(findex)
! !USES:
    use LIS_coreMod,only : LIS_rc
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
    integer, intent(in) :: findex

! !DESCRIPTION: 
!  This routines reads the runtime configurations for using the
!  SCAN station data. Using the metadata provided for the 
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_scan](\ref{readcrd_scan}) \newline
!     reads the runtime options specified for SCAN station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP
    integer :: i
    integer :: n 

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the USDA/SCAN forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(scan_struc(LIS_rc%nnest))
    call readcrd_scan()

    do n=1, LIS_rc%nnest
       scan_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, scan_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1
    do n=1,LIS_rc%nnest
       allocate(scan_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(scan_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       scan_struc(n)%metdata1 = 0
       scan_struc(n)%metdata2 = 0

       write(LIS_logunit,*) '--------------------------------------------------------'
       write(LIS_logunit,*) 'Opening SCAN Meta Data'
       open(78,file=scan_struc(n)%metadata,status='old')
       read(78,100) scan_struc(n)%nstns, scan_struc(n)%undef
       allocate(scan_struc(n)%stnid(scan_struc(n)%nstns))
       allocate(scan_struc(n)%stnlat(scan_struc(n)%nstns))
       allocate(scan_struc(n)%stnlon(scan_struc(n)%nstns))
       write(LIS_logunit,*) 'Number of stations: ',scan_struc(n)%nstns
       write(LIS_logunit,*) 'Undef value... ',scan_struc(n)%undef
       
       do i=1,scan_struc(n)%nstns
          read(78,101) scan_struc(n)%stnid(i),scan_struc(n)%stnlat(i),&
               scan_struc(n)%stnlon(i)
          write(LIS_logunit,*) i,scan_struc(n)%stnid(i),&
               scan_struc(n)%stnlat(i),scan_struc(n)%stnlon(i)
          if((scan_struc(n)%stnlat(i).gt.LIS_rc%gridDesc(n,7)) &
               .or.(scan_struc(n)%stnlon(i).gt.LIS_rc%gridDesc(n,8)) &
               .or.(scan_struc(n)%stnlat(i).lt.LIS_rc%gridDesc(n,4)) &
               .or.(scan_struc(n)%stnlon(i).lt.LIS_rc%gridDesc(n,5))) then
             write(LIS_logunit,*) 'Station ',scan_struc(n)%stnid(i),'(',scan_struc(n)%stnlat(i),',',&
                  scan_struc(n)%stnlon(i),')',&
                  'is not within bounds..'
             write(LIS_logunit,*) 'stopping program...'
             call LIS_endrun
          endif
       enddo
       write(LIS_logunit,*) '--------------------------------------------------------'
       close(78)
       
100    format(I3,F9.3)
101    format(I4,2(F8.2))
       
       allocate(scan_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),scan_struc(n)%nstns))
       call compute_stnwts(scan_struc(n)%nstns,LIS_rc%gridDesc,&
            scan_struc(n)%stnlat,scan_struc(n)%stnlon,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n),scan_struc(n)%stnwt)
    enddo
    
  end subroutine init_SCAN

end module scan_forcingMod
