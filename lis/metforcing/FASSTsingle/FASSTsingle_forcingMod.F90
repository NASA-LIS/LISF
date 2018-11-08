!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module FASSTsingle_forcingMod
!BOP
! !MODULE: FASSTsingle_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various FASSTsingle stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt FASSTsingle\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 13 Apr 2007: Bailing Li, Initial Specification
! 05 Oct 2010: David Mocko, Updated for FASST single-point test case
! 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_FASSTsingle !defines the native resolution of 
  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: FASSTsingle_struc
!EOP

  type, public         :: FASSTsingle_type_dec
     character*40         :: FASSTsinglefile
     real                 :: undef
     real*8               :: starttime,FASSTsingletime1,FASSTsingletime2
     integer              :: findtime1,findtime2,nstns
     logical              :: startRead
     integer              :: icnti,ncols,pindex,met_count
     real(kind=8)         :: lat,mlong,elev,timeoffset,timstep,mflag,iheight
     real(kind=8)         :: metm(35),metm_back(35),metdata4(35)
     real, allocatable        :: stnwt(:,:)
  end type FASSTsingle_type_dec

  type(FASSTsingle_type_dec), allocatable :: FASSTsingle_struc(:)

contains

!BOP
!
! !ROUTINE: init_FASSTsingle
! \label{init_FASSTsingle}
! 
! !INTERFACE:
  subroutine init_FASSTsingle(findex)
! !USES:
    use LIS_coreMod,only     : LIS_rc
    use LIS_logMod, only     : LIS_logunit,LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time,LIS_doy2date

    implicit none
    integer, intent(in) :: findex
    integer(kind=4),parameter:: maxcol = 35

! !DESCRIPTION:
!  This routines reads the runtime configurations for using the
!  FASSTsingle station data. Using the metadata provided for the
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used.
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_FASSTsingle](\ref{readcrd_FASSTsingle}) \newline
!     reads the runtime options specified for FASSTsingle station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP

    real    :: gmt
    integer :: i,n,styr,sdoy,stmo,std,sth,stm,sts,doy
    integer :: icnti,ncols,pindex,met_count
    real(kind=8) :: lat,mlong,elev,timeoffset,timstep,mflag,iheight
    character*25 :: skipline
    logical :: file_exists

    !EMK...Added temporary arrays
    real :: tmp_stnlat(1), tmp_stnlon(1)
    real, allocatable :: tmp_stnwt(:,:)

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the FASST-single forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(FASSTsingle_struc(LIS_rc%nnest))
    call readcrd_FASSTsingle()
    do n = 1,LIS_rc%nnest
       write(LIS_logunit,*)                                          &
            '--------------------------------------------------------'
       write(LIS_logunit,*) 'Opening FASSTsingle forcing file: ',     &
            FASSTsingle_struc(n)%FASSTsinglefile
       inquire(file=FASSTsingle_struc(n)%FASSTsinglefile,exist=file_exists)
       if (.not.file_exists) then
          write(LIS_logunit,*) 'Filename not found; stopping program'
          call LIS_endrun
       endif
       open(78,file=FASSTsingle_struc(n)%FASSTsinglefile,status='old')
       FASSTsingle_struc(n)%nstns = 1
       FASSTsingle_struc(n)%undef = -999.0
       ! Read in start time here
       read(78,*) icnti,ncols,lat,mlong,elev,pindex,met_count,timeoffset, &
            timstep,mflag,iheight
       FASSTsingle_struc(n)%icnti      = icnti
       FASSTsingle_struc(n)%ncols      = ncols
       FASSTsingle_struc(n)%lat        = lat
       FASSTsingle_struc(n)%mlong      = mlong
       FASSTsingle_struc(n)%elev       = elev
       FASSTsingle_struc(n)%pindex     = pindex
       FASSTsingle_struc(n)%met_count  = met_count
       FASSTsingle_struc(n)%timeoffset = timeoffset
       FASSTsingle_struc(n)%timstep    = timstep
       FASSTsingle_struc(n)%mflag      = mflag
       FASSTsingle_struc(n)%iheight    = iheight
       read(78,100) skipline
       read(78,100) skipline
       read(78,'(4I4)') styr,sdoy,sth,stm
       call LIS_doy2date(styr,sdoy,stmo,std)
       write(LIS_logunit,*) 'Number of stations: ',                  &
            FASSTsingle_struc(n)%nstns
       write(LIS_logunit,103) 'Starting time:',styr,'/',stmo,'/',    &
            std,' ',sth,':',stm
       write(LIS_logunit,*) 'Observation time interval (hr) ',timstep
       sts = 0
       call LIS_date2time(FASSTsingle_struc(n)%starttime,doy,gmt,    &
            styr,stmo,std,sth,stm,sts)

       !         do i = 1,FASSTsingle_struc(n)%nstns
       !            if ((FASSTsingle_struc(n)%lat.gt.LIS_rc%gridDesc(n,7))     &
       !            .or.(FASSTsingle_struc(n)%mlong.gt.LIS_rc%gridDesc(n,8))   &
       !            .or.(FASSTsingle_struc(n)%lat.lt.LIS_rc%gridDesc(n,4))     &
       !            .or.(FASSTsingle_struc(n)%mlong.lt.LIS_rc%gridDesc(n,5))) then
       !               write(LIS_logunit,*) 'Station ',                        &
       !                                    FASSTsingle_struc(n)%lat,',',      &
       !                                    FASSTsingle_struc(n)%mlong,')',    &
       !                                    'is not within bounds..'
       !               write(LIS_logunit,*) 'stopping program...'
       !               call LIS_endrun
       !            endif
       !         enddo
       write(LIS_logunit,*)                                          &
            '--------------------------------------------------------'
       close(78)

100    format(A25,F6.2)
103    format(1X,A14,I8,A1,I2,A1,I2,A1,I2,A1,I2)

       allocate(FASSTsingle_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),&
            FASSTsingle_struc(n)%nstns))
       !EMK...Fixed argument type for lat and long (must be single precision)
       !call compute_stnwts(FASSTsingle_struc(n)%nstns,LIS_rc%gridDesc,&
       !     FASSTsingle_struc(n)%lat,FASSTsingle_struc(n)%mlong,  &
       !     LIS_rc%lnc(n)*LIS_rc%lnr(n),FASSTsingle_struc(n)%stnwt)
       allocate(tmp_stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),1))
       tmp_stnlat(1) = real(FASSTsingle_struc(n)%lat)
       tmp_stnlon(1) = real(FASSTsingle_struc(n)%mlong)
       call compute_stnwts(FASSTsingle_struc(n)%nstns,LIS_rc%gridDesc,&
            tmp_stnlat,tmp_stnlon,  &
            LIS_rc%lnc(n)*LIS_rc%lnr(n),tmp_stnwt)
       FASSTsingle_struc(n)%stnwt(:,1) = dble(tmp_stnwt(:,1))
    enddo

  end subroutine init_FASSTsingle

end module FASSTsingle_forcingMod

