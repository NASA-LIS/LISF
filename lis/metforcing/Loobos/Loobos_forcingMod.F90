!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module Loobos_forcingMod
!BOP
! !MODULE: Loobos_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various Loobos stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt Loobos\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 05 Oct 2010: David Mocko, Updated for Loobos test case
! 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_Loobos !defines the native resolution of 
                                       !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: Loobos_struc
!EOP

  type, public         :: Loobos_type_dec
     real                 :: ts
     character(len=LIS_CONST_PATH_LEN) :: Loobosfile
     real                 :: undef
     real*8               :: starttime,Loobostime1,Loobostime2
     integer              :: findtime1,findtime2,nstns
     logical              :: startRead
     character*9, allocatable :: stnid(:)
     real, allocatable        :: stnwt(:,:)
     real, allocatable        :: stnlat(:),stnlon(:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type Loobos_type_dec
  
  type(Loobos_type_dec), allocatable :: Loobos_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_Loobos
! \label{init_Loobos}
! 
! !INTERFACE:
  subroutine init_Loobos(findex)
! !USES:
    use LIS_coreMod,only     : LIS_rc
    use LIS_logMod, only     : LIS_logunit,LIS_endrun, &
         LIS_getNextUnitNumber,  &
         LIS_releaseUnitNumber
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    
    implicit none
    integer, intent(in) :: findex
    
! !DESCRIPTION:
!  This routines reads the runtime configurations for using the
!  Loobos station data. Using the metadata provided for the
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used.
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_Loobos](\ref{readcrd_Loobos}) \newline
!     reads the runtime options specified for Loobos station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP

    real    :: gmt
    integer :: i,n,styr,stmo,std,sth,stm,sts,tinc,doy
    character*25 :: skipline
    logical :: file_exists
    integer :: ftn

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the Loobos forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
    
    allocate(Loobos_struc(LIS_rc%nnest))
    call readcrd_Loobos()

    do n=1, LIS_rc%nnest
       Loobos_struc(n)%ts = 1800
       call LIS_update_timestep(LIS_rc, n, Loobos_struc(n)%ts)
    enddo


    LIS_rc%met_nf(findex) = 9 !number of met variables in Loobos

    ftn = LIS_getNextUnitNumber()
    do n = 1,LIS_rc%nnest

       allocate(Loobos_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(Loobos_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       Loobos_struc(n)%metdata1 = 0
       Loobos_struc(n)%metdata2 = 0

       write(LIS_logunit,*)                                          &
            '--------------------------------------------------------'
       write(LIS_logunit,*) 'Opening Loobos forcing file: ',     &
            trim(Loobos_struc(n)%Loobosfile)
       inquire(file=trim(Loobos_struc(n)%Loobosfile),&
            exist=file_exists)
       if (.not.file_exists) then
          write(LIS_logunit,*) 'Filename not found; stopping program'
          call LIS_endrun
       endif
       open(ftn,file=trim(Loobos_struc(n)%Loobosfile),status='old')
       Loobos_struc(n)%nstns = 1
       Loobos_struc(n)%undef = -999.0
       allocate(Loobos_struc(n)%stnid(Loobos_struc(n)%nstns))
       allocate(Loobos_struc(n)%stnlat(Loobos_struc(n)%nstns))
       allocate(Loobos_struc(n)%stnlon(Loobos_struc(n)%nstns))
       
       styr = 1996
       stmo = 12
       std  = 31
       sth  = 23
       stm  = 0 
       Loobos_struc(n)%stnlat(1) = 52.168
       Loobos_struc(n)%stnlon(1) = 5.744
       tinc = 30

       read(ftn,100) skipline
       read(ftn,100) skipline
       read(ftn,100) skipline
       read(ftn,100) skipline
       read(ftn,100) skipline

       Loobos_struc(n)%stnid(1) = 'Loobos'
       write(LIS_logunit,*) 'Number of stations: ',                  &
            Loobos_struc(n)%nstns
       write(LIS_logunit,*) 'Undef value: ',                         &
            Loobos_struc(n)%undef
       write(LIS_logunit,103) 'Starting time:',styr,'/',stmo,'/',    &
            std,' ',sth,':',stm
       write(LIS_logunit,*) 'Observation time interval (min) ',tinc
       sts = 0
       call LIS_date2time(Loobos_struc(n)%starttime,doy,gmt,     &
            styr,stmo,std,sth,stm,sts)
       
       do i = 1,Loobos_struc(n)%nstns
          write(LIS_logunit,*) Loobos_struc(n)%stnid(i),         &
               Loobos_struc(n)%stnlat(i),        &
               Loobos_struc(n)%stnlon(i)
          if ((Loobos_struc(n)%stnlat(i).gt.LIS_rc%gridDesc(n,7)) &
               .or.(Loobos_struc(n)%stnlon(i).gt.LIS_rc%gridDesc(n,8)) &
               .or.(Loobos_struc(n)%stnlat(i).lt.LIS_rc%gridDesc(n,4)) &
               .or.(Loobos_struc(n)%stnlon(i).lt.LIS_rc%gridDesc(n,5))) then
             write(LIS_logunit,*) 'Station ',                        &
                  Loobos_struc(n)%stnid(i),'(',  &
                  Loobos_struc(n)%stnlat(i),',', &
                  Loobos_struc(n)%stnlon(i),')', &
                  'is not within bounds..'
             write(LIS_logunit,*) 'stopping program...'
             call LIS_endrun
          endif
       enddo
       write(LIS_logunit,*)                                          &
            '--------------------------------------------------------'
       close(ftn)
       
100    format(A25,F6.2)
101    format(A25,1X,I4,4I2)
102    format(A25,I4)
103    format(1X,A14,I8,A1,I2,A1,I2,A1,I2,A1,I2)
       
       allocate(Loobos_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),&
            Loobos_struc(n)%nstns))
       call compute_stnwts(Loobos_struc(n)%nstns,LIS_rc%gridDesc,&
            Loobos_struc(n)%stnlat,Loobos_struc(n)%stnlon,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n),Loobos_struc(n)%stnwt)
    enddo
    call LIS_releaseUnitNumber(ftn)
    
  end subroutine init_Loobos

end module Loobos_forcingMod

