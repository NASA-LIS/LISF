!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module Bondville_forcingMod
!BOP
! !MODULE: Bondville_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various Bondville stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt Bondville\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 05 Oct 2010: David Mocko, Updated for Bondville test case
! 26 Oct 2018: David Mocko, Updated for Noah-MP-4.0.1 HRLDAS test case
! 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_Bondville !defines the native resolution of 
                                       !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: Bondville_struc
!EOP

  type, public         :: Bondville_type_dec
     real                 :: ts
     character*80         :: Bondvillefile
     integer              :: mp
     real                 :: undef
     real*8               :: starttime,Bondvilletime1,Bondvilletime2
     integer              :: findtime1,findtime2,nstns
     logical              :: startRead
     character*9, allocatable :: stnid(:)
     real, allocatable        :: stnwt(:,:)
     real, allocatable        :: stnlat(:),stnlon(:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type Bondville_type_dec
  
  type(Bondville_type_dec), allocatable :: Bondville_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_Bondville
! \label{init_Bondville}
! 
! !INTERFACE:
  subroutine init_Bondville(findex)
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
!  Bondville station data. Using the metadata provided for the
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used.
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_Bondville](\ref{readcrd_Bondville}) \newline
!     reads the runtime options specified for Bondville station data
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
       write(LIS_logunit,*) '[ERR] Currently the Bondville forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
    
    allocate(Bondville_struc(LIS_rc%nnest))
    call readcrd_Bondville()

    do n=1, LIS_rc%nnest
       Bondville_struc(n)%ts = 1800
       call LIS_update_timestep(LIS_rc, n, Bondville_struc(n)%ts)
    enddo


    LIS_rc%met_nf(findex) = 9 !number of met variables in Bondville

    ftn = LIS_getNextUnitNumber()
    do n = 1,LIS_rc%nnest

       allocate(Bondville_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(Bondville_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       Bondville_struc(n)%metdata1 = 0
       Bondville_struc(n)%metdata2 = 0

       write(LIS_logunit,*)                                          &
            '--------------------------------------------------------'
       write(LIS_logunit,*) 'Opening Bondville forcing file: ',     &
            trim(Bondville_struc(n)%Bondvillefile)
       inquire(file=trim(Bondville_struc(n)%Bondvillefile),&
            exist=file_exists)
       if (.not.file_exists) then
          write(LIS_logunit,*) 'Filename not found; stopping program'
          call LIS_endrun
       endif
       Bondville_struc(n)%nstns = 1
       Bondville_struc(n)%undef = -999.0
       allocate(Bondville_struc(n)%stnid(Bondville_struc(n)%nstns))
       allocate(Bondville_struc(n)%stnlat(Bondville_struc(n)%nstns))
       allocate(Bondville_struc(n)%stnlon(Bondville_struc(n)%nstns))
       open(ftn,file=trim(Bondville_struc(n)%Bondvillefile),status='old')

! If MP=0, use the older "bondville.dat" forcing file.
!    This option is to support older testcases.
! If MP.ne.0, use the new Noah-MP-4.0.1 HRLDAS "bondville.dat" file.
       if (Bondville_struc(n)%MP.eq.0) then
          read(ftn,100) skipline
! Read in start time here
          read(ftn,101) skipline,styr,stmo,std,sth,stm
          read(ftn,100) skipline
          read(ftn,100) skipline
          read(ftn,100) skipline,Bondville_struc(n)%stnlat(1)
          read(ftn,100) skipline,Bondville_struc(n)%stnlon(1)
          read(ftn,102) skipline,tinc
       else
          do i = 1,4
             read(ftn,100) skipline
          enddo
          read(ftn,110) skipline,Bondville_struc(n)%stnlat(1)
          read(ftn,110) skipline,Bondville_struc(n)%stnlon(1)
          tinc = 1800
! Read in start time here
          do i = 1,49
             read(ftn,101) skipline
          enddo
          read(ftn,111) styr,stmo,std,sth,stm
       endif

       tinc = tinc / 60
       Bondville_struc(n)%stnid(1) = 'Bondville'
       write(LIS_logunit,*) 'Number of stations: ',                  &
            Bondville_struc(n)%nstns
       write(LIS_logunit,*) 'Undef value: ',                         &
            Bondville_struc(n)%undef
       write(LIS_logunit,103) 'Starting time:',styr,'/',stmo,'/',    &
            std,' ',sth,':',stm
       write(LIS_logunit,*) 'Observation time interval (min) ',tinc
       sts = 0
       call LIS_date2time(Bondville_struc(n)%starttime,doy,gmt,     &
            styr,stmo,std,sth,stm,sts)
       
       do i = 1,Bondville_struc(n)%nstns
          write(LIS_logunit,*) Bondville_struc(n)%stnid(i),         &
               Bondville_struc(n)%stnlat(i),        &
               Bondville_struc(n)%stnlon(i)
          if ((Bondville_struc(n)%stnlat(i).gt.LIS_rc%gridDesc(n,7)) &
               .or.(Bondville_struc(n)%stnlon(i).gt.LIS_rc%gridDesc(n,8)) &
               .or.(Bondville_struc(n)%stnlat(i).lt.LIS_rc%gridDesc(n,4)) &
               .or.(Bondville_struc(n)%stnlon(i).lt.LIS_rc%gridDesc(n,5))) then
             write(LIS_logunit,*) 'Station ',                        &
                  Bondville_struc(n)%stnid(i),'(',  &
                  Bondville_struc(n)%stnlat(i),',', &
                  Bondville_struc(n)%stnlon(i),')', &
                  'is not within bounds..'
             write(LIS_logunit,*) 'stopping program...'
             call LIS_endrun
          endif
       enddo
       write(LIS_logunit,*)                                          &
            '--------------------------------------------------------'
       close(ftn)
       
 100   format(A25,F6.2)
 101   format(A25,1X,I4,4I2)
 102   format(A25,I4)
 103   format(1X,A14,I8,A1,I2,A1,I2,A1,I2,A1,I2)

 110   format(A25,1X,F6.2)
 111   format(I4,4I3)

       allocate(Bondville_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),&
            Bondville_struc(n)%nstns))
       call compute_stnwts(Bondville_struc(n)%nstns,LIS_rc%gridDesc,&
            Bondville_struc(n)%stnlat,Bondville_struc(n)%stnlon,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n),Bondville_struc(n)%stnwt)
    enddo
    call LIS_releaseUnitNumber(ftn)
    
  end subroutine init_Bondville

end module Bondville_forcingMod

