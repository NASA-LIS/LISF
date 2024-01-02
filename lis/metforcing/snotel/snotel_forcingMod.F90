!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module snotel_forcingMod
!BOP
! !MODULE: snotel_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various SNOTEL stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt snotel\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 08Jun2010: Yuqiong Liu:  Initial Specification
! 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_SNOTEL      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: snotel_struc
!EOP

  type, public ::  snotel_type_dec
     real           :: ts
     character(len=LIS_CONST_PATH_LEN) :: snoteldir 
     character*100  :: metadata 
     character*100  :: coorddata 
     real          :: undef
     real*8        :: starttime,snoteltime1,snoteltime2
     integer       :: findtime1,findtime2,nstns
     logical       :: startRead
     character*2, allocatable :: statename(:)
     character*6, allocatable :: stnid(:)
     real, allocatable :: stnlat(:),stnlon(:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type snotel_type_dec

  type(snotel_type_dec), allocatable :: snotel_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_SNOTEL
! \label{init_SNOTEL}
! 
! !INTERFACE:
  subroutine init_SNOTEL(findex)
! !USES:
    use LIS_coreMod,only : LIS_rc
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
    integer, intent(in) :: findex

! !DESCRIPTION: 
!  This routines reads the runtime configurations for using the
!  SNOTEL station data. Using the metadata provided for the 
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_snotel](\ref{readcrd_snotel}) \newline
!     reads the runtime options specified for SNOTEL station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP
    integer :: i,j
    integer :: n 
    integer :: siteid
    logical :: file_exists
    character*30 :: stnname
    character*6 :: stnid
    character*2 :: statename
    real    :: stnlat, stnlon
    integer, allocatable :: idx(:)  !whether the site is in local domain
    integer     :: nstns_total !total number of snotel sites to be processed

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the SNOTEL forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
    allocate(snotel_struc(LIS_rc%nnest))
    call readcrd_snotel()

    do n=1, LIS_rc%nnest
       snotel_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, snotel_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1

    do n=1,LIS_rc%nnest
       allocate(snotel_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(snotel_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       snotel_struc(n)%metdata1 = 0
       snotel_struc(n)%metdata2 = 0

       write(LIS_logunit,*) '--------------------------------------------------------'
       inquire(file=snotel_struc(n)%metadata,exist=file_exists)
       if (.not.file_exists) then
            write(LIS_logunit,*) 'file does not exist: ', snotel_struc(n)%metadata
            call LIS_endrun            
       endif 
       write(LIS_logunit,*) 'Opening SNOTEL Meta Data'
       open(78,file=snotel_struc(n)%metadata,status='old')
       read(78,*)   !skip the header
       read(78,100) nstns_total, snotel_struc(n)%undef
       write(LIS_logunit,*) 'Total number of stations: ',nstns_total
       write(LIS_logunit,*) 'Undef value... ',snotel_struc(n)%undef
       close(78)

       inquire(file=snotel_struc(n)%coorddata,exist=file_exists)
       if (.not.file_exists) then
            write(LIS_logunit,*) 'file does not exist: ', snotel_struc(n)%coorddata
            call LIS_endrun
       endif
       open(78,file=snotel_struc(n)%coorddata,form='formatted')
       write(LIS_logunit,*) 'opening coord file'

! first check how many stations are located in the local domain
       allocate(idx(nstns_total))
       snotel_struc(n)%nstns = 0;
       do i=1,nstns_total
          read(78,101) statename, stnname,stnid,siteid,stnlat, stnlon
          if((stnlat.le.LIS_rc%gridDesc(n,7)) &
               .and.(stnlon.le.LIS_rc%gridDesc(n,8)) &
               .and.(stnlat.ge.LIS_rc%gridDesc(n,4)) &
               .and.(stnlon.ge.LIS_rc%gridDesc(n,5))) then
              idx(i) = 1;
              snotel_struc(n)%nstns =  snotel_struc(n)%nstns + 1
          else
              idx(i) = 0;
          endif
       enddo
       close(78)
       write(LIS_logunit,*) 'number of snotel sites in local domain: ', snotel_struc(n)%nstns
       write(LIS_logunit,*) '--------------------------------------------------------'

       allocate(snotel_struc(n)%statename(snotel_struc(n)%nstns))
       allocate(snotel_struc(n)%stnid(snotel_struc(n)%nstns))
       allocate(snotel_struc(n)%stnlat(snotel_struc(n)%nstns))
       allocate(snotel_struc(n)%stnlon(snotel_struc(n)%nstns))

! open coord file again to read the coordinates
       open(78,file=snotel_struc(n)%coorddata,form='formatted')
       write(LIS_logunit,*) 'opening coord file'
       j=0
       do i=1,nstns_total
          read(78,101) statename, stnname,stnid,siteid,stnlat, stnlon
          if(idx(i) .eq. 1) then
              j=j+1
              snotel_struc(n)%statename(j) = statename
              snotel_struc(n)%stnid(j) = stnid
              snotel_struc(n)%stnlat(j) = stnlat
              snotel_struc(n)%stnlon(j) = stnlon
              
              write(LIS_logunit,*) i,snotel_struc(n)%stnid(j),snotel_struc(n)%stnlat(j), snotel_struc(n)%stnlon(j)
          endif
       enddo

       deallocate(idx)
       close(78)
       
100    format(I2,F8.1)
101    format(A2,A30,A6,I14,2F10.3)
       
!       allocate(snotel_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),snotel_struc(n)%nstns))
!       call compute_stnwts(snotel_struc(n)%nstns,LIS_rc%gridDesc,&
!            snotel_struc(n)%stnlat,snotel_struc(n)%stnlon,&
!            LIS_rc%lnc(n)*LIS_rc%lnr(n),snotel_struc(n)%stnwt)
    enddo
    
  end subroutine init_SNOTEL

end module snotel_forcingMod
