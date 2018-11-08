!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module coop_forcingMod
!BOP
! !MODULE: coop_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various COOP stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt coop\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 13Jul2011: Yuqiong Liu:  Initial Specification
! 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_COOP      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: coop_struc
!EOP

  type, public ::  coop_type_dec
     real           :: ts
     character*100  :: coopdir 
     character*100  :: metadata 
     character*100  :: coorddata 
     real          :: undef
     real*8        :: starttime,cooptime1,cooptime2
     integer       :: findtime1,findtime2,nstns
     integer       :: nstates  !number of states where coop sites are located 
     logical       :: startRead
     character*2, allocatable :: statename(:)
     integer, allocatable :: stnid(:)
     real, allocatable :: stnlat(:),stnlon(:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type coop_type_dec

  type(coop_type_dec), allocatable :: coop_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_COOP
! \label{init_COOP}
! 
! !INTERFACE:
  subroutine init_COOP(findex)
! !USES:
    use LIS_coreMod,only : LIS_rc
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
    integer, intent(in) :: findex

! !DESCRIPTION: 
!  This routines reads the runtime configurations for using the
!  COOP station data. Using the metadata provided for the 
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_coop](\ref{readcrd_coop}) \newline
!     reads the runtime options specified for COOP station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP
    integer :: n,i,j
    integer :: siteid
    logical :: file_exists
    character*39 :: sitename
    integer, allocatable :: idx(:) !whether the site is in local domain
    real    :: lat1, lon1
    integer :: nstns_total !total number of coop sites to be processed

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the CO-OP forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(coop_struc(LIS_rc%nnest))
    call readcrd_coop()

    do n=1, LIS_rc%nnest
       coop_struc(n)%ts = 3600 
       call LIS_update_timestep(LIS_rc, n, coop_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1
    do n=1,LIS_rc%nnest

       allocate(coop_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(coop_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       coop_struc(n)%metdata1 = 0
       coop_struc(n)%metdata2 = 0

       write(LIS_logunit,*) '--------------------------------------------------------'
       inquire(file=coop_struc(n)%metadata,exist=file_exists)
       if (.not.file_exists) then
            write(LIS_logunit,*) 'file does not exist: ', coop_struc(n)%metadata
            call LIS_endrun            
       endif 
       write(LIS_logunit,*) 'Opening COOP Meta Data'
       open(78,file=coop_struc(n)%metadata,status='old')
       read(78,*)   !skip the header
       read(78,100) nstns_total, coop_struc(n)%undef
       read(78,*)
       read(78,*) coop_struc(n)%nstates
       read(78,*)
       allocate(coop_struc(n)%statename(coop_struc(n)%nstates))
       allocate(idx(nstns_total))
       do i=1,coop_struc(n)%nstates
           read(78,fmt='(A2)') coop_struc(n)%statename(i)
       enddo
       write(LIS_logunit,*) 'Total number of stations: ',nstns_total
       write(LIS_logunit,*) 'Undef value... ',coop_struc(n)%undef
       close(78)

       inquire(file=coop_struc(n)%coorddata,exist=file_exists)
       if (.not.file_exists) then
            write(LIS_logunit,*) 'file does not exist: ', coop_struc(n)%coorddata
            call LIS_endrun
       endif

       open(78,file=coop_struc(n)%coorddata,form='formatted')
       write(LIS_logunit,*) 'opening coord file'

! first check how many sites are located in the local domain
       coop_struc(n)%nstns = 0;
       do i=1,nstns_total
          read(78,101) siteid, sitename, lat1, lon1
          if((lat1.le.LIS_rc%gridDesc(n,7)) &
               .and.(lon1.le.LIS_rc%gridDesc(n,8)) &
               .and.(lat1.ge.LIS_rc%gridDesc(n,4)) &
               .and.(lon1.ge.LIS_rc%gridDesc(n,5))) then
              idx(i) = 1;
              coop_struc(n)%nstns =  coop_struc(n)%nstns + 1
          else
              idx(i) = 0;
          endif
       enddo
       close(78)
       write(LIS_logunit,*) 'number of coop sites in local domain: ', coop_struc(n)%nstns  
       write(LIS_logunit,*) '--------------------------------------------------------'
           
 !      allocate(coop_struc(n)%stnname(coop_struc(n)%nstns))
       allocate(coop_struc(n)%stnid(coop_struc(n)%nstns))
       allocate(coop_struc(n)%stnlat(coop_struc(n)%nstns))
       allocate(coop_struc(n)%stnlon(coop_struc(n)%nstns))

! open coord file again to read the coordinates
       open(78,file=coop_struc(n)%coorddata,form='formatted')
       write(LIS_logunit,*) 'opening coord file'
       j=0
       do i=1,nstns_total
          read(78,101) siteid, sitename, lat1, lon1
          if (idx(i) .eq. 1) then
             j=j+1
             coop_struc(n)%stnid(j) = siteid
!             coop_struc(n)%stnname(j) = trim(sitename)
             coop_struc(n)%stnlat(j) = lat1
             coop_struc(n)%stnlon(j) = lon1
             write(LIS_logunit,*) i,coop_struc(n)%stnid(j),coop_struc(n)%stnlat(j), coop_struc(n)%stnlon(j)
          endif
       enddo

       deallocate(idx)
       close(78)
       
100    format(I3,F8.1)
101    format(I6,A39,2F10.3)   
!Yuqiong: do not compute station weights; assign stn observations directly to corresponding
!         grid cell; other grid cells assigned to undefined value.
!       allocate(coop_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),coop_struc(n)%nstns))
!       call compute_stnwts(coop_struc(n)%nstns,LIS_rc%gridDesc,&
!            coop_struc(n)%stnlat,coop_struc(n)%stnlon,&
!            LIS_rc%lnc(n)*LIS_rc%lnr(n),coop_struc(n)%stnwt)
    enddo
    
  end subroutine init_COOP

end module coop_forcingMod
