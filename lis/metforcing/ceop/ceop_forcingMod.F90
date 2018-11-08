!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ceop_forcingMod
!BOP
! !MODULE: ceop_forcingMod
! 
! !DESCRIPTION: 
!  Contains routines and data structures that are used for the 
!  implementation of the station data from various CEOP stations. 
!  The stations report estimates of meteorological forcing terms, 
!  which is spatially interpolated using the inverse distance 
!  weighting scheme (IDW). 
! 
!  The implementation in LIS has the derived data type {\tt ceop\_struc}
!  that includes the variables to specify the runtime options, and the
!  calculation of weights for spatial interpolation.
!
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_CEOP      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ceop_struc
!EOP

  type, public :: ceop_type_dec
     real          :: ts
     character*40  :: ceopdir 
     character*40  :: metadata 
     real          :: undef
     real*8        :: starttime,ceoptime1,ceoptime2
     integer       :: findtime1,findtime2,nstns
     logical       :: startRead
     integer       :: location
     real, allocatable :: stnwt(:,:)
     real, allocatable :: stnlat(:),stnlon(:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type ceop_type_dec

  type(ceop_type_dec), allocatable :: ceop_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_CEOP
! \label{init_CEOP}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_CEOP(findex)
! !USES:
    use LIS_coreMod,only : LIS_rc
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
    integer, intent(in) :: findex

! !DESCRIPTION: 
!  This routines reads the runtime configurations for using the
!  CEOP station data. Using the metadata provided for the 
!  stations, this routine invokes the call to compute the
!  interpolation weights to be later used. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_ceop](\ref{readcrd_ceop}) \newline
!     reads the runtime options specified for CEOP station data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[compute\_stnwts](\ref{compute_stnwts}) \newline
!    computes the weights for spatial interpolation
!  \end{description}
!EOP

    integer :: i, styr,stmo,std,sth,stm,sts,tinc
    character*20, allocatable :: stnid(:)
    real    :: gmt
    integer :: doy
    integer :: n 

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the CEOP forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(ceop_struc(LIS_rc%nnest))
    call readcrd_ceop()
    
    do n=1, LIS_rc%nnest
       ceop_struc(n)%ts = 60*60 
       call LIS_update_timestep(LIS_rc, n, ceop_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9 

    do n=1,LIS_rc%nnest

       allocate(ceop_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(ceop_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       ceop_struc(n)%metdata1 = 0
       ceop_struc(n)%metdata2 = 0

       write(LIS_logunit,*) '--------------------------------------------------------'
       write(LIS_logunit,*) 'Opening CEOP Meta Data'
       open(78,file=ceop_struc(n)%metadata,status='old')
       read(78,100) ceop_struc(n)%nstns, ceop_struc(n)%undef, styr,stmo,std,sth,stm,sts,&
            tinc
       allocate(stnid(ceop_struc(n)%nstns))
       allocate(ceop_struc(n)%stnlat(ceop_struc(n)%nstns))
       allocate(ceop_struc(n)%stnlon(ceop_struc(n)%nstns))
       write(LIS_logunit,*) 'Number of stations: ',ceop_struc(n)%nstns
       write(LIS_logunit,*) 'Undef value... ',ceop_struc(n)%undef
       write(LIS_logunit,*) 'Starting time..',styr,'/',stmo,'/',std,' ',&
            sth,':',stm,':',sts
       
       call LIS_date2time(ceop_struc(n)%starttime,doy,gmt,styr,stmo,std,&
            sth,stm,sts)
       write(LIS_logunit,*) 'Observation time interval (min) ',tinc
       do i=1,ceop_struc(n)%nstns
          read(78,101) stnid(i),ceop_struc(n)%stnlat(i),ceop_struc(n)%stnlon(i)
          write(LIS_logunit,*) i,stnid(i),ceop_struc(n)%stnlat(i),ceop_struc(n)%stnlon(i)
          if((ceop_struc(n)%stnlat(i).gt.LIS_rc%gridDesc(n,7)) &
               .or.(ceop_struc(n)%stnlon(i).gt.LIS_rc%gridDesc(n,8)) &
               .or.(ceop_struc(n)%stnlat(i).lt.LIS_rc%gridDesc(n,4)) &
               .or.(ceop_struc(n)%stnlon(i).lt.LIS_rc%gridDesc(n,5))) then
             write(LIS_logunit,*) 'Station ',stnid(i),'(',ceop_struc(n)%stnlat(i),',',&
                  ceop_struc(n)%stnlon(i),')',&
                  'is not within bounds..'
             write(LIS_logunit,*) 'stopping program...'
             call LIS_endrun
          endif
       enddo
       write(LIS_logunit,*) '--------------------------------------------------------'
       close(78)
       
100    format(I3,F9.3,I5,6(I3))
101    format(A20,2(F8.2))
       
       allocate(ceop_struc(n)%stnwt(LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%nstns))
       call compute_stnwts(ceop_struc(n)%nstns,LIS_rc%gridDesc,&
            ceop_struc(n)%stnlat,ceop_struc(n)%stnlon,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n),ceop_struc(n)%stnwt)
       deallocate(stnid)
    enddo
    
  end subroutine init_CEOP

end module ceop_forcingMod
