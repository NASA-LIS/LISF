!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module agrrad_forcingMod
!BOP
! !MODULE: agrrad_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the AGRMET radiation data. The AGRMET
!  data is in polar stereographic projection at 48km resolution.
!
!  The implementation in LIS has the derived data type {\tt agrrad\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ts]
!    Frequency of the AGRRAD forcing data.
!  \item[agrdir]
!    Directory containing the input data.
!  \item[agrtime1]
!    The nearest, previous hourly instance of the incoming 
!    data (as a real time).
!  \item[agrtime2]
!    The nearest, next hourly instance of the incoming 
!    data (as a real time).
!  \item[changetime1]
!    The time to switch the input resolution to T170
!  \item[changetime2]
!    The time to switch the input resolution to T170
!  \item[findtime1, findtime2]
!    Flags to indicate which time is to be read for 
!    temporal interpolation.
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[fillflag1]
!    Flag used to control initializing smask1.
!  \item[smask1]
!    Array used to mark mismatches between the AGRRAD forcing and
!    the land/sea mask.
!  \end{description}
!
! !USES: 
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_agrrad      !defines the native resolution of 
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: agrrad_struc

!EOP

  type, public        :: agrrad_type_dec
     real             :: ts
     real*8           :: agrtime1,agrtime2
     real*8           :: changetime1,changetime2
     integer          :: findtime1,findtime2
     character*40     :: agrdir
     integer, allocatable :: n111(:)
     integer, allocatable :: n121(:)
     integer, allocatable :: n211(:)
     integer, allocatable :: n221(:)
     real, allocatable    :: w111(:),w121(:)
     real, allocatable    :: w211(:),w221(:)
     integer, allocatable :: smask1(:,:)
     logical           :: fillflag1

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type agrrad_type_dec

  type(agrrad_type_dec), allocatable :: agrrad_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_agrrad
! \label{init_agrrad}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_agrrad(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for AGRMET 
!  radiation data. The grid description arrays are based 
!  on the decoding schemes used by NCEP and followed in 
!  the LIS interpolation schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_agrrad](\ref{readcrd_agrrad}) \newline
!     reads the runtime options specified for AGRMET radiation data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    integer :: n
    real    :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the AGRMET radiation forcing'
       write(LIS_logunit,*) '[ERR]  reader is not set up to run in forecast'
       write(LIS_logunit,*) '[ERR]  mode.  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(agrrad_struc(LIS_rc%nnest))

    call readcrd_agrrad()

    do n=1, LIS_rc%nnest
       agrrad_struc(n)%ts = 3*60*60 
       call LIS_update_timestep(LIS_rc, n, agrrad_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 2

    do n=1,LIS_rc%nnest

       allocate(agrrad_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(agrrad_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       agrrad_struc(n)%metdata1 = 0
       agrrad_struc(n)%metdata2 = 0

       allocate(agrrad_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       agrrad_struc(n)%fillflag1 = .true.

       agrrad_struc(n)%agrtime1 = 3000.0
       agrrad_struc(n)%agrtime2 = 0.0
       
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 720
       gridDesci(3) = 361
       gridDesci(4) = -89.75
       gridDesci(5) = -179.75
       gridDesci(7) = 89.75
       gridDesci(8) = 179.75
       gridDesci(6) = 128
       gridDesci(9) = 0.50
       gridDesci(10) = 0.50
       gridDesci(20) = 64
       
       if(LIS_rc%met_interp(findex).eq."bilinear") then 
          allocate(agrrad_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(agrrad_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(agrrad_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(agrrad_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(agrrad_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(agrrad_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(agrrad_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(agrrad_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               agrrad_struc(n)%n111,agrrad_struc(n)%n121,&
               agrrad_struc(n)%n211,agrrad_struc(n)%n221,&
               agrrad_struc(n)%w111,agrrad_struc(n)%w121,&
               agrrad_struc(n)%w211,agrrad_struc(n)%w221)
       endif

       yr1 = 2002     !grid update time
       mo1 = 5
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(agrrad_struc(n)%changetime1,updoy,upgmt,&
            yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2002     !grid update time
       mo1 = 11
       da1 = 6
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(agrrad_struc(n)%changetime2,updoy,upgmt,&
            yr1,mo1,da1,hr1,mn1,ss1 )

    enddo
  end subroutine init_agrrad

end module agrrad_forcingMod
