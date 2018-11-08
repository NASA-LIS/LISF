!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module capa_forcingMod
!BOP
! !MODULE: capa_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of Canadian Precipitation Analysis (CaPA)
! 
!  The implementation in LIS has the derived data type {\tt capa\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ts]
!    Frequency of the CaPA forcing data.
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[capadir]
!    Directory containing the input data
!  \item[capatime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[mi]
!    Number of points in the input grid
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES: 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_capa      !defines the native resolution of 
                           !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: capa_struc

!EOP
 
  type, public        :: capa_type_dec
     real             :: ts
     integer          :: ncold
     integer          :: nrold
     character*40     :: capadir
     real*8           :: capatime
     integer          :: mi
     integer, allocatable :: n112(:,:)
     integer, allocatable :: n122(:,:)
     integer, allocatable :: n212(:,:)
     integer, allocatable :: n222(:,:)
     real,    allocatable :: w112(:,:),w122(:,:)
     real,    allocatable :: w212(:,:),w222(:,:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type capa_type_dec

  type(capa_type_dec), allocatable :: capa_struc(:)
contains
  
!BOP
!
! !ROUTINE: defineNativeCAPA
! \label{defineNativeCAPA}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_capa(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex

! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for CAPA
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_capa](\ref{readcrd_capa}) \newline
!     reads the runtime options specified for CAPA data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    real    :: gridDesci(50)
    integer :: n

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the CaPA forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(capa_struc(LIS_rc%nnest))

    call readcrd_capa()

    do n=1, LIS_rc%nnest
       capa_struc(n)%ts = 6*60*60 
       call LIS_update_timestep(LIS_rc, n, capa_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1

    do n=1,LIS_rc%nnest
       allocate(capa_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(capa_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       capa_struc(n)%metdata1 = 0
       capa_struc(n)%metdata2 = 0

       if ( LIS_rc%met_interp(findex) == "budget-bilinear" ) then
          gridDesci = 0
!------------------------------------------------------------------------ 
! CAPA product start time to use JonG's data. 1979010100 for 00Z-06Z.
!------------------------------------------------------------------------
       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
          gridDesci = 0
          gridDesci(1) = 5
          gridDesci(2) = 1201
          gridDesci(3) = 776
          gridDesci(4) = 14.695
          gridDesci(5) = -146.683
          gridDesci(6) = 136 ! returned from getgb
          gridDesci(7) = -69.0
          gridDesci(8) = 10.000
          gridDesci(9) = 10.000
          gridDesci(10) = 0 
          gridDesci(11) = 0
          gridDesci(13) = 1 !NCEP/AFWA style of defining PS grid
          gridDesci(20) = 0.0

          capa_struc(n)%ncold = 1201
          capa_struc(n)%nrold = 776

          allocate(capa_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(capa_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(capa_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(capa_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(capa_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(capa_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(capa_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(capa_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               capa_struc(n)%n112,capa_struc(n)%n122,&
               capa_struc(n)%n212,capa_struc(n)%n222,&
               capa_struc(n)%w112,capa_struc(n)%w122,&
               capa_struc(n)%w212,capa_struc(n)%w222)
       else
          write(LIS_logunit,*) 'CAPA precip implementation only supports '
          write(LIS_logunit,*) 'conservative implementation .... '
          write(LIS_logunit,*) 'program stopping ...'
          call LIS_endrun()
       endif
    enddo
  end subroutine init_capa

end module capa_forcingMod
