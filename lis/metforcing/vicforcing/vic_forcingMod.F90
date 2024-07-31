!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "vic411_atmos_forcing.h"
module vic_forcingMod
!BOP
! !MODULE: vic_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of VIC formatted forcing data.
!
!  The implementation in LIS has the derived data type {\tt vicforcing\_struc}
!  that includes the variables that specify the runtime options.
!
!  They are described below: 
!  \begin{description}
!  \item[vicdir]
!    Directory containing the input data
!  \item[forcingInterval]
!    Interval of the forcing data [hours].  This value should be either
!    VIC's model time-step for the energy-balance mode or the snow time-step
!    for the water-balance mode.
!  \item[NC]
!    specifies the number of columns (longitude grid-cells) for the gridded
!    VIC data (not used when read\_gridded == 0)
!  \item[NR]
!    specifies the number of rows (latitude grid-cells) for the gridded
!    VIC data (not used when read\_gridded == 0)
!  \item[slat]
!    specifies the lower left latitude for the gridded VIC forcing domain
!  \item[slon]
!    specifies the lower left longitude for the gridded VIC forcing domain
!  \item[elat]
!    specifies the upper right latitude for the gridded VIC forcing domain
!  \item[elon]
!    specifies the upper right longitude for the gridded VIC forcing domain
!  \item[dlon]
!    specifies the domain resolution (dx) for the gridded VIC forcing domain
!  \item[dlat]
!    specifies the domain resolution (dy) for the gridded VIC forcing domain
!  \item[time1]
!    The nearest, previous forcing interval instance of the incoming 
!    data (as a real time). 
!  \item[time2]
!    The nearest, next forcing interval instance of the incoming 
!    data (as a real time).
!  \end{description}
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: defineNativevicforcing      !defines the native resolution of 
                                        !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: vicforcing_struc
!EOP

  type, public :: vicforcing_type_dec
     character(len=LIS_CONST_PATH_LEN) :: vicdir
     integer       :: forcingInterval
     integer       :: NC
     integer       :: NR
     real          :: slat
     real          :: slon
     real          :: elat
     real          :: elon
     real          :: dlon
     real          :: dlat
     real*8        :: time1
     real*8        :: time2

     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)
     
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type vicforcing_type_dec
  
  type(vicforcing_type_dec), allocatable :: vicforcing_struc(:)
contains
  
!BOP
!
! !ROUTINE: defineNativevicforcing
!  \label{defineNativevicforcing}
! !INTERFACE:
   subroutine defineNativevicforcing(findex)
! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain

! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for the VIC formatted
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}. Based on the specified data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights.
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readviccrd](\ref{readviccrd}) \newline
!     reads the runtime options specified for VIC formatted data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
   implicit none

   integer,  intent(in)     :: findex

   real :: gridDesci(LIS_rc%nnest, 50)
   integer :: n 

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the VIC-specific forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

   allocate(vicforcing_struc(LIS_rc%nnest))

   call readviccrd()

   LIS_rc%met_nf(findex) = NUM_ATMOS_FORCING

   gridDesci = 0 

   do n=1,LIS_rc%nnest

      allocate(vicforcing_struc(n)%metdata1(LIS_rc%met_nf(findex),&
           LIS_rc%ngrid(n)))
      allocate(vicforcing_struc(n)%metdata2(LIS_rc%met_nf(findex),&
           LIS_rc%ngrid(n)))

      vicforcing_struc(n)%metdata1 = 0
      vicforcing_struc(n)%metdata2 = 0

      gridDesci(n,1) = 0
      gridDesci(n,2) = vicforcing_struc(n)%NC
      gridDesci(n,3) = vicforcing_struc(n)%NR
      gridDesci(n,4) = vicforcing_struc(n)%slat
      gridDesci(n,5) = vicforcing_struc(n)%slon
      gridDesci(n,6) = 128
      gridDesci(n,7) = vicforcing_struc(n)%elat
      gridDesci(n,8) = vicforcing_struc(n)%elon
      gridDesci(n,9) = vicforcing_struc(n)%dlon
      gridDesci(n,10) = vicforcing_struc(n)%dlat
      gridDesci(n,20) = 0

!Setting up weights for Interpolation
      if(LIS_rc%met_interp(findex).eq."bilinear") then 

         allocate(vicforcing_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci(n,:),&
              vicforcing_struc(n)%n111,vicforcing_struc(n)%n121,&
              vicforcing_struc(n)%n211,vicforcing_struc(n)%n221,&
              vicforcing_struc(n)%w111,vicforcing_struc(n)%w121,&
              vicforcing_struc(n)%w211,vicforcing_struc(n)%w221)
      elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 

         allocate(vicforcing_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(vicforcing_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci(n,:),&
              vicforcing_struc(n)%n111,vicforcing_struc(n)%n121,&
              vicforcing_struc(n)%n211,vicforcing_struc(n)%n221,&
              vicforcing_struc(n)%w111,vicforcing_struc(n)%w121,&
              vicforcing_struc(n)%w211,vicforcing_struc(n)%w221)

         allocate(vicforcing_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(vicforcing_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(vicforcing_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(vicforcing_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(vicforcing_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(vicforcing_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(vicforcing_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(vicforcing_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

         call conserv_interp_input(n,gridDesci(n,:),&
              vicforcing_struc(n)%n112,vicforcing_struc(n)%n122,&
              vicforcing_struc(n)%n212,vicforcing_struc(n)%n222,&
              vicforcing_struc(n)%w112,vicforcing_struc(n)%w122,&
              vicforcing_struc(n)%w212,vicforcing_struc(n)%w222)
      endif
   enddo

   end subroutine defineNativevicforcing
end module vic_forcingMod

