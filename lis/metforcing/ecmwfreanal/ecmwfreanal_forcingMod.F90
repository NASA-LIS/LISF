!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ecmwfreanal_forcingMod
!BOP
! !MODULE: ecmwfreanal_forcingMod
! 
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the bias-corrected ECMWF atmospheric reanalysis
!  data (Berg et al. 2003). The data is global 0.5 degree dataset in latlon
!  projection, and at 6 hourly intervals. The derived
!  data type {\tt ecmwfreanal\_struc}
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[fmodeltime1]
!    The nearest, previous 6 hour instance of the incoming 
!    data (as a real time). 
!  \item[fmodeltime2]
!    The nearest, next 6 hour instance of the incoming 
!    data (as a real time).
!  \item[remask1d]
!    The data used to remask the input data to the LIS mask. 
!  \item[ecmwfreanaldir]
!    Directory containing the input data
!  \item[emaskfile]
!    File containing the 0.5 deg land-sea mask used in the input data. 
!  \item[elevfile]
!    File with the elevation definition for the input data. 
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n12,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
! 
!  Berg, A.A., J.S. Famiglietti, J.P. Walker, and P.R. Houser, 2003: 
!  Impact of bias correction to reanalysis products on simulations of
!  North American soil moisture and hydrological fluxes. Journal of 
!  Geophysical Research, 108, 4490, DOI: 10.1029/2002JD003334. \newline
!
! !USES: 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_ECMWFREANAL      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ecmwfreanal_struc

!EOP

  type, public :: ecmwfreanal_type_dec
     real                   :: ts
     integer                :: ncold, nrold   !AWIPS 212 dimensions
     real*8                 :: fmodeltime1,fmodeltime2
     integer                :: remask1d(720*360)
     character*100          :: ecmwfreanaldir
     character*100           :: emaskfile !1/2deg ECMWFREANAL Land-Sea Mask File
     character*100          :: elevfile
     integer                :: mi
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
     integer            :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type ecmwfreanal_type_dec

  type(ecmwfreanal_type_dec), allocatable :: ecmwfreanal_struc(:) 
contains
  
!BOP
!
! !ROUTINE: init_ECMWFREANAL
! \label{init_ECMWFREANAL}
!
! !REVISION HISTORY: 
! 26Jan2004: Sujay Kumar; Init_ial Specification
! 
! !INTERFACE:
  subroutine init_ECMWFREANAL(findex)
! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_update_timestep
   use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !AGRUMENTS: 
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for ECMWFREANAL 
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_ecmwfreanal](\ref{readcrd_ecmwfreanal}) \newline
!     reads the runtime options specified for ECMWFREANAL data
!   \item[readmask\_ecmwfreanal](\ref{readmask_ecmwfreanal}) \newline
!     reads the 0.5 degree land sea mask
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_ecmwfreanal\_elev](\ref{read_ecmwfreanal_elev}) \newline
!    reads the native elevation of the ecmwfreanal data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
    real     :: gridDesci(50)
    integer  :: n 

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the ECMWF-reanalysis forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    gridDesci = 0.0
    allocate(ecmwfreanal_struc(LIS_rc%nnest))

    do n=1, LIS_rc%nnest
       ecmwfreanal_struc(n)%ts = 6*3600 
       call LIS_update_timestep(LIS_rc, n, ecmwfreanal_struc(n)%ts)
    enddo

    call readcrd_ecmwfreanal()

    LIS_rc%met_nf(findex) = 9

    do n=1,LIS_rc%nnest
       
       allocate(ecmwfreanal_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(ecmwfreanal_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       ecmwfreanal_struc(n)%metdata1 = 0
       ecmwfreanal_struc(n)%metdata2 = 0

       ecmwfreanal_struc(n)%findtime1 = 0 
       ecmwfreanal_struc(n)%findtime2 = 0

       gridDesci(1) = 0
       gridDesci(2) = ecmwfreanal_struc(n)%ncold
       gridDesci(3) = ecmwfreanal_struc(n)%nrold
       gridDesci(4) = -89.750
       gridDesci(5) = -179.750
       gridDesci(6) = 128
       gridDesci(7) = 89.750
       gridDesci(8) = 179.750
       gridDesci(9) = 0.500
       gridDesci(10) = 0.500
       gridDesci(20) = 64
       call readmask_ecmwfreanal(n)
       ecmwfreanal_struc(n)%mi = ecmwfreanal_struc(n)%ncold*&
            ecmwfreanal_struc(n)%nrold
!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(ecmwfreanal_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               ecmwfreanal_struc(n)%n111,ecmwfreanal_struc(n)%n121,&
               ecmwfreanal_struc(n)%n211,ecmwfreanal_struc(n)%n221,&
               ecmwfreanal_struc(n)%w111,ecmwfreanal_struc(n)%w121,&
               ecmwfreanal_struc(n)%w211,ecmwfreanal_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(ecmwfreanal_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               ecmwfreanal_struc(n)%n111,ecmwfreanal_struc(n)%n121,&
               ecmwfreanal_struc(n)%n211,ecmwfreanal_struc(n)%n221,&
               ecmwfreanal_struc(n)%w111,ecmwfreanal_struc(n)%w121,&
               ecmwfreanal_struc(n)%w211,ecmwfreanal_struc(n)%w221)

          allocate(ecmwfreanal_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               ecmwfreanal_struc(n)%n112,ecmwfreanal_struc(n)%n122,&
               ecmwfreanal_struc(n)%n212,ecmwfreanal_struc(n)%n222,&
               ecmwfreanal_struc(n)%w112,ecmwfreanal_struc(n)%w122,&
               ecmwfreanal_struc(n)%w212,ecmwfreanal_struc(n)%w222)
       endif
       if(trim(LIS_rc%met_ecor(findex)).ne."none") then 
          call read_ecmwfreanal_elev(n,findex)
       endif
    enddo

  end subroutine init_Ecmwfreanal
end module ecmwfreanal_forcingMod
