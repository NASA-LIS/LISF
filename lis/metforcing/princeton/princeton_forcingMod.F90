!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module princeton_forcingMod
!BOP
! !MODULE: princeton_forcingMod
! 
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the 50-yr Princeton
!  data (Sheffield et al. 2006). The data is global 1 degree dataset in latlon
!  projection, and at 3 hourly intervals. The derived
!  data type {\tt princeton\_struc}
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
!  \item[princetondir]
!    Directory containing the input data
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
!  Sheffield, J., G. Goteti, and E. F. Wood, 2006: Development of a 50-yr 
!  high-resolution global dataset of meteorological forcings for land surface 
!  modeling, J. Climate, 19 (13), 3088-3111 \newline
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_PRINCETON      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: princeton_struc

!EOP

  type, public :: princeton_type_dec
     real                   :: ts
     integer                :: ncold, nrold   
     character(len=LIS_CONST_PATH_LEN) :: princetondir
     character*100          :: elevfile
     character*100          :: version
     integer                :: mi
     real*8                 :: princetontime1,princetontime2
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
     integer                :: findtime1, findtime2

     integer           :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type princeton_type_dec

  type(princeton_type_dec), allocatable :: princeton_struc(:) 

contains
  
!BOP
!
! !ROUTINE: init_PRINCETON
! \label{init_PRINCETON}
!
! !REVISION HISTORY: 
! 26Jan2007: Hiroko Kato; Initial Specification
! 15 May 2017: Bailing Li; Added changes for reading in version 2.2 data
! 22 Oct 2018: Daniel Sarmiento; Added changes to support version 3 data
! 
! !INTERFACE:
  subroutine init_PRINCETON(findex)
! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_update_timestep
   use LIS_logMod,     only : LIS_logunit, LIS_endrun
   use LIS_forecastMod

   implicit none
   integer, intent(in)  :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for PRINCETON 
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_princeton](\ref{readcrd_princeton}) \newline
!     reads the runtime options specified for PRINCETON data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_princeton\_elev](\ref{read_princeton_elev}) \newline
!    reads the native elevation of the princeton data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
    real     :: gridDesci(50)
    integer  :: n 

    gridDesci = 0.0
    allocate(princeton_struc(LIS_rc%nnest))
    call readcrd_princeton()

    do n=1, LIS_rc%nnest
       princeton_struc(n)%ts = 3*3600
       call LIS_update_timestep(LIS_rc, n, princeton_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9 

    do n=1,LIS_rc%nnest
       !Set dataset dimensions
       if (princeton_struc(n)%version == "2" .OR. princeton_struc(n)%version == "2.2") then
          princeton_struc(:)%ncold = 360
          princeton_struc(:)%nrold = 180
       elseif (princeton_struc(n)%version == "3") then 
          !Dimensions of driver data changed from versions 2.x to version 3
          princeton_struc(:)%ncold = 1440
          princeton_struc(:)%nrold = 600
       endif

       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then

          if(mod(LIS_rc%nensem(n),&
               LIS_forecast_struc(1)%niterations).ne.0) then
             write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
             write(LIS_logunit,*) '[ERR] of the number of iterations '
             write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
             write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
             call LIS_endrun()
          endif

          allocate(princeton_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(princeton_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

          princeton_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          princeton_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          princeton_struc(n)%nIter = LIS_forecast_struc(1)%niterations

       ! Regular retrospective or non-forecast mode:
       else
          allocate(princeton_struc(n)%metdata1(1,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(princeton_struc(n)%metdata2(1,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

          princeton_struc(n)%st_iterid = 1
          princeton_struc(n)%en_iterId = 1
          princeton_struc(n)%nIter = 1

       endif

       princeton_struc(n)%metdata1 = 0
       princeton_struc(n)%metdata2 = 0

       gridDesci(1) = 0
       gridDesci(2) = princeton_struc(n)%ncold
       gridDesci(3) = princeton_struc(n)%nrold
       !Define driver data domains
       if (princeton_struc(n)%version == "2" .OR. princeton_struc(n)%version == "2.2") then 
          gridDesci(4) = -89.50
          gridDesci(5) = -179.50
          gridDesci(6) = 128
          gridDesci(7) = 89.50
          gridDesci(8) = 179.50
          gridDesci(9) = 1.00
          gridDesci(10) = 1.00
       elseif (princeton_struc(n)%version == "3") then
          gridDesci(4) = -59.875
          gridDesci(5) = -179.875
          gridDesci(6) = 128
          gridDesci(7) = 89.875
          gridDesci(8) = 179.875
          gridDesci(9) = 0.25
          gridDesci(10) = 0.25
       endif
       gridDesci(20) = 0
       princeton_struc(n)%mi = princeton_struc(n)%ncold*princeton_struc(n)%nrold

     ! Setting up weights for Interpolation
       select case( LIS_rc%met_interp(findex) )

         case( "bilinear" )
          allocate(princeton_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               princeton_struc(n)%n111,princeton_struc(n)%n121,&
               princeton_struc(n)%n211,princeton_struc(n)%n221,&
               princeton_struc(n)%w111,princeton_struc(n)%w121,&
               princeton_struc(n)%w211,princeton_struc(n)%w221)

         case( "budget-bilinear" )
          allocate(princeton_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(princeton_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               princeton_struc(n)%n111,princeton_struc(n)%n121,&
               princeton_struc(n)%n211,princeton_struc(n)%n221,&
               princeton_struc(n)%w111,princeton_struc(n)%w121,&
               princeton_struc(n)%w211,princeton_struc(n)%w221)

          allocate(princeton_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(princeton_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(princeton_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(princeton_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(princeton_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(princeton_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(princeton_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(princeton_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               princeton_struc(n)%n112,princeton_struc(n)%n122,&
               princeton_struc(n)%n212,princeton_struc(n)%n222,&
               princeton_struc(n)%w112,princeton_struc(n)%w122,&
               princeton_struc(n)%w212,princeton_struc(n)%w222)

       case default
         write(LIS_logunit,*) "[ERR] User-input issue with PRINCETON forcing ..."
         write(LIS_logunit,*) " -- Currently only supported interpolation options include:"
         write(LIS_logunit,*) "  - bilinear or budget-bilinear - "
         write(LIS_logunit,*) " Program stopping ... "
         call LIS_endrun
       end select

       if( LIS_rc%met_ecor(findex).ne."none" ) then 
          call read_princeton_elev(n,findex)
       endif

    enddo

  end subroutine init_Princeton
end module princeton_forcingMod
