!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module petusgs_forcingMod

!BOP
! !MODULE: petusgs_forcingMod
! 
! !REVISION HISTORY:
! 05Mar2012: K. Arsenault;  Added FEWSNET PET dataset as a forcing option
! 25Oct2013: K. Arsenault;  Added PET USGS to LIS7
! 02May2014; K. Arsenault;  Added climatological PET file option
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the potential evapotranspiration (PET) data 
!  from the US Geological Survey (USGS) method for estimating PET from
!  NCEP GDAS 1.0 files and fields.  The files are global and summed daily
!  (6Z GDAS forecast to next day forecast), reflecting a 00Z timestamp. 
!
!  The implementation in LIS has the derived data type {\tt petusgs\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[petdir]
!    Directory containing the input data
!  \item[pettype]
!    Option to select PET file type: climatology or current (retrospective)
!  \item[pettime]
!    The nearest (daily) instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution to T170
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
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_petusgs     !defines the native resolution of 
                                    !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: petusgs_struc
!EOP

  type, public :: petusgs_type_dec

     character(len=LIS_CONST_PATH_LEN) :: petdir    ! USGS PET Forcing Directory
     character*20      :: pettype     ! USGS PET File type (climatology|current)
     real*8            :: pettime
     real*8            :: griduptime1
     logical           :: gridchange1
     integer           :: mi
     integer, allocatable :: n111(:,:)
     integer, allocatable :: n121(:,:)
     integer, allocatable :: n211(:,:)
     integer, allocatable :: n221(:,:)
     real, allocatable    :: w111(:,:),w121(:,:)
     real, allocatable    :: w211(:,:),w221(:,:)

     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type petusgs_type_dec

  type(petusgs_type_dec), allocatable :: petusgs_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_petusgs
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar;  Initial Specification
! 11Mar2012: K. Arsenault; Expanded for use with USGSP PET data
! 
! !INTERFACE:
  subroutine init_petusgs(findex)

! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
    use LIS_forecastMod

    implicit none
    integer, intent(in) :: findex

    integer :: n 
    real    :: gridDesci(50)
    integer :: updoy,r1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt

!-----------------------------
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for USGS PET 
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readpetusgscrd](\ref{readpetusgscrd}) \newline
!     reads the runtime options specified for USGS PET data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    allocate(petusgs_struc(LIS_rc%nnest))
    call readpetusgscrd()

    LIS_rc%met_nf(findex) = 1  ! number of met variables in USGS PET forcing

    do n=1,LIS_rc%nnest   

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

         petusgs_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
         petusgs_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
         petusgs_struc(n)%nIter = LIS_forecast_struc(1)%niterations

         allocate(petusgs_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(petusgs_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))

       ! Regular retrospective or non-forecast mode:
       else 
         petusgs_struc(n)%st_iterid = 1
         petusgs_struc(n)%en_iterId = 1
         petusgs_struc(n)%nIter = 1

         allocate(petusgs_struc(n)%metdata1(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(petusgs_struc(n)%metdata2(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
       endif
    
       petusgs_struc(n)%metdata1 = 0
       petusgs_struc(n)%metdata2 = 0

       gridDesci     = 0
       gridDesci(1)  = 0
       gridDesci(2)  = 360
       gridDesci(3)  = 181
       gridDesci(4)  = -90.00
       gridDesci(5)  = -180.00
       gridDesci(6)  = 128
       gridDesci(7)  = 90.00
       gridDesci(8)  = 179.00
       gridDesci(9)  = 1.00
       gridDesci(10) = 1.00
       gridDesci(20) = 64

    !- Timestep call:
!       call LIS_update_timestep( LIS_rc, n, petusgs_struc(n)%ts )
       call LIS_update_timestep( LIS_rc, n, 86400.0 )

    !- Setting up weights for spatial interpolation:

       allocate(petusgs_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petusgs_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petusgs_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petusgs_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petusgs_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petusgs_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petusgs_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(petusgs_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       
       if( LIS_rc%met_interp(findex) == "bilinear" ) then

          call bilinear_interp_input( n,gridDesci,&
             petusgs_struc(n)%n111,petusgs_struc(n)%n121, &
             petusgs_struc(n)%n211,petusgs_struc(n)%n221, &
             petusgs_struc(n)%w111,petusgs_struc(n)%w121, &
             petusgs_struc(n)%w211,petusgs_struc(n)%w221 )

       elseif( LIS_rc%met_interp(findex) == "budget-bilinear" ) then

          call conserv_interp_input( n,gridDesci,&
             petusgs_struc(n)%n111,petusgs_struc(n)%n121,petusgs_struc(n)%n211,&
             petusgs_struc(n)%n221,petusgs_struc(n)%w111,petusgs_struc(n)%w121,&
             petusgs_struc(n)%w211,petusgs_struc(n)%w221 )


       elseif( LIS_rc%met_interp(findex) == "neighbor" ) then

          call neighbor_interp_input(n, gridDesci,&
             petusgs_struc(n)%n111 )

       else
          write(LIS_logunit,*) 'This interp. option not supported for PET USGS'
          write(LIS_logunit,*) 'Program stopping ... '
          call LIS_endrun()

       endif 

       
!  ... If USGS updates the resolution to the GDAS-based PET files, the grid changes will update here:
!       yr1 = 2012     !grid update time
!       mo1 = 10
!       da1 = 29
!       hr1 = 12
!       mn1 = 0; ss1 = 0
!       call LIS_date2time(petusgs_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
!       petusgs_struc(n)%gridchange1 = .true.

    enddo
  end subroutine init_petusgs

end module petusgs_forcingMod
