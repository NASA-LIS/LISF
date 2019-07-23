!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module imerg_forcingMod
!BOP
! !MODULE: imerg_forcingMod
! 
! !REVISION HISTORY:
! 05Jan2006: Yudong Tian; start from Jon G.'s code for older versions of LIS. 
! 05Mar2015: Jonathan Case; modified for GPM IMERG product.
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  Global Precipitation Measurement's (GPM)'s Integrated Multi-satellitE Retrievals for GPM (IMERG) product.
!  IMERG products are produced as global fields at 0.1-deg resolution in separate 30-minute files in HDF5 format. 
!  As of spring 2015, data are non-missing only between 60S and 60N lat.
!
!  The implementation in LIS has the derived data type {\tt imerg\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[imergdir]
!    Directory containing the input data
!  \item[imergtime]
!    The nearest, hourly instance of the incoming 
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
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_imerg      !defines the native resolution of 
                            !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: imerg_struc
!EOP

  type, public :: imerg_type_dec
     real    :: ts
     integer :: ncold, nrold     ! IMERG dimensions
     character*100 :: imergdir   ! IMERG Forcing Directory
     character*5 :: imergver     ! IMERG version (V06B set as default)
     character*5 :: imergprd     ! IMERG product (early, late, final)
     real*8  :: imergtime
     real*8  :: griduptime1
     logical :: gridchange1
     integer :: mi
     integer, allocatable :: n112(:,:)
     integer, allocatable :: n122(:,:)
     integer, allocatable :: n212(:,:)
     integer, allocatable :: n222(:,:)
     real, allocatable    :: w112(:,:),w122(:,:)
     real, allocatable    :: w212(:,:),w222(:,:)
     real, allocatable    :: metdata1(:,:,:) 
     real, allocatable    :: metdata2(:,:,:) 
     integer              :: nIter, st_iterid,en_iterid  ! Forecast parameters
  end type imerg_type_dec

  type(imerg_type_dec), allocatable :: imerg_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_imerg
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 06Jan2006: Yudong Tian; modification for LISv4.2
! 
! !INTERFACE:
  subroutine init_imerg(findex)
! !USES: 
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only : LIS_logunit, LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_FORC_AttributesMod
    use LIS_forecastMod

    implicit none
    integer, intent(in) :: findex

    integer :: n 
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    real    :: gridDesci(50)
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for IMERG
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_imerg](\ref{readcrd_imerg}) \newline
!     reads the runtime options specified for IMERG data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
    allocate(imerg_struc(LIS_rc%nnest))

    call readcrd_imerg()

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    do n=1, LIS_rc%nnest
       imerg_struc(n)%ts = 1800
       call LIS_update_timestep(LIS_rc, n, imerg_struc(n)%ts)

       imerg_struc(n)%ncold = 3600
       imerg_struc(n)%nrold = 1800
    enddo
    
    LIS_rc%met_nf(findex) = 2

    gridDesci = 0
    gridDesci(1) = 0
    gridDesci(2) = 3600
    gridDesci(3) = 1800
    gridDesci(4) =  -89.95
    gridDesci(5) = -179.95
    gridDesci(6) = 128
    gridDesci(7) =  89.95
    gridDesci(8) = 179.95
    gridDesci(9) =  0.1
    gridDesci(10) = 0.1
    gridDesci(20) = 64

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

         imerg_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
         imerg_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
         imerg_struc(n)%nIter = LIS_forecast_struc(1)%niterations

         allocate(imerg_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(imerg_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))

       ! Regular retrospective or non-forecast mode:
       else
         imerg_struc(n)%st_iterid = 1
         imerg_struc(n)%en_iterId = 1
         imerg_struc(n)%nIter = 1

         allocate(imerg_struc(n)%metdata1(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(imerg_struc(n)%metdata2(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
       endif

       imerg_struc(n)%metdata1 = 0
       imerg_struc(n)%metdata2 = 0

       allocate(imerg_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(imerg_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(imerg_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(imerg_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(imerg_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(imerg_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(imerg_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(imerg_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       
       call conserv_interp_input( n, gridDesci, &
            imerg_struc(n)%n112,imerg_struc(n)%n122,&
            imerg_struc(n)%n212,imerg_struc(n)%n222,&
            imerg_struc(n)%w112,imerg_struc(n)%w122,&
            imerg_struc(n)%w212,imerg_struc(n)%w222)
       
       yr1 = 2012     !grid update time
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(imerg_struc(n)%griduptime1,&
            updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
       imerg_struc(n)%gridchange1 = .true.
    enddo
  end subroutine init_imerg

end module imerg_forcingMod
