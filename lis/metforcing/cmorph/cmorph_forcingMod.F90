!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module cmorph_forcingMod
!BOP
! !MODULE: cmorph_forcingMod
! 
! !REVISION HISTORY:
! 05Jan2006: Yudong Tian; start from Jon G.'s code for older versions of LIS. 
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  Climate Prediction Center (CPC)'s MORPHing technique 
!  (CMORPH). CMORPH products are produced as global fields from 
!  60 N-60S at 8km resolution at 30 minute intervals. 
!
!  The implementation in LIS has the derived data type {\tt cmorph\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[cmordir]
!    Directory containing the input data
!  \item[cmortime]
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
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_cmorph      !defines the native resolution of 
                                   !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: cmorph_struc
!EOP

  type, public :: cmorph_type_dec
     real    :: ts
     integer :: ncold, nrold   !AWIPS 212 dimensions
     character(len=LIS_CONST_PATH_LEN) :: cmorphdir   !CMOR Forcing Directory
     real*8  :: cmorphtime
     real*8  :: griduptime1
     logical :: gridchange1
     integer :: mi
     integer, allocatable :: n112(:,:)
     integer, allocatable :: n122(:,:)
     integer, allocatable :: n212(:,:)
     integer, allocatable :: n222(:,:)
     real, allocatable ::  w112(:,:),w122(:,:)
     real, allocatable ::  w212(:,:),w222(:,:)

     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 
  end type cmorph_type_dec

  type(cmorph_type_dec), allocatable :: cmorph_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_cmorph
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 06Jan2006: Yudong Tian; modification for LISv4.2
! 
! !INTERFACE:
  subroutine init_cmorph(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
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
!  Defines the native resolution of the input forcing for CMORPH
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_cmorph](\ref{readcrd_cmorph}) \newline
!     reads the runtime options specified for CMORPH data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
    allocate(cmorph_struc(LIS_rc%nnest))
    call readcrd_cmorph()

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    do n=1, LIS_rc%nnest
       cmorph_struc(n)%ts = 3600 !check 
       call LIS_update_timestep(LIS_rc, n, cmorph_struc(n)%ts)
    enddo
    
    LIS_rc%met_nf(findex) = 2

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

         cmorph_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
         cmorph_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
         cmorph_struc(n)%nIter = LIS_forecast_struc(1)%niterations

         allocate(cmorph_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
         allocate(cmorph_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       ! Regular retrospective or non-forecast mode:
       else
         cmorph_struc(n)%st_iterid = 1
         cmorph_struc(n)%en_iterId = 1
         cmorph_struc(n)%nIter = 1

         allocate(cmorph_struc(n)%metdata1(1,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
         allocate(cmorph_struc(n)%metdata2(1,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       endif

       cmorph_struc(n)%metdata1 = 0
       cmorph_struc(n)%metdata2 = 0
     
       gridDesci = 0
       gridDesci(1) = 0
       gridDesci(2) = 4948
       gridDesci(3) = 1649
       gridDesci(4) = -59.963614
       gridDesci(5) = -179.963622
       gridDesci(6) = 128
       gridDesci(7) = 59.963614
       gridDesci(8) = 179.963622
       gridDesci(9) = 0.072756669
       gridDesci(10) = 0.072771377
       gridDesci(20) = 64

       allocate(cmorph_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(cmorph_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(cmorph_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(cmorph_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       allocate(cmorph_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
       
       call conserv_interp_input(n,gridDesci,&
            cmorph_struc(n)%n112,cmorph_struc(n)%n122,cmorph_struc(n)%n212,&
            cmorph_struc(n)%n222,cmorph_struc(n)%w112,cmorph_struc(n)%w122,&
            cmorph_struc(n)%w212,cmorph_struc(n)%w222)
       
       yr1 = 2012     !grid update time
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(cmorph_struc(n)%griduptime1,&
            updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
       cmorph_struc(n)%gridchange1 = .true.
    enddo
  end subroutine init_cmorph

end module cmorph_forcingMod
