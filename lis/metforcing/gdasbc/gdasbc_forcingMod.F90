!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdasbc_forcingMod
!BOP
! !MODULE: gdasbc_forcingMod
! 
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
! 
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDASBC      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdasbc_struc
!EOP

  type, public ::  gdasbc_type_dec 
     real          :: ts
     integer       :: ncold, nrold   ! AWIPS 212 dimensions
     character*50  :: gdasbc_filesrc
     character(len=LIS_CONST_PATH_LEN) :: gdasbcdir ! NLDAS-2 Forcing Directory
     real*8        :: gdasbctime1,gdasbctime2

     real                   :: gridDesc(50)
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

     integer, allocatable   :: n113(:)
     integer                :: findtime1, findtime2

     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type gdasbc_type_dec

  type(gdasbc_type_dec), allocatable :: gdasbc_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_GDASBC
! \label{init_GDASBC}
!
! !INTERFACE:
  subroutine init_GDASBC(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit,LIS_endrun
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use map_utils,      only : proj_latlon
    use LIS_forecastMod

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NLDAS-2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gdasbc](\ref{readcrd_gdasbc}) \newline
!     reads the runtime options specified for NLDAS-2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_gdasbc\_elev](\ref{read_gdasbc_elev}) \newline
!    reads the native elevation of the NLDAS-2 data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
    
    integer :: n
    
    allocate(gdasbc_struc(LIS_rc%nnest))
    call readcrd_gdasbc()

    do n=1, LIS_rc%nnest
       gdasbc_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, gdasbc_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1

  ! Set NLDAS-2 grid dimensions and extent information:
    gdasbc_struc(:)%ncold = 3600
    gdasbc_struc(:)%nrold = 1500

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

          allocate(gdasbc_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))
          allocate(gdasbc_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))

          gdasbc_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          gdasbc_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          gdasbc_struc(n)%nIter = LIS_forecast_struc(1)%niterations

       ! Regular retrospective or non-forecast mode:
       else
          allocate(gdasbc_struc(n)%metdata1(1,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))
          allocate(gdasbc_struc(n)%metdata2(1,&
                   LIS_rc%met_nf(findex),&
                   LIS_rc%ngrid(n)))

          gdasbc_struc(n)%st_iterid = 1
          gdasbc_struc(n)%en_iterId = 1
          gdasbc_struc(n)%nIter = 1

       endif

       gdasbc_struc(n)%metdata1 = 0
       gdasbc_struc(n)%metdata2 = 0

       gdasbc_struc(n)%gridDesc = 0        
       gdasbc_struc(n)%findtime1 = 0 
       gdasbc_struc(n)%findtime2 = 0 

       gdasbc_struc(n)%gridDesc(1) = 0
       gdasbc_struc(n)%gridDesc(2) = gdasbc_struc(n)%ncold
       gdasbc_struc(n)%gridDesc(3) = gdasbc_struc(n)%nrold
       gdasbc_struc(n)%gridDesc(4) = -59.95
       gdasbc_struc(n)%gridDesc(5) = -179.95
       gdasbc_struc(n)%gridDesc(6) = 128
       gdasbc_struc(n)%gridDesc(7) = 89.95
       gdasbc_struc(n)%gridDesc(8) = -179.95
       gdasbc_struc(n)%gridDesc(9) = 0.1
       gdasbc_struc(n)%gridDesc(10) = 0.1
       gdasbc_struc(n)%gridDesc(20) = 64

     ! Check for grid and interp option selected:
       if( gdasbc_struc(n)%gridDesc(9)  == LIS_rc%gridDesc(n,9) .and. &
           gdasbc_struc(n)%gridDesc(10) == LIS_rc%gridDesc(n,10).and. &
           LIS_rc%gridDesc(n,1) == proj_latlon .and. &
           LIS_rc%met_interp(findex) .ne. "neighbor" ) then
         write(LIS_logunit,*) "[ERR] The GDASBC grid was selected for the"
         write(LIS_logunit,*) "[ERR] LIS run domain; however, 'bilinear', 'budget-bilinear',"
         write(LIS_logunit,*) "[ERR] or some other unknown option was selected to spatially"
         write(LIS_logunit,*) "[ERR] downscale the grid, which will cause errors during runtime."
         write(LIS_logunit,*) "[ERR] Program stopping ..."
         call LIS_endrun()
       endif

       gdasbc_struc(n)%mi = gdasbc_struc(n)%ncold*gdasbc_struc(n)%nrold

     ! Setting up weights for spatial interpolation:
       select case( LIS_rc%met_interp(findex) )
        case( "bilinear" )
          allocate(gdasbc_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,&
               gdasbc_struc(n)%gridDesc(:),&
               gdasbc_struc(n)%n111,gdasbc_struc(n)%n121,gdasbc_struc(n)%n211,&
               gdasbc_struc(n)%n221,gdasbc_struc(n)%w111,gdasbc_struc(n)%w121,&
               gdasbc_struc(n)%w211,gdasbc_struc(n)%w221)

        case( "budget-bilinear" )

          allocate(gdasbc_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasbc_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gdasbc_struc(n)%gridDesc(:),&
               gdasbc_struc(n)%n111,gdasbc_struc(n)%n121,&
               gdasbc_struc(n)%n211,gdasbc_struc(n)%n221,&
               gdasbc_struc(n)%w111,gdasbc_struc(n)%w121,&
               gdasbc_struc(n)%w211,gdasbc_struc(n)%w221)

          allocate(gdasbc_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasbc_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasbc_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasbc_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasbc_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasbc_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasbc_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasbc_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gdasbc_struc(n)%gridDesc(:),&
               gdasbc_struc(n)%n112,gdasbc_struc(n)%n122,&
               gdasbc_struc(n)%n212,gdasbc_struc(n)%n222,&
               gdasbc_struc(n)%w112,gdasbc_struc(n)%w122,&
               gdasbc_struc(n)%w212,gdasbc_struc(n)%w222)

        case( "neighbor" )

          allocate(gdasbc_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call neighbor_interp_input(n,gdasbc_struc(n)%gridDesc(:),&
               gdasbc_struc(n)%n113)
        case default
          write(LIS_logunit,*) "[ERR] Interpolation option not specified for GDASBC"
          write(LIS_logunit,*) "[ERR] Program stopping ..."
          call LIS_endrun()
       end select     
    enddo

  end subroutine init_GDASBC
end module gdasbc_forcingMod
