!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module rdhm356_forcingMod

!BOP
! !MODULE: rdhm356_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation and temperature data from the
!  National Center for Environmental Prediction (NCEP) dmip II  
!  (dmip2) Doppler Radar+gage product.  The dmip II is a national
!  level product, on an hourly interval, and supplements mainly
!  the which base forcing precipitation (e.g., NLDAS) is being used.
! 
!  The implementation in LDT has the derived data type {\tt rdhm356\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncol]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrow]
!    Number of rows (along the north south dimension) for the input data
!  \item[rdhm356dir]
!    Directory containing the input data
!  \item[rdhm356time]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \end{description}
!
! !USES: 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_rdhm356      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: rdhm356_struc_precip, rdhm356_struc_temper, const_wind 

!EOP

  type, public ::  rdhm356_type_dec
     real               :: ts
     integer            :: ncol                 ! Number of cols
     integer            :: nrow                 ! Number of rows
     character*40       :: rdhm356dir           ! dmip II Directory
     real*8             :: scale 
     real*8             :: rdhm356time1         ! Time of first incoming file
     real*8             :: rdhm356time2         ! Time of second incoming file
     real*8             :: griduptime1          ! Designated time of dmipII grid change
     logical            :: gridchange1          ! Flag for when grid change occurs
     integer            :: mi                   ! Number of points in the input grid

     ! added by Shugong Wang 
     real               :: lower_left_hrapx
     real               :: lower_left_hrapy
     real               :: upper_right_hrapx
     real               :: upper_right_hrapy
     real               :: hrap_resolution 
     real               :: undef_value 
     
     ! == Arrays for Bilinear Interpolation option (=1)
     integer, allocatable   :: n111(:), n121(:)
     integer, allocatable   :: n211(:), n221(:)
     real, allocatable      :: w111(:), w121(:) 
     real, allocatable      :: w211(:), w221(:) 
     ! == Arrays for Budget Bilinear Interpolation option (=2)
     integer, allocatable   :: n112(:,:), n122(:,:) ! Spatial Interpolation weights
     integer, allocatable   :: n212(:,:), n222(:,:) ! " "
     real,    allocatable   :: w112(:,:), w122(:,:) ! " "
     real,    allocatable   :: w212(:,:), w222(:,:) ! " "

     integer, allocatable   :: n113(:)


  end type rdhm356_type_dec

  type(rdhm356_type_dec), allocatable :: rdhm356_struc_precip(:)       
! Precip Main Pointer array (number of domains)
  type(rdhm356_type_dec), allocatable :: rdhm356_struc_temper(:)       
! Temperature Main Pointer array (number of domains)
  real, allocatable :: const_wind(:)                                            
! constant wind speed (m/s)

contains
  
!BOP
!
! !ROUTINE: init_rdhm356
! \label{init_rdhm356}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 25May2006: Kristi Arsenault; Data and code implementation
! 20Jul2006: Brian Cosgrove; Conversion to DMIP2
! 03May2010: Soni Yatheendradas; Precip and Temper. input grids now can have 
!            different extents and different from the run-domain extent, as 
!            per the new input grids posted onto the DMIP2 website for Sierra 
!            Nevada (hardcoding input grid extents for now, will make flexible 
!            later by making a call to read_xmrg2 and read_xmrgtemp2 from here) 
! 
! 17Dec2013: Shugong Wang; (1) replace old stlye XMRG reader with LDT_XMRG_Reader
!                          (2) 
! !INTERFACE:
  subroutine init_rdhm356(findex)

! !USES: 
   use LDT_coreMod,    only : LDT_rc, LDT_domain
   use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
   use LDT_xmrg_reader

   implicit none
   integer,  intent(in) :: findex
!
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for dmip2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_rdhm356](\ref{readcrd_rdhm356}) \newline
!     reads the runtime options specified for dmip2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!     computes the neighbor, weights for bilinear interpolation
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format - time of grid change 
!  \end{description}
!
!EOP
    real :: gridDesciPrecip(LDT_rc%nnest, 20)
    real :: gridDesciTemper(LDT_rc%nnest, 20)
    integer :: n
    real :: hrap_x_precip, hrap_y_precip, rlon_precip, rlat_precip
    real :: hrap_x_temper, hrap_y_temper, rlon_temper, rlat_temper
    real :: upgmt

    real :: rlat_precip_ll, rlat_temper_ll
    real :: rlon_precip_ll, rlon_temper_ll
    integer :: t 

    ! - Forcing Data structure -- Allocate for different LDT_domains
    allocate ( rdhm356_struc_precip(LDT_rc%nnest) )
    allocate ( rdhm356_struc_temper(LDT_rc%nnest) )
    allocate ( const_wind(LDT_rc%nnest) ) 

    ! - Retrieve forcing configuration from lis.config
    call readcrd_rdhm356()

    do n=1, LDT_rc%nnest
       ! rdhm356_struc_precip(n)%ts = 3600, ts is read from lis configuration file by rdhm356_readcrd 
       call LDT_update_timestep(LDT_rc, n, rdhm356_struc_precip(n)%ts)
    enddo
    
    LDT_rc%met_ts(findex) = rdhm356_struc_precip(1)%ts
    LDT_rc%met_nf(findex) = 3
    LDT_rc%met_zterp(findex) = .false.
    LDT_rc%met_proj(findex)  = "polar"

    do n=1,LDT_rc%nnest

       rdhm356_struc_precip(n)%ncol = nint((rdhm356_struc_precip(n)%upper_right_hrapx - &
                                            rdhm356_struc_precip(n)%lower_left_hrapx) / &
                                            rdhm356_struc_precip(n)%hrap_resolution) + 1 
       rdhm356_struc_precip(n)%nrow = nint((rdhm356_struc_precip(n)%upper_right_hrapy - &
                                            rdhm356_struc_precip(n)%lower_left_hrapy) / &
                                            rdhm356_struc_precip(n)%hrap_resolution) + 1
       
       LDT_rc%met_nc(findex) = rdhm356_struc_precip(n)%ncol
       LDT_rc%met_nr(findex) = rdhm356_struc_precip(n)%nrow

       hrap_x_precip = rdhm356_struc_precip(n)%lower_left_hrapx
       hrap_y_precip = rdhm356_struc_precip(n)%lower_left_hrapy
        

       rdhm356_struc_temper(n)%ncol = nint((rdhm356_struc_temper(n)%upper_right_hrapx - &
                                            rdhm356_struc_temper(n)%lower_left_hrapx) / &
                                            rdhm356_struc_temper(n)%hrap_resolution) + 1 
       rdhm356_struc_temper(n)%nrow = nint((rdhm356_struc_temper(n)%upper_right_hrapy - &
                                            rdhm356_struc_temper(n)%lower_left_hrapy) / &
                                            rdhm356_struc_temper(n)%hrap_resolution) + 1
       hrap_x_temper = rdhm356_struc_temper(n)%lower_left_hrapx
       hrap_y_temper = rdhm356_struc_temper(n)%lower_left_hrapy


       call hrap_to_latlon(hrap_x_precip, hrap_y_precip, rlon_precip, rlat_precip)
       call hrap_to_latlon(hrap_x_temper, hrap_y_temper, rlon_temper, rlat_temper)
       
       call hrap_to_latlon(hrap_x_precip,  &
                           hrap_y_precip,  &
                           rlon_precip_ll, &
                           rlat_precip_ll)
       call hrap_to_latlon(hrap_x_temper,  &
                           hrap_y_temper,  &
                           rlon_temper_ll, &
                           rlat_temper_ll)

       gridDesciPrecip = 0
       gridDesciPrecip(n,1) = 8                            ! Projection type (hrap) 
       gridDesciPrecip(n,2) = rdhm356_struc_precip(n)%ncol ! X-dir amount of points
       gridDesciPrecip(n,3) = rdhm356_struc_precip(n)%nrow ! y-dir amount of points
       gridDesciPrecip(n,4) = rlat_precip_ll               ! Starting latitude point 
       gridDesciPrecip(n,5) = rlon_precip_ll               ! Starting longitude point 
       gridDesciPrecip(n,6) = 8 
       gridDesciPrecip(n,7) = 0                            ! Orientation (was 64) sw
       gridDesciPrecip(n,8) = 4.7625                       ! X-spacing length (kms) 
       gridDesciPrecip(n,9) = 4.7625                       ! Y-spacing length (kms)
       gridDesciPrecip(n,10) = 60                          ! True lat  
       gridDesciPrecip(n,11) = -105.0                      ! Standard longitude (???)
       gridDesciPrecip(n,13) = 0
       gridDesciPrecip(n,20) = 64.0

       LDT_rc%met_gridDesc(findex,1:20) = gridDesciPrecip(n,1:20)

       gridDesciTemper = 0
       gridDesciTemper(n,1) = 8                            ! Projection type (hrap)
       gridDesciTemper(n,2) = rdhm356_struc_temper(n)%ncol ! X-dir amount of points
       gridDesciTemper(n,3) = rdhm356_struc_temper(n)%nrow ! y-dir amount of points
       gridDesciTemper(n,4) = rlat_temper_ll               ! Starting latitude point
       gridDesciTemper(n,5) = rlon_temper_ll               ! Starting longitude point 
       gridDesciTemper(n,6) = 8 
       gridDesciTemper(n,7) = 0                            ! Orientation (was 64) sw 
       gridDesciTemper(n,8) = 4.7625                       ! X-spacing length (kms) 
       gridDesciTemper(n,9) = 4.7625                       ! Y-spacing length (kms)
       gridDesciTemper(n,10) = 60                          ! True lat  
       gridDesciTemper(n,11) = -105.0                      ! Standard longitude (???)
       gridDesciTemper(n,13) = 0
       gridDesciTemper(n,20) = 64.0


       rdhm356_struc_precip(n)%mi = rdhm356_struc_precip(n)%ncol * rdhm356_struc_precip(n)%nrow
       rdhm356_struc_temper(n)%mi = rdhm356_struc_temper(n)%ncol * rdhm356_struc_temper(n)%nrow


    !- If only processing parameters, then return to main routine calls ...
       if( LDT_rc%runmode == "LSM parameter processing" ) return

         ! === BILINEAR INTERPOLATION ==== 
       if (trim(LDT_rc%met_gridtransform(findex)) .eq. "bilinear") then
          allocate(rdhm356_struc_precip(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_precip(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_precip(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_precip(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_precip(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_precip(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_precip(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_precip(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n, gridDesciPrecip(n,:), &
               rdhm356_struc_precip(n)%n111, rdhm356_struc_precip(n)%n121, &
               rdhm356_struc_precip(n)%n211, rdhm356_struc_precip(n)%n221, &
               rdhm356_struc_precip(n)%w111, rdhm356_struc_precip(n)%w121, &
               rdhm356_struc_precip(n)%w211, rdhm356_struc_precip(n)%w221 )


          allocate(rdhm356_struc_temper(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_temper(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_temper(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_temper(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_temper(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_temper(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_temper(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(rdhm356_struc_temper(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n, gridDesciTemper(n,:), &
               rdhm356_struc_temper(n)%n111, rdhm356_struc_temper(n)%n121, &
               rdhm356_struc_temper(n)%n211, rdhm356_struc_temper(n)%n221, &
               rdhm356_struc_temper(n)%w111, rdhm356_struc_temper(n)%w121, &
               rdhm356_struc_temper(n)%w211, rdhm356_struc_temper(n)%w221 )

       elseif ( trim(LDT_rc%met_gridtransform(findex)) .eq. "budget-bilinear" ) then 
          ! === BUDGET BILINEAR INTERPOLATION ==== 
          allocate(rdhm356_struc_precip(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_precip(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_precip(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_precip(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_precip(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_precip(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_precip(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_precip(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input( n,gridDesciPrecip(n,:), &
                 rdhm356_struc_precip(n)%n112, rdhm356_struc_precip(n)%n122, &
                 rdhm356_struc_precip(n)%n212, rdhm356_struc_precip(n)%n222, &
                 rdhm356_struc_precip(n)%w112, rdhm356_struc_precip(n)%w122, &
                 rdhm356_struc_precip(n)%w212, rdhm356_struc_precip(n)%w222 )

          allocate(rdhm356_struc_temper(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_temper(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_temper(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_temper(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_temper(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_temper(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_temper(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(rdhm356_struc_temper(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input( n,gridDesciTemper(n,:), &
                 rdhm356_struc_temper(n)%n112, rdhm356_struc_temper(n)%n122, &
                 rdhm356_struc_temper(n)%n212, rdhm356_struc_temper(n)%n222, &
                 rdhm356_struc_temper(n)%w112, rdhm356_struc_temper(n)%w122, &
                 rdhm356_struc_temper(n)%w212, rdhm356_struc_temper(n)%w222 )

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then     ! Nearest Neighbor

          allocate(rdhm356_struc_precip(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesciPrecip(n,:),&
               rdhm356_struc_precip(n)%n113)

          allocate(rdhm356_struc_temper(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesciTemper(n,:),&
               rdhm356_struc_temper(n)%n113)

       end if   ! End interp option statement
    enddo
!    n = 1 
!    open(unit=1001, file='shit.txt', status='unknown');
!    do t=1, LDT_rc%lnc(n)*LDT_rc%lnr(n)
!        write(1001,'(2F14.6,4I8, 4F14.6)') rdhm356_struc_precip(n)%rlat1(t), rdhm356_struc_precip(n)%rlon1(t), &
!                      rdhm356_struc_precip(n)%n111(t), rdhm356_struc_precip(n)%n121(t), &
!                      rdhm356_struc_precip(n)%n211(t), rdhm356_struc_precip(n)%n221(t), &
!                      rdhm356_struc_precip(n)%w111(t), rdhm356_struc_precip(n)%w121(t), &
!                      rdhm356_struc_precip(n)%w211(t), rdhm356_struc_precip(n)%w221(t)
!    enddo
!    close(1001)
  end subroutine init_rdhm356

end module rdhm356_forcingMod
