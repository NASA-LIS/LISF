!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module agrradps_forcingMod
!BOP
! !MODULE: agrradps_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the AGRMET radiation data.
!  The AGRMET data is in polar stereographic projection at 48km resolution
!  in two hemispheres, NH and SH at hourly.
! 
!  The implementation in LIS has the derived data type {\tt agrradps\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[agrtime1]
!    The nearest, previous hourly instance of the incoming 
!    data (as a real time).
!  \item[agrtime2]
!    The nearest, next hourly instance of the incoming 
!    data (as a real time).
!  \item[agrpsdir]
!    Directory containing the input data
!  \item[gridspan]
!    value indicating the span of the LIS domain 
!    (1-NH only, 2-SH only, 3-span includes both NH and SH)
!  \item[mo1]
!    number of elements in the NH for the LIS grid
!  \item[mo2]
!    number of elements in the SH for the LIS grid
!  \item[hemi\_nc]
!    number of columns (in the east-west dimension) for each
!    hemisphere
!  \item[hemi\_nr]
!    number of rows (in the north-south dimension) for each
!    hemisphere   
!  \item[shemi]
!    value indicating the starting index of hemisphere loops
!    (1-NH only,2-SH only, 1-span includes both NH and SH)
!  \item[nhemi]
!    value indicating the ending index of hemisphere loops
!    (1-NH only,2-SH only, 2-span includes both NH and SH)
!  \item[mi]
!    number of elements in the input grid
!  \item[imax]
!    x-dimension size of the input grid
!  \item[jmax]
!    y-dimension size of the input grid
!  \item[interp]
!    variable specifying if spatial interpolation needs to be
!    done
!  \item[rlat1\_nh]
!    Array containing the latitudes of the input grid for each corresponding
!    grid point in LIS (to be used for bilinear interpolation)
!  \item[rlon1\_nh]
!    Array containing the longitudes of the input grid for each corresponding
!    grid point in LIS (to be used for bilinear interpolation)
!  \item[n111\_nh,n121\_nh,n211\_nh,n221\_nh]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111\_nh,w121\_nh,w211\_nh,w221\_nh]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[rlat1\_sh]
!    Array containing the latitudes of the input grid for each corresponding
!    grid point in LIS (to be used for bilinear interpolation)
!  \item[rlon1\_sh]
!    Array containing the longitudes of the input grid for each corresponding
!    grid point in LIS (to be used for bilinear interpolation)
!  \item[n111\_sh,n121\_sh,n211\_sh,n221\_sh]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111\_sh,w121\_sh,w211\_sh,w221\_sh]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[smask]
!    landmask to used to fill any gaps in bilinear interpolation
!    due to mismatches in the LIS and AGRMET masks
!  \item[fillflag1]
!    flag to check for filling gaps due to LIS and AGRMET
!    mask mismatches, for bilinear interpolation
!  \item[fillflag2]
!    flag to check for filling gaps due to LIS and AGRMET
!    mask mismatches, for neighbor interpolation
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_agrradps    !defines the native resolution of 
                                    !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: agrradps_struc

!EOP
 
  type agrradps_type_dec
     real                 :: ts
     real*8               :: agrtime1,agrtime2
     character(len=LIS_CONST_PATH_LEN) :: agrpsdir
     integer              :: gridspan
     integer              :: mo1
     integer              :: mo2
     integer              :: hemi_nc(2)
     integer              :: hemi_nr(2)
     integer              :: shemi
     integer              :: nhemi
     integer              :: mi
     integer              :: imax
     integer              :: jmax
     integer              :: interp
     real, allocatable        :: rlat1_nh(:)
     real, allocatable        :: rlon1_nh(:)
     integer, allocatable     :: n111_nh(:)
     integer, allocatable     :: n121_nh(:)
     integer, allocatable     :: n211_nh(:)
     integer, allocatable     :: n221_nh(:)
     real, allocatable        :: w111_nh(:),w121_nh(:)
     real, allocatable        :: w211_nh(:),w221_nh(:)
     real, allocatable        :: rlat1_sh(:)
     real, allocatable        :: rlon1_sh(:)
     integer, allocatable     :: n111_sh(:)
     integer, allocatable     :: n121_sh(:)
     integer, allocatable     :: n211_sh(:)
     integer, allocatable     :: n221_sh(:)
     real, allocatable        :: w111_sh(:),w121_sh(:)
     real, allocatable        :: w211_sh(:),w221_sh(:)
     integer, allocatable     :: smask(:,:)
     logical              :: fillflag1
     logical              :: fillflag2
     integer              :: ncold
     integer              :: nrold  
     integer              :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type agrradps_type_dec

  type(agrradps_type_dec), allocatable :: agrradps_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_agrradps
! \label{init_agrradps}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_agrradps(findex)
! !USES: 
    use LIS_coreMod, only: LIS_rc
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only :LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex

!
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for AGRRADPS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_agrradps](\ref{readcrd_agrradps}) \newline
!     reads the runtime options specified for AGRRADPS data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP


    real :: gridDesci(50)
    real :: gridDesco(50)
    real, parameter :: xmeshl = 47.625
    real, parameter :: xpnmcaf = 257
    real, parameter :: ypnmcaf = 257
    real :: xmesh, orient,xi1,xj1
    real :: alat1,alon1
    integer :: ihemi
    integer :: updoy, yr1,mon1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n 


    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the AGRMET radiation forcing'
       write(LIS_logunit,*) '[ERR]  reader is not set up to run in forecast'
       write(LIS_logunit,*) '[ERR]  mode.  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(agrradps_struc(LIS_rc%nnest))
    call readcrd_agrradps()

    do n=1, LIS_rc%nnest
       agrradps_struc(n)%ts = 3*60*60 
       call LIS_update_timestep(LIS_rc, n, agrradps_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 2
    do n=1,LIS_rc%nnest

       allocate(agrradps_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(agrradps_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       agrradps_struc(n)%metdata1 = 0
       agrradps_struc(n)%metdata2 = 0

       agrradps_struc(n)%findtime1 = 0 
       agrradps_struc(n)%findtime2 = 0
!------------------------------------------------------------------------ 
! AGRMET product 
!------------------------------------------------------------------------
       agrradps_struc(n)%agrtime1 = 1000.0
       agrradps_struc(n)%agrtime2 = 0.0

       agrradps_struc(n)%interp = 1
       agrradps_struc(n)%imax = 512
       agrradps_struc(n)%jmax = 512
       agrradps_struc(n)%mi = agrradps_struc(n)%imax* &
                                     agrradps_struc(n)%jmax

! The grid span (whether it spans across the equator or is limited
! to either hemisphere is determined. This information is currently 
! determined if for a run domain defined in the lat/lon mode. 
! currently we can run on a regular LIS lat/lon or the AFWA lat/lon

       if(LIS_rc%gridDesc(n,1) .eq.0) then !latlon domain
          allocate(agrradps_struc(n)%smask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          agrradps_struc(n)%fillflag1 = .true.
          agrradps_struc(n)%fillflag2 = .true.
          if(LIS_rc%gridDesc(n,4).ge.0.and.LIS_rc%gridDesc(n,7).ge.0) then
             agrradps_struc(n)%gridspan = 1 ! NH only
             agrradps_struc(n)%shemi = 1
             agrradps_struc(n)%nhemi = 1
          elseif(LIS_rc%gridDesc(n,4).le.0.and.LIS_rc%gridDesc(n,7).le.0) then
             agrradps_struc(n)%gridspan = 2 ! SH only 
             agrradps_struc(n)%shemi = 2
             agrradps_struc(n)%nhemi = 2
          else
             agrradps_struc(n)%gridspan = 3 ! NH and SH
             agrradps_struc(n)%shemi = 1
             agrradps_struc(n)%nhemi = 2
          endif
          agrradps_struc(n)%hemi_nc = nint((LIS_rc%gridDesc(n,8)- &
               LIS_rc%gridDesc(n,5))/LIS_rc%gridDesc(n,9))+1
          if(agrradps_struc(n)%gridspan.eq.1) then
             agrradps_struc(n)%hemi_nr(1) = nint((LIS_rc%gridDesc(n,7)-&
                  LIS_rc%gridDesc(n,4))/LIS_rc%gridDesc(n,10))+1
             agrradps_struc(n)%hemi_nr(2) = 0
             agrradps_struc(n)%mo1 = agrradps_struc(n)%hemi_nc(1)*&
                  agrradps_struc(n)%hemi_nr(1)
             agrradps_struc(n)%mo2 = 0
          elseif(agrradps_struc(n)%gridspan.eq.2) then
             agrradps_struc(n)%hemi_nr(1) = 0
             agrradps_struc(n)%hemi_nr(2) = nint((LIS_rc%gridDesc(n,7)-&
                  LIS_rc%gridDesc(n,4))/LIS_rc%gridDesc(n,10)+1)
             agrradps_struc(n)%mo1 = 0
             agrradps_struc(n)%mo2 = agrradps_struc(n)%hemi_nc(2)*&
                  agrradps_struc(n)%hemi_nr(2)
          else
             agrradps_struc(n)%hemi_nr(1) = nint((LIS_rc%gridDesc(n,7)-&
                  LIS_rc%gridDesc(n,10)/2)/LIS_rc%gridDesc(n,10)+1)
             agrradps_struc(n)%hemi_nr(2) = nint((-LIS_rc%gridDesc(n,10)/2-&
                  LIS_rc%gridDesc(n,4))/LIS_rc%gridDesc(n,10)+1)
             agrradps_struc(n)%mo1 = agrradps_struc(n)%hemi_nc(1)*&
                  agrradps_struc(n)%hemi_nr(1)
             agrradps_struc(n)%mo2 = agrradps_struc(n)%hemi_nc(2)*&
                  agrradps_struc(n)%hemi_nr(2)
          endif

          gridDesco = 0
          gridDesci = 0
          do ihemi = agrradps_struc(n)%shemi,agrradps_struc(n)%nhemi
             gridDesco(1) = 0
             gridDesco(2) = agrradps_struc(n)%hemi_nc(ihemi)
             gridDesco(3) = agrradps_struc(n)%hemi_nr(ihemi)
!             gridDesco(5) = LIS_rc%gridDesc(n,5)
!             gridDesco(8) = LIS_rc%gridDesc(n,8)
! flip longitude sign
             gridDesco(5) = -1.0*LIS_rc%gridDesc(n,5)
             gridDesco(8) = -1.0*LIS_rc%gridDesc(n,8)
             gridDesco(6) = LIS_rc%gridDesc(n,6)
             gridDesco(9) = LIS_rc%gridDesc(n,9)
             gridDesco(10) = LIS_rc%gridDesc(n,10)
             gridDesco(20) = 255
             if(agrradps_struc(n)%gridspan.eq.1.or.&
                  agrradps_struc(n)%gridspan.eq.2) then
                gridDesco(4) = LIS_rc%gridDesc(n,4)
                gridDesco(7) = LIS_rc%gridDesc(n,7)
             else
                if(ihemi.eq.1) then
                   gridDesco(4) = LIS_rc%gridDesc(n,9)/2
                   gridDesco(7) = LIS_rc%gridDesc(n,7)
                else
                   gridDesco(4) = LIS_rc%gridDesc(n,4)
                   gridDesco(7) = -LIS_rc%gridDesc(n,9)/2
                endif
             endif

             if(ihemi.eq.1) then
                xmesh = xmeshl
                orient = 80.0
             else
                xmesh = -xmeshl
                orient = 260.0
             endif
             xj1 = float(1)-ypnmcaf
             xi1 = float(1)-xpnmcaf

             call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)

             gridDesci = 0
             gridDesci(1) = 5
             gridDesci(2) = agrradps_struc(n)%imax
             gridDesci(3) = agrradps_struc(n)%jmax
             gridDesci(4) = alat1
             gridDesci(5) = alon1
             gridDesci(6) = 8
             gridDesci(7) = orient
             gridDesci(8) = xmesh
             gridDesci(9) = xmesh
             gridDesci(10) = 0.0
             if(ihemi .eq.2) then
                gridDesci(10) = 128
             endif
             gridDesci(11) = 64
             gridDesci(13) = 1  !global grid
             gridDesci(20) = 0

             if(ihemi.eq.1) then
             allocate(agrradps_struc(n)%rlat1_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%rlon1_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%n111_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%n121_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%n211_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%n221_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%w111_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%w121_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%w211_nh(agrradps_struc(n)%mo1))
             allocate(agrradps_struc(n)%w221_nh(agrradps_struc(n)%mo1))
              call bilinear_interp_input_withgrid(gridDesci,gridDesco,&
                  agrradps_struc(n)%mo1,agrradps_struc(n)%rlat1_nh,&
                  agrradps_struc(n)%rlon1_nh,agrradps_struc(n)%n111_nh,&
                  agrradps_struc(n)%n121_nh,agrradps_struc(n)%n211_nh,&
                  agrradps_struc(n)%n221_nh,agrradps_struc(n)%w111_nh,&
                  agrradps_struc(n)%w121_nh,agrradps_struc(n)%w211_nh,&
                  agrradps_struc(n)%w221_nh)
             elseif(ihemi.eq.2) then
             allocate(agrradps_struc(n)%rlat1_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%rlon1_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%n111_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%n121_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%n211_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%n221_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%w111_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%w121_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%w211_sh(agrradps_struc(n)%mo2))
             allocate(agrradps_struc(n)%w221_sh(agrradps_struc(n)%mo2))
              call bilinear_interp_input_withgrid(gridDesci,gridDesco,&
                  agrradps_struc(n)%mo2,agrradps_struc(n)%rlat1_sh,&
                  agrradps_struc(n)%rlon1_sh,agrradps_struc(n)%n111_sh,&
                  agrradps_struc(n)%n121_sh,agrradps_struc(n)%n211_sh,&
                  agrradps_struc(n)%n221_sh,agrradps_struc(n)%w111_sh,&
                  agrradps_struc(n)%w121_sh,agrradps_struc(n)%w211_sh,&
                  agrradps_struc(n)%w221_sh)
             endif
          enddo
       else
           write(LIS_logunit,*) 'ERROR: projection not supported here '
           call LIS_endrun
       endif
    enddo
  end subroutine init_agrradps

end module agrradps_forcingMod
