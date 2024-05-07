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
!  The implementation in LDT has the derived data type {\tt agrradps\_struc}
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
!    value indicating the span of the LDT domain 
!    (1-NH only, 2-SH only, 3-span includes both NH and SH)
!  \item[mo1]
!    number of elements in the NH for the LDT grid
!  \item[mo2]
!    number of elements in the SH for the LDT grid
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
!  \item[n111\_nh,n121\_nh,n211\_nh,n221\_nh]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111\_nh,w121\_nh,w211\_nh,w221\_nh]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n111\_sh,n121\_sh,n211\_sh,n221\_sh]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111\_sh,w121\_sh,w211\_sh,w221\_sh]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[smask]
!    landmask to used to fill any gaps in bilinear interpolation
!    due to mismatches in the LDT and AGRMET masks
!  \item[fillflag1]
!    flag to check for filling gaps due to LDT and AGRMET
!    mask mismatches, for bilinear interpolation
!  \item[fillflag2]
!    flag to check for filling gaps due to LDT and AGRMET
!    mask mismatches, for neighbor interpolation
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
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
     character(len=LDT_CONST_PATH_LEN) :: agrpsdir
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
     integer, allocatable     :: n111_nh(:)
     integer, allocatable     :: n121_nh(:)
     integer, allocatable     :: n211_nh(:)
     integer, allocatable     :: n221_nh(:)
     real, allocatable        :: w111_nh(:),w121_nh(:)
     real, allocatable        :: w211_nh(:),w221_nh(:)
     integer, allocatable     :: n111_sh(:)
     integer, allocatable     :: n121_sh(:)
     integer, allocatable     :: n211_sh(:)
     integer, allocatable     :: n221_sh(:)
     real, allocatable        :: w111_sh(:),w121_sh(:)
     real, allocatable        :: w211_sh(:),w221_sh(:)
     integer, allocatable     :: smask(:,:)
     logical              :: fillflag1
     logical              :: fillflag2
     integer              :: findtime1, findtime2
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
    use LDT_coreMod,    only : LDT_rc
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
    use LDT_logMod,     only : LDT_logunit, LDT_endrun

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex

!
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for AGRRADPS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readagrradpscrd](\ref{readagrradpscrd}) \newline
!     reads the runtime options specified for AGRRADPS data
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP

    integer :: n 
    real    :: gridDesci(20)
    real    :: gridDesco(20)
    real, parameter :: xmeshl = 47.625
    real, parameter :: xpnmcaf = 257
    real, parameter :: ypnmcaf = 257
    real    :: xmesh, orient,xi1,xj1
    real    :: alat1,alon1
    integer :: ihemi
    integer :: updoy, yr1,mon1,da1,hr1,mn1,ss1
    real    :: upgmt
! _________________________________________

    allocate(agrradps_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing AGRMET Radiation forcing grid ... "

!    call readcrd_agrradps()

    agrradps_struc%imax = 512
    agrradps_struc%jmax = 512
    agrradps_struc%mi   = agrradps_struc(1)%imax* &
                          agrradps_struc(1)%jmax
    LDT_rc%met_nc(findex) = agrradps_struc(1)%imax
    LDT_rc%met_nr(findex) = agrradps_struc(1)%jmax

    LDT_rc%met_nf(findex) = 2
    LDT_rc%met_ts(findex) = 3*60*60

    do n = 1, LDT_rc%nnest

       agrradps_struc(n)%ts = 3*60*60 
       call LDT_update_timestep(LDT_rc, n, agrradps_struc(n)%ts)

       agrradps_struc(n)%findtime1 = 0 
       agrradps_struc(n)%findtime2 = 0
       agrradps_struc(n)%agrtime1 = 1000.0
       agrradps_struc(n)%agrtime2 = 0.0

       agrradps_struc(n)%interp = 1

! The grid span (whether it spans across the equator or is limited
! to either hemisphere is determined. This information is currently 
! determined if for a run domain defined in the lat/lon mode. 
! currently we can run on a regular LDT lat/lon or the AFWA lat/lon

       if( LDT_rc%gridDesc(n,1) .eq. 0 ) then    ! latlon domain
          allocate(agrradps_struc(n)%smask(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          agrradps_struc(n)%fillflag1 = .true.
          agrradps_struc(n)%fillflag2 = .true.
          if(LDT_rc%gridDesc(n,4).ge.0 .and. LDT_rc%gridDesc(n,7).ge.0) then
             agrradps_struc(n)%gridspan = 1 ! NH only
             agrradps_struc(n)%shemi = 1
             agrradps_struc(n)%nhemi = 1
          elseif(LDT_rc%gridDesc(n,4).le.0 .and. LDT_rc%gridDesc(n,7).le.0) then
             agrradps_struc(n)%gridspan = 2 ! SH only 
             agrradps_struc(n)%shemi = 2
             agrradps_struc(n)%nhemi = 2
          else
             agrradps_struc(n)%gridspan = 3 ! NH and SH
             agrradps_struc(n)%shemi = 1
             agrradps_struc(n)%nhemi = 2
          endif
          agrradps_struc(n)%hemi_nc = nint((LDT_rc%gridDesc(n,8)- &
               LDT_rc%gridDesc(n,5))/LDT_rc%gridDesc(n,9))+1
          if(agrradps_struc(n)%gridspan.eq.1) then
             agrradps_struc(n)%hemi_nr(1) = nint((LDT_rc%gridDesc(n,7)-&
                  LDT_rc%gridDesc(n,4))/LDT_rc%gridDesc(n,10))+1
             agrradps_struc(n)%hemi_nr(2) = 0
             agrradps_struc(n)%mo1 = agrradps_struc(n)%hemi_nc(1)*&
                  agrradps_struc(n)%hemi_nr(1)
             agrradps_struc(n)%mo2 = 0
          elseif(agrradps_struc(n)%gridspan.eq.2) then
             agrradps_struc(n)%hemi_nr(1) = 0
             agrradps_struc(n)%hemi_nr(2) = nint((LDT_rc%gridDesc(n,7)-&
                  LDT_rc%gridDesc(n,4))/LDT_rc%gridDesc(n,10)+1)
             agrradps_struc(n)%mo1 = 0
             agrradps_struc(n)%mo2 = agrradps_struc(n)%hemi_nc(2)*&
                  agrradps_struc(n)%hemi_nr(2)
          else
             agrradps_struc(n)%hemi_nr(1) = nint((LDT_rc%gridDesc(n,7)-&
                  LDT_rc%gridDesc(n,10)/2)/LDT_rc%gridDesc(n,10)+1)
             agrradps_struc(n)%hemi_nr(2) = nint((-LDT_rc%gridDesc(n,10)/2-&
                  LDT_rc%gridDesc(n,4))/LDT_rc%gridDesc(n,10)+1)
             agrradps_struc(n)%mo1 = agrradps_struc(n)%hemi_nc(1)*&
                  agrradps_struc(n)%hemi_nr(1)
             agrradps_struc(n)%mo2 = agrradps_struc(n)%hemi_nc(2)*&
                  agrradps_struc(n)%hemi_nr(2)
          endif

          gridDesco = 0
          do ihemi = agrradps_struc(n)%shemi, agrradps_struc(n)%nhemi

             gridDesco(1) = 0
             gridDesco(2) = agrradps_struc(n)%hemi_nc(ihemi)
             gridDesco(3) = agrradps_struc(n)%hemi_nr(ihemi)
           ! flip longitude sign
             gridDesco(5) = -1.0*LDT_rc%gridDesc(n,5)
             gridDesco(8) = -1.0*LDT_rc%gridDesc(n,8)
             gridDesco(6) = LDT_rc%gridDesc(n,6)
             gridDesco(9) = LDT_rc%gridDesc(n,9)
             gridDesco(10) = LDT_rc%gridDesc(n,10)
             gridDesco(20) = 255
             if(agrradps_struc(n)%gridspan.eq.1.or.&
                  agrradps_struc(n)%gridspan.eq.2) then
                gridDesco(4) = LDT_rc%gridDesc(n,4)
                gridDesco(7) = LDT_rc%gridDesc(n,7)
             else
                if(ihemi.eq.1) then
                   gridDesco(4) = LDT_rc%gridDesc(n,9)/2
                   gridDesco(7) = LDT_rc%gridDesc(n,7)
                else
                   gridDesco(4) = LDT_rc%gridDesc(n,4)
                   gridDesco(7) = -LDT_rc%gridDesc(n,9)/2
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

          !- AGRMET Radiation Grid description:
             LDT_rc%met_proj(findex)        = "polar"
             LDT_rc%met_gridDesc(findex,1)  = 5
             LDT_rc%met_gridDesc(findex,2)  = agrradps_struc(n)%imax
             LDT_rc%met_gridDesc(findex,3)  = agrradps_struc(n)%jmax
             LDT_rc%met_gridDesc(findex,4)  = alat1
             LDT_rc%met_gridDesc(findex,5)  = alon1
             LDT_rc%met_gridDesc(findex,6)  = 8
             LDT_rc%met_gridDesc(findex,7)  = orient
             LDT_rc%met_gridDesc(findex,8)  = xmesh
             LDT_rc%met_gridDesc(findex,9)  = xmesh
             LDT_rc%met_gridDesc(findex,10) = 0.0
             if(ihemi .eq.2) then
                LDT_rc%met_gridDesc(findex,10) = 128
             endif
             LDT_rc%met_gridDesc(findex,11) = 64
             LDT_rc%met_gridDesc(findex,13) = 1     ! global grid
             LDT_rc%met_gridDesc(findex,20) = 64
 
             gridDesci(:) = LDT_rc%met_gridDesc(findex,:)

#if 0
             if( ihemi.eq.1 ) then
               allocate(agrradps_struc(n)%n111_nh(agrradps_struc(n)%mo1))
               allocate(agrradps_struc(n)%n121_nh(agrradps_struc(n)%mo1))
               allocate(agrradps_struc(n)%n211_nh(agrradps_struc(n)%mo1))
               allocate(agrradps_struc(n)%n221_nh(agrradps_struc(n)%mo1))
               allocate(agrradps_struc(n)%w111_nh(agrradps_struc(n)%mo1))
               allocate(agrradps_struc(n)%w121_nh(agrradps_struc(n)%mo1))
               allocate(agrradps_struc(n)%w211_nh(agrradps_struc(n)%mo1))
               allocate(agrradps_struc(n)%w221_nh(agrradps_struc(n)%mo1))
               call bilinear_interp_input(n, gridDesci(:),&
                  agrradps_struc(n)%n111_nh,&
                  agrradps_struc(n)%n121_nh,agrradps_struc(n)%n211_nh,&
                  agrradps_struc(n)%n221_nh,agrradps_struc(n)%w111_nh,&
                  agrradps_struc(n)%w121_nh,agrradps_struc(n)%w211_nh,&
                  agrradps_struc(n)%w221_nh)

             elseif( ihemi.eq.2 ) then
               allocate(agrradps_struc(n)%n111_sh(agrradps_struc(n)%mo2))
               allocate(agrradps_struc(n)%n121_sh(agrradps_struc(n)%mo2))
               allocate(agrradps_struc(n)%n211_sh(agrradps_struc(n)%mo2))
               allocate(agrradps_struc(n)%n221_sh(agrradps_struc(n)%mo2))
               allocate(agrradps_struc(n)%w111_sh(agrradps_struc(n)%mo2))
               allocate(agrradps_struc(n)%w121_sh(agrradps_struc(n)%mo2))
               allocate(agrradps_struc(n)%w211_sh(agrradps_struc(n)%mo2))
               allocate(agrradps_struc(n)%w221_sh(agrradps_struc(n)%mo2))
               call bilinear_interp_input(n, gridDesci(:),&
                  agrradps_struc(n)%n111_sh,&
                  agrradps_struc(n)%n121_sh,agrradps_struc(n)%n211_sh,&
                  agrradps_struc(n)%n221_sh,agrradps_struc(n)%w111_sh,&
                  agrradps_struc(n)%w121_sh,agrradps_struc(n)%w211_sh,&
                  agrradps_struc(n)%w221_sh)
             endif
#endif
          enddo   ! end ihemi loop

       else
          write(LDT_logunit,*) 'ERROR: projection not supported here '
          call LDT_endrun
       endif
    enddo

  end subroutine init_agrradps

end module agrradps_forcingMod
