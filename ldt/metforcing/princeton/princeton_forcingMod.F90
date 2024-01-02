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
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
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
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n12,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
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
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

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
     real               :: ts
     integer            :: nc, nr   !AWIPS 212 dimensions
     character(len=LDT_CONST_PATH_LEN)      :: princetondir
     character(len=LDT_CONST_PATH_LEN)      :: elevfile
     integer            :: mi
     integer            :: nmif
     real*8             :: princetontime1,princetontime2
     integer            :: findtime1, findtime2
     real               :: gridDesc(20)

     integer, allocatable :: n111(:)
     integer, allocatable :: n121(:)
     integer, allocatable :: n211(:)
     integer, allocatable :: n221(:)
     real, allocatable    :: w111(:),w121(:)
     real, allocatable    :: w211(:),w221(:)     
     integer, allocatable :: n112(:,:)
     integer, allocatable :: n122(:,:)
     integer, allocatable :: n212(:,:)
     integer, allocatable :: n222(:,:)
     real, allocatable    :: w112(:,:),w122(:,:)
     real, allocatable    :: w212(:,:),w222(:,:)
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
! 
! !INTERFACE:
  subroutine init_PRINCETON(findex)

! !USES: 
    use LDT_coreMod,    only : LDT_rc
    use LDT_timeMgrMod, only : LDT_update_timestep
    use LDT_logMod,     only : LDT_endrun, LDT_logunit
 
    implicit none
    integer, intent(in)     :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for PRINCETON 
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes \ref{interp}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readprincetoncrd](\ref{readprincetoncrd}) \newline
!     reads the runtime options specified for PRINCETON data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_princetonelev_ldtproc](\ref{read_princetonelev_ldtproc}) \newline
!    reads the native elevation of the princeton data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
   integer :: n 
   real, allocatable :: elev(:,:)
   real, allocatable :: dummy(:,:)

   allocate(princeton_struc(LDT_rc%nnest))

   write(LDT_logunit,fmt=*) "MSG: Initializing PRINCETON forcing grid ... "

! - Read LDT config NLDAS-1 entries:
    call readcrd_princeton(findex)

  ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex) = "latlon"
    LDT_rc%met_nc(findex) = 360
    LDT_rc%met_nr(findex) = 180
     
    princeton_struc%nc = 360
    princeton_struc%nr = 180

 !- PRINCETON Grid description:
    princeton_struc%gridDesc(1)  = 0
    princeton_struc%gridDesc(2)  = princeton_struc(1)%nc
    princeton_struc%gridDesc(3)  = princeton_struc(1)%nr
    princeton_struc%gridDesc(4)  = -89.50
    princeton_struc%gridDesc(5)  = -179.50
    princeton_struc%gridDesc(6)  = 128
    princeton_struc%gridDesc(7)  = 89.50
    princeton_struc%gridDesc(8)  = 179.50
    princeton_struc%gridDesc(9)  = 1.0
    princeton_struc%gridDesc(10) = 1.0
    princeton_struc%gridDesc(20) = 0

    princeton_struc%mi = princeton_struc%nc*princeton_struc%nr

    LDT_rc%met_gridDesc(findex,1:20) = princeton_struc(1)%gridDesc(1:20)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 9    ! Number of forcing fields
    LDT_rc%met_ts(findex) = 3*3600
    LDT_rc%met_zterp(findex) = .true.

    do n=1,LDT_rc%nnest

       princeton_struc(n)%ts = 3*3600
       call LDT_update_timestep(LDT_rc, n, princeton_struc(n)%ts)

     ! Setting up weights for Interpolation
       if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
          allocate(princeton_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n,  LDT_rc%met_gridDesc(findex,:), &
               princeton_struc(n)%n111,princeton_struc(n)%n121,&
               princeton_struc(n)%n211,princeton_struc(n)%n221,&
               princeton_struc(n)%w111,princeton_struc(n)%w121,&
               princeton_struc(n)%w211,princeton_struc(n)%w221)

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
          allocate(princeton_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(princeton_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n,  LDT_rc%met_gridDesc(findex,:),&
               princeton_struc(n)%n111,princeton_struc(n)%n121,&
               princeton_struc(n)%n211,princeton_struc(n)%n221,&
               princeton_struc(n)%w111,princeton_struc(n)%w121,&
               princeton_struc(n)%w211,princeton_struc(n)%w221)
          allocate(princeton_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(princeton_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(princeton_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(princeton_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(princeton_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(princeton_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(princeton_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(princeton_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n,  LDT_rc%met_gridDesc(findex,:),&
               princeton_struc(n)%n112,princeton_struc(n)%n122,&
               princeton_struc(n)%n212,princeton_struc(n)%n222,&
               princeton_struc(n)%w112,princeton_struc(n)%w122,&
               princeton_struc(n)%w212,princeton_struc(n)%w222)
       endif

#if 0
     ! Read PRINCETON Elevation file in:
       if(trim(LDT_rc%met_ecor(findex)).ne."none") then 
          allocate( elev(LDT_rc%lnc(n),LDT_rc%lnr(n)) )
          allocate( dummy(princeton_struc(n)%nc,princeton_struc(n)%nr) )
          call read_princeton_elev( n, findex, elev, dummy )

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_domain(n)%gindex(c,r).ne.-1) then
                   LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = elev(c,r)
                endif
             enddo
          enddo
          deallocate( elev, dummy )
       endif
#endif

    enddo   ! End nest loop

  end subroutine init_Princeton

end module princeton_forcingMod
