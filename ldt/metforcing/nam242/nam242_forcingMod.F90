!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nam242_forcingMod
!BOP
! !MODULE: nam242_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Data
!  Assimilation System (NAM) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. NAM forcing variables are produced
!  on a quadratic gaussian grid. LDT uses the 00, 03, 06 and as 
!  needed, the 09 forecasts. The forecasts are produced at 6 hr intervals. 
!
!  The implementation in LDT has the derived data type {\tt nam242\_struc} that
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the NAM data
!  \item[ts]
!    Frequency in seconds of the forcing data
!  \item[namtime1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[namtime2]
!    The nearest, next 3 hour instance of the incoming 
!    data (as a real time).
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for 
!    temporal interpolation.
!  \item[namdir]
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
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \end{description}
!
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_nam242    !defines the native resolution of 
                           !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nam242_struc
!EOP

  type, public        :: nam242_type_dec
     integer          :: nc, nr
     integer          :: nmif
     character(len=LDT_CONST_PATH_LEN)    :: namdir   !NAM Forcing Directory
     character(len=LDT_CONST_PATH_LEN)    :: elevfile
     real             :: ts
     real*8           :: namtime1,namtime2
     integer          :: findtime1,findtime2
     integer          :: mi
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
  end type nam242_type_dec
  
  type(nam242_type_dec), allocatable :: nam242_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_nam242
!  \label{init_nam242}
!
! !REVISION HISTORY: 
!     Sep 2012: NOHRSC/NOAA: Initial specification
! 
! !INTERFACE:
  subroutine init_nam242(findex)
! !USES: 
   use LDT_coreMod,    only : LDT_rc, LDT_domain
   use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
   use LDT_logMod,     only : LDT_logunit,LDT_endrun

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex

! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NAM
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}). Based on the NAM data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the NAM resolution switches are also 
!  defined in this routine. 
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_nam242](\ref{readcrd_nam242}) \newline
!     reads the runtime options specified for NAM data
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP

    integer :: n
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1,tdoy
    real    :: upgmt, tgmt
    real    :: gridDesci(20)
    real, allocatable :: elev(:,:)
    real, allocatable :: dummy(:,:)

    allocate(nam242_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing NAM-242 forcing grid ... "

    call readcrd_nam242(findex)

    do n=1, LDT_rc%nnest
       nam242_struc(n)%ts = 3*60*60 
       call LDT_update_timestep(LDT_rc, n, nam242_struc(n)%ts)

     ! TEMPORARY FIX: NAM242 only runs correctly when model ts < 3hr
     ! Set local - 1 hour timestep (to replicate model timestep):
       call LDT_update_timestep(LDT_rc, n, 3600.)
    enddo

    LDT_rc%met_proj(findex) = "polar"

    nam242_struc%nc = 553
    nam242_struc%nr = 425
    LDT_rc%met_nc(findex) = nam242_struc(1)%nc
    LDT_rc%met_nr(findex) = nam242_struc(1)%nr

    gridDesci(1) = 5
    gridDesci(2) = nam242_struc(1)%nc
    gridDesci(3) = nam242_struc(1)%nr
    gridDesci(4) = 30
    gridDesci(5) = -173
    gridDesci(6) = 0       ! not used
    gridDesci(7) = -135    ! Dagang question
    gridDesci(8) = 11.25
    gridDesci(9) = 11.25
    gridDesci(10) = 60     ! Dagang question
    gridDesci(11) = -135   ! Dagang question
    gridDesci(20) = 0      ! Dagang question

    LDT_rc%met_gridDesc(findex,:) = gridDesci(:)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    nam242_struc(:)%nmif  = 9
    LDT_rc%met_nf(findex) = 9   ! Number of forcing variables
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true.

    do n=1,LDT_rc%nnest

       nam242_struc(n)%mi = nam242_struc(n)%nc*nam242_struc(n)%nr
       
     ! Setting up weights for Interpolation
       if(LDT_rc%met_gridtransform(findex).eq."bilinear") then 
          allocate(nam242_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(:),&
               nam242_struc(n)%n111,nam242_struc(n)%n121,&
               nam242_struc(n)%n211,nam242_struc(n)%n221,&
               nam242_struc(n)%w111,nam242_struc(n)%w121,&
               nam242_struc(n)%w211,nam242_struc(n)%w221)

       elseif(LDT_rc%met_gridtransform(findex).eq."budget-bilinear") then 
          allocate(nam242_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nam242_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(:),&
               nam242_struc(n)%n111,nam242_struc(n)%n121,&
               nam242_struc(n)%n211,nam242_struc(n)%n221,&
               nam242_struc(n)%w111,nam242_struc(n)%w121,&
               nam242_struc(n)%w211,nam242_struc(n)%w221)

          allocate(nam242_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nam242_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nam242_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nam242_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nam242_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nam242_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nam242_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nam242_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(:),&
               nam242_struc(n)%n112,nam242_struc(n)%n122,&
               nam242_struc(n)%n212,nam242_struc(n)%n222,&
               nam242_struc(n)%w112,nam242_struc(n)%w122,&
               nam242_struc(n)%w212,nam242_struc(n)%w222)
       endif

!       if ( LDT_rc%met_ecor(findex).eq."lapse-rate" ) then 
!          call read_nam242_elev(n, findex, 1)
!       endif
#if 0
     ! Read NAM242 Elevation file in:
       if(trim(LDT_rc%met_ecor(findex)).ne."none") then
          allocate( elev(LDT_rc%lnc(n),LDT_rc%lnr(n)) )
          allocate( dummy(nam242_struc(n)%nc,nam242_struc(n)%nr) )
          call read_nam242_elev( n, findex, elev, dummy )

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

    enddo
  end subroutine init_nam242

end module nam242_forcingMod

