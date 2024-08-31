!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nldas2_forcingMod
!BOP
! !MODULE: nldas2_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data used in the North American
!  Land Data Assimilation System Phase II.  The variables are produced 
!  at 0.125 degree spatial resolution, and at hourly intervals.  For more
!  details please view the forcing files manual available at the 
!  following URL:
!
!  http://ldas.gsfc.nasa.gov//nldas/NLDAS2forcing.php
! 
!  The implemenatation in LDT has the derived data type {\tt nldas2\_struc}
!  that includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[nldas2time1]
!    The nearest, previous hourly instance of the incoming 
!    data (as a real time). 
!  \item[nldas2time2]
!    The nearest, next hourly instance of the incoming 
!    data (as a real time).
!  \item[nldas2dir]
!    Directory containing the input data
!  \item[nldas2\_filesrc]
!    Center(GES-DISC|NCEP)-based NLDAS-2 filename source option
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for nearest neighbor interpolation. 
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !REVISION HISTORY: 
! 02 Feb 2004: Sujay Kumar; Initial Specification
! 24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
! 14 Mar 2014: David Mocko: Added CAPE and PET forcing from NLDAS-2
! 
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_NLDAS2      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nldas2_struc
!EOP

  type, public ::  nldas2_type_dec 

     real          :: ts
     integer       :: nc, nr         ! AWIPS 212 dimensions
     character*50  :: nldas2_filesrc
     character(len=LDT_CONST_PATH_LEN) :: nldas2dir        ! NLDAS-2 Forcing Directory
     character(len=LDT_CONST_PATH_LEN) :: file_elevdiff
     character(len=LDT_CONST_PATH_LEN) :: file_narrelev

     real*8        :: nldas2time1,nldas2time2
     integer       :: findtime1, findtime2
     integer       :: model_level_data 
     integer       :: model_level_press 
     integer       :: model_pcp_data 
     integer       :: model_dswrf_data 
     
     real,  allocatable     :: orig_ediff(:)

     real                   :: gridDesc(20)
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

  end type nldas2_type_dec

  type(nldas2_type_dec), allocatable :: nldas2_struc(:)
!EOP

contains
  
!BOP
!
! !ROUTINE: init_NLDAS2
! \label{init_NLDAS2}
!
! !INTERFACE:
  subroutine init_NLDAS2(findex)
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_timeMgrMod, only : LDT_update_timestep
    use LDT_logMod,     only : LDT_logunit,LDT_endrun
!    use LDT_spatialDownscalingMod, only : LDT_init_pcpclimo_native
    use map_utils,      only : proj_latlon

    implicit none

    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NLDAS-2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_nldas2](\ref{readcrd_nldas2}) \newline
!     reads the runtime options specified for NLDAS-2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_nldas2\_elev](\ref{read_nldas2_elev}) \newline
!    reads the native elevation of the NLDAS-2 data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
    
    integer :: n
    real, allocatable :: elev(:,:)
    real, allocatable :: elevdiff(:,:)

    
    allocate(nldas2_struc(LDT_rc%nnest))

    write(unit=LDT_logunit,fmt=*)"[INFO] Initializing NLDAS-2 forcing grid ... "

! - Read LDT config NLDAS-2 entries:
    call readcrd_nldas2(findex)

    do n=1, LDT_rc%nnest
       nldas2_struc(n)%ts = 3600
       call LDT_update_timestep(LDT_rc, n, nldas2_struc(n)%ts)
    enddo

  ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex) = "latlon"
    LDT_rc%met_nc(findex) = 464
    LDT_rc%met_nr(findex) = 224

    nldas2_struc(:)%nc = 464
    nldas2_struc(:)%nr = 224

    nldas2_struc%gridDesc(1) = 0
    nldas2_struc%gridDesc(2) = nldas2_struc(1)%nc
    nldas2_struc%gridDesc(3) = nldas2_struc(1)%nr
    nldas2_struc%gridDesc(4) = 25.0625
    nldas2_struc%gridDesc(5) = -124.9375
    nldas2_struc%gridDesc(6) = 128
    nldas2_struc%gridDesc(7) = 52.9375
    nldas2_struc%gridDesc(8) = -67.0625
    nldas2_struc%gridDesc(9) = 0.125
    nldas2_struc%gridDesc(10) = 0.125
    nldas2_struc%gridDesc(20) = 64

    LDT_rc%met_gridDesc(findex,1:20) = nldas2_struc(1)%gridDesc(1:20)

    nldas2_struc%mi = nldas2_struc%nc*nldas2_struc%nr

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 13
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true.

    do n=1,LDT_rc%nnest

       nldas2_struc(n)%findtime1 = 0 
       nldas2_struc(n)%findtime2 = 0 

       if( nldas2_struc(n)%gridDesc(9)  == LDT_rc%gridDesc(n,9) .and. &
           nldas2_struc(n)%gridDesc(10) == LDT_rc%gridDesc(n,10).and. &
           LDT_rc%gridDesc(n,1) == proj_latlon .and. &
           LDT_rc%met_gridtransform(findex) .ne. "neighbor" ) then
         write(LDT_logunit,*) "[ERR]  The NLDAS-2 0.125 deg grid was selected for the"
         write(LDT_logunit,*) "  LDT run domain; however, 'bilinear', 'budget-bilinear',"
         write(LDT_logunit,*) "  or some other unknown option was selected to spatially"
         write(LDT_logunit,*) "  downscale the grid, which will cause errors during runtime."
         write(LDT_logunit,*) "Program stopping ..."
         call LDT_endrun()
       endif

     ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )
        case( "bilinear" )
          allocate(nldas2_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,&
               nldas2_struc(n)%gridDesc(:),&
               nldas2_struc(n)%n111,nldas2_struc(n)%n121,nldas2_struc(n)%n211,&
               nldas2_struc(n)%n221,nldas2_struc(n)%w111,nldas2_struc(n)%w121,&
               nldas2_struc(n)%w211,nldas2_struc(n)%w221)

        case( "budget-bilinear" )

          allocate(nldas2_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas2_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,nldas2_struc(n)%gridDesc(:),&
               nldas2_struc(n)%n111,nldas2_struc(n)%n121,&
               nldas2_struc(n)%n211,nldas2_struc(n)%n221,&
               nldas2_struc(n)%w111,nldas2_struc(n)%w121,&
               nldas2_struc(n)%w211,nldas2_struc(n)%w221)

          allocate(nldas2_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas2_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas2_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas2_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas2_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas2_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas2_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas2_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,nldas2_struc(n)%gridDesc(:),&
               nldas2_struc(n)%n112,nldas2_struc(n)%n122,&
               nldas2_struc(n)%n212,nldas2_struc(n)%n222,&
               nldas2_struc(n)%w112,nldas2_struc(n)%w122,&
               nldas2_struc(n)%w212,nldas2_struc(n)%w222)

        case( "neighbor" )
          allocate(nldas2_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          
          call neighbor_interp_input(n,nldas2_struc(n)%gridDesc(:),&
               nldas2_struc(n)%n113)

        case default
          write(LDT_logunit,*) 'Interpolation option not specified for NLDAS2'
          write(LDT_logunit,*) 'Program stopping...'
          call LDT_endrun()
       end select 

#if 0
       if(LDT_rc%met_ecor(findex).ne."none") then 
          allocate(nldas2_struc(n)%orig_ediff(&
                    nldas2_struc(n)%nc*nldas2_struc(n)%nr))

          allocate( elev(LDT_rc%lnc(n),LDT_rc%lnr(n)) )
          allocate( elevdiff(nldas2_struc(n)%nc,nldas2_struc(n)%nr) )

          call read_nldas2_elev(n, findex, elev, elevdiff)

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_domain(n)%gindex(c,r).ne.-1) then
                   LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = elev(c,r)
                   nldas2_struc(n)%orig_ediff(LDT_domain(n)%gindex(c,r)) = elevdiff(c,r)
                endif
             enddo
          enddo
          deallocate( elev, elevdiff )
       endif

       if(LDT_rc%pcp_downscale(findex).ne.0) then
          call LDT_init_pcpclimo_native(n,findex,&
               nint(nldas2_struc(n)%gridDesc(2)),&
          
       endif
#endif
    enddo

  end subroutine init_NLDAS2

end module nldas2_forcingMod
