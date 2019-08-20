!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nldas1_forcingMod
!BOP
! !MODULE: nldas1_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data used in the North American
!  Land Data Assimilation System (NLDAS; Cosgrove et al.(2003)). 
!  The variables are produced at 0.125 degree spatial resolution, 
!  and at hourly intervals, over the continental United States. 
!  The SW$\downarrow$ is generated from a 1/2$^\circ$ product
!  derived at the University of Maryland from NOAA's Geostationary Operational
!  Environmental Satellites (GOES). The LW$\downarrow$ is 
!  derived from 3 hourly NCEP Eta Data Assimilation System (EDAS) output
!  fields~\cite{rogers}, and from 3 hourly and 6 hourly Eta mesoscale 
!  model forecast fields when EDAS data is unavailable. 
!
!  This dataset represents NLDAS Phase 1 (NLDAS-1); for more see:
!       http://ldas.gsfc.nasa.gov/nldas/NLDAS1forcing.php
!
!  Cosgrove, B. Real-time and retrospective forcing in the North 
!  American Land Data Assimilation System (NLDAS) project. Journal of 
!  Geophysical Research, 108 (D22), 8842, DOI:10.1029/2002JD003118
! 
!  The implemenatation in LDT has the derived data type {\tt nldas1\_struc}
!  that includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[nldas1time1]
!    The nearest, previous hourly instance of the incoming 
!    data (as a real time). 
!  \item[nldas1time2]
!    The nearest, next hourly instance of the incoming 
!    data (as a real time).
!  \item[nldas1dir]
!    Directory containing the input data
!  \item[nldas1\_filesrc]
!    Center(GES-DISC|NCEP)-based NLDAS-1 filename source option
!  \item[elevfile]
!    File with the elevation definition for the input data. 
!  \item[mi]
!    Number of points in the input grid
!  \item[orig\_ediff]
!    Original NLDAS-1 elevation difference field 
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
!    Array containing the neighbor information of the input grid 
!    for each grid point in LDT, for nearest neighbor interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !REVISION HISTORY: 
!  02 Feb 2004: Sujay Kumar; Initial Specification
!  20 Oct 2007: Kristi Arsenault; Changed to EDAS elev file for new diff
!  15 Feb 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions
!  14 Mar 2014: David Mocko: Added NLDAS-1 precipitation and radiation
!                               field choices to the ldt.config file
!                            Changed flag from "1"/"2" to "NCEP"/"GES-DISC".
!
! !USES: 
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_nldas1      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nldas1_struc
!EOP

  type, public ::  nldas1_type_dec 
     real          :: ts
     integer       :: nc, nr           ! AWIPS 212 dimensions
     integer       :: mi
     character*50  :: nldas1_filesrc
     character*100 :: nldas1dir        ! NLDAS-1 Forcing Directory
     character*100 :: file_elevdiff    ! Original Elevdiff File
     character*100 :: file_edaselev
     character*40  :: prec_field,swdn_field
     real*8        :: nldas1time1,nldas1time2
     integer       :: findtime1, findtime2
     real          :: gridDesc(20)

     integer            :: applymask
     real, allocatable  :: conusmask(:,:)
     character*100      :: conusmaskfile
     real               :: mask_gd(6)

     real               :: edas_gridDesc(6)
     real, allocatable  :: orig_ediff(:)

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
  end type nldas1_type_dec

  type(nldas1_type_dec), allocatable :: nldas1_struc(:)
!EOP

contains
  
!BOP
!
! !ROUTINE: init_NLDAS1
! \label{init_NLDAS1} 
! 
! !INTERFACE:
  subroutine init_NLDAS1(findex)
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_timeMgrMod, only : LDT_update_timestep
    use LDT_logMod,     only : LDT_logunit, LDT_getNextUnitNumber, &
         LDT_releaseUnitNumber, LDT_endrun
    use map_utils,      only : ij_to_latlon, proj_latlon

    implicit none
    
    integer, intent(in) :: findex
    
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NLDAS-1
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_nldas1](\ref{readcrd_nldas1}) \newline
!     reads the runtime options specified for NLDAS-1 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!    computes the neighbor, weights for nearest neighbor interpolation
!   \item[read\_nldas1\_elev](\ref{read_nldas1_elev}) \newline
!    reads the native EDAS elevation and elevation difference data to be
!    used for topographic adjustments to the metforcing data.
!  \end{description}
!EOP
    
    real, allocatable :: rlat(:,:), rlon(:,:)
    integer :: n
    integer :: ftn
    integer :: nc_dom
    integer :: line1, line2, line
    integer :: c,r
    real, allocatable :: elev(:,:)
    real, allocatable :: elevdiff(:,:)
    
    allocate(nldas1_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing NLDAS-1 forcing grid ... "

! - Read LDT config NLDAS-1 entries:
    call readcrd_nldas1( findex )

    do n=1, LDT_rc%nnest
       nldas1_struc(n)%ts = 3600
       call LDT_update_timestep(LDT_rc, n, nldas1_struc(n)%ts)
    enddo

  ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex) = "latlon"
    LDT_rc%met_nc(findex) = 464
    LDT_rc%met_nr(findex) = 224

    nldas1_struc(:)%nc = 464
    nldas1_struc(:)%nr = 224

 !- NLDAS-1 Grid description:
    nldas1_struc(:)%gridDesc(1) = 0
    nldas1_struc(:)%gridDesc(2) = nldas1_struc(1)%nc
    nldas1_struc(:)%gridDesc(3) = nldas1_struc(1)%nr
    nldas1_struc(:)%gridDesc(4) = 25.0625
    nldas1_struc(:)%gridDesc(5) = -124.9375
    nldas1_struc(:)%gridDesc(6) = 128
    nldas1_struc(:)%gridDesc(7) = 52.9375
    nldas1_struc(:)%gridDesc(8) = -67.0625
    nldas1_struc(:)%gridDesc(9) = 0.125
    nldas1_struc(:)%gridDesc(10) = 0.125
    nldas1_struc(:)%gridDesc(20) = 64

    LDT_rc%met_gridDesc(findex,1:20) = nldas1_struc(1)%gridDesc(1:20)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 9    ! Number of forcing fields
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true.

    do n=1,LDT_rc%nnest

       if( nldas1_struc(n)%gridDesc(9) == LDT_rc%gridDesc(n,9) .and. &
           nldas1_struc(n)%gridDesc(10) == LDT_rc%gridDesc(n,10).and. &
           LDT_rc%gridDesc(n,1) == proj_latlon .and. &
           LDT_rc%met_gridtransform(findex) .ne. "neighbor" ) then
         write(LDT_logunit,*) "[WARN] The NLDAS-1 0.125 deg grid was selected for the"
         write(LDT_logunit,*) "  LDT run domain; however, 'bilinear', 'budget-bilinear',"
         write(LDT_logunit,*) "  or some other unknown option was selected to spatially"
         write(LDT_logunit,*) "  downscale the grid, which will cause errors during runtime."
         write(LDT_logunit,*) "Program stopping ..."
         call LDT_endrun()
       endif

       nldas1_struc(n)%mi = nldas1_struc(n)%nc*nldas1_struc(n)%nr

     ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )

        case( "bilinear" )
          allocate(nldas1_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,nldas1_struc(n)%gridDesc(:),&
               nldas1_struc(n)%n111,nldas1_struc(n)%n121,nldas1_struc(n)%n211,&
               nldas1_struc(n)%n221,nldas1_struc(n)%w111,nldas1_struc(n)%w121,&
               nldas1_struc(n)%w211,nldas1_struc(n)%w221)

        case( "budget-bilinear" )
          allocate(nldas1_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(nldas1_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n,nldas1_struc(n)%gridDesc(:),&
               nldas1_struc(n)%n111,nldas1_struc(n)%n121,&
               nldas1_struc(n)%n211,nldas1_struc(n)%n221,&
               nldas1_struc(n)%w111,nldas1_struc(n)%w121,&
               nldas1_struc(n)%w211,nldas1_struc(n)%w221)

          allocate(nldas1_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas1_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas1_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas1_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n,nldas1_struc(n)%gridDesc(:),&
               nldas1_struc(n)%n112,nldas1_struc(n)%n122,&
               nldas1_struc(n)%n212,nldas1_struc(n)%n222,&
               nldas1_struc(n)%w112,nldas1_struc(n)%w122,&
               nldas1_struc(n)%w212,nldas1_struc(n)%w222)

        case( "neighbor" )
          allocate(nldas1_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n,nldas1_struc(n)%gridDesc(:),&
               nldas1_struc(n)%n113)

        case default
          write(LDT_logunit,*) 'Interpolation option not specified for NLDAS1'
          write(LDT_logunit,*) 'Program stopping...'
          call LDT_endrun()
       end select

#if 0
    !- If Elevation Correction desired, remove existing NLDAS-1 1/8 degree
    !   elevation correction that is currently built into forcing:

       if ( LDT_rc%met_ecor(findex) .ne."none" ) then 

       !- Read in original NLDAS-1 EDAS-GTOPO Elevation Difference file::
          allocate( nldas1_struc(n)%orig_ediff(&
                     nldas1_struc(n)%nc*nldas1_struc(n)%nr))
          allocate( elev(LDT_rc%lnc(n),LDT_rc%lnr(n)) )
          allocate( elevdiff(nldas1_struc(n)%nc,nldas1_struc(n)%nr) )

          call read_nldas1_elev(n, findex, elev, elevdiff)

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_domain(n)%gindex(c,r).ne.-1) then
                   LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = elev(c,r)
                   nldas1_struc(n)%orig_ediff(LDT_domain(n)%gindex(c,r)) = elevdiff(c,r)
                endif
             enddo
          enddo
          deallocate( elev, elevdiff )
       endif
       
     ! Apply NLDAS-CONUS Landmask:
       if( nldas1_struc(n)%applymask == 1 ) then 

          allocate(nldas1_struc(n)%conusmask(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(rlat(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          allocate(rlon(LDT_rc%lnc(n),LDT_rc%lnr(n)))
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                call ij_to_latlon(LDT_domain(n)%ldtproj,float(c),float(r),&
                     rlat(c,r),rlon(c,r))                
             enddo
          enddo
          
          write(LDT_logunit,*) 'Reading the CONUS mask: ',&
               trim(nldas1_struc(n)%conusmaskfile)

          ftn = LDT_getNextUnitNumber()
          open(ftn,file=trim(nldas1_struc(n)%conusmaskfile),form='unformatted',&
               recl=4,access='direct',status='old')

          nc_dom = nint((nldas1_struc(n)%mask_gd(4)-nldas1_struc(n)%mask_gd(2))/&
                         nldas1_struc(n)%mask_gd(5))+1
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                line1 = nint((rlat(c,r)-nldas1_struc(n)%mask_gd(1))/&
                     nldas1_struc(n)%mask_gd(6))+1
                line2 = nint((rlon(c,r)-nldas1_struc(n)%mask_gd(2))/&
                     nldas1_struc(n)%mask_gd(5))+1
                line = (line1-1)*nc_dom+line2
                if(rlat(c,r).ge.nldas1_struc(n)%mask_gd(1).and.&
                     rlat(c,r).le.nldas1_struc(n)%mask_gd(3).and.&
                     rlon(c,r).ge.nldas1_struc(n)%mask_gd(2).and.&
                     rlon(c,r).le.nldas1_struc(n)%mask_gd(4)) then 
                   read(ftn,rec=line) nldas1_struc(n)%conusmask(c,r)
                else
                   nldas1_struc(n)%conusmask(c,r) = 0.0
                endif
             enddo
          enddo
          call LDT_releaseUnitNumber(ftn)
          deallocate(rlat)
          deallocate(rlon)
       endif   ! End Apply Mask condition
#endif

    enddo   ! End nest loop

!    print*, 'writing '
!    open(100,file='conusmask.bin',form='unformatted')
!    write(100) nldas1_struc(1)%conusmask
!    close(100)
!    stop

  end subroutine init_NLDAS1

end module nldas1_forcingMod
