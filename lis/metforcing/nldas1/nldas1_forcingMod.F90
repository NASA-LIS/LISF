!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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
!  The implemenatation in LIS has the derived data type {\tt nldas1\_struc}
!  that includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
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
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Array containing the neighbor information of the input grid 
!    for each grid point in LIS, for nearest neighbor interpolation.
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
!                               field choices to the lis.config file
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
     real         :: ts
     integer      :: ncold, nrold    ! AWIPS 212 dimensions
     character*50 :: nldas1_filesrc
     character*80 :: nldas1dir       ! NLDAS-1 Forcing Directory
     character*80 :: ediff_file      ! Original Elevdiff File
     character*80 :: elevfile
     character*40 :: prec_field,swdn_field
     real*8       :: nldas1time1,nldas1time2

     integer :: mi
     integer            :: applymask
     real, allocatable      :: conusmask(:,:)
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
     integer            :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
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
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun
    use map_utils,      only : ij_to_latlon, proj_latlon

    implicit none
    
    integer, intent(in) :: findex
    
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for NLDAS-1
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
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
!   \item[read\_edas\_elev](\ref{read_edas_elev}) \newline
!    reads the native elevation of the EDAS data to be used
!    for topographic adjustments to the metforcing data
!  \end{description}
!EOP
    
    real :: gridDesci(LIS_rc%nnest, 50)
    real, allocatable :: rlat(:,:), rlon(:,:)
    integer :: n
    integer :: ftn
    integer :: nc_dom
    integer :: line1, line2, line
    integer :: c,r

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the NLDAS-1 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
    
    allocate(nldas1_struc(LIS_rc%nnest))
    call readcrd_nldas1()

    do n=1, LIS_rc%nnest
       nldas1_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, nldas1_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9

    nldas1_struc(:)%ncold = 464
    nldas1_struc(:)%nrold = 224

    gridDesci = 0 

    do n=1,LIS_rc%nnest

       allocate(nldas1_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(nldas1_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       nldas1_struc(n)%metdata1 = 0
       nldas1_struc(n)%metdata2 = 0

       gridDesci(n,1) = 0
       gridDesci(n,2) = nldas1_struc(n)%ncold
       gridDesci(n,3) = nldas1_struc(n)%nrold
       gridDesci(n,4) = 25.0625
       gridDesci(n,5) = -124.9375
       gridDesci(n,6) = 128
       gridDesci(n,7) = 52.9375
       gridDesci(n,8) = -67.0625
       gridDesci(n,9) = 0.125
       gridDesci(n,10) = 0.125
       gridDesci(n,20) = 64

       if( gridDesci(n,9)  == LIS_rc%gridDesc(n,9) .and. &
           gridDesci(n,10) == LIS_rc%gridDesc(n,10).and. &
           LIS_rc%gridDesc(n,1) == proj_latlon .and. &
           LIS_rc%met_interp(findex) .ne. "neighbor" ) then
         write(LIS_logunit,*) "WARNING MSG:  The NLDAS1 grid was selected for the"
         write(LIS_logunit,*) "  LIS run domain; however, 'bilinear', 'budget-bilinear',"
         write(LIS_logunit,*) "  or some other unknown option was selected to spatially"
         write(LIS_logunit,*) "  downscale the grid, which will cause errors during runtime."
         write(LIS_logunit,*) "Program stopping ..."
         call LIS_endrun()
       endif

       nldas1_struc(n)%mi = nldas1_struc(n)%ncold*nldas1_struc(n)%nrold

!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(nldas1_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))


          call bilinear_interp_input(n,gridDesci(n,:),&
               nldas1_struc(n)%n111,nldas1_struc(n)%n121,nldas1_struc(n)%n211,&
               nldas1_struc(n)%n221,nldas1_struc(n)%w111,nldas1_struc(n)%w121,&
               nldas1_struc(n)%w211,nldas1_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(nldas1_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas1_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               nldas1_struc(n)%n111,nldas1_struc(n)%n121,&
               nldas1_struc(n)%n211,nldas1_struc(n)%n221,&
               nldas1_struc(n)%w111,nldas1_struc(n)%w121,&
               nldas1_struc(n)%w211,nldas1_struc(n)%w221)

          allocate(nldas1_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas1_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas1_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas1_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas1_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               nldas1_struc(n)%n112,nldas1_struc(n)%n122,&
               nldas1_struc(n)%n212,nldas1_struc(n)%n222,&
               nldas1_struc(n)%w112,nldas1_struc(n)%w122,&
               nldas1_struc(n)%w212,nldas1_struc(n)%w222)

       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 

          allocate(nldas1_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               nldas1_struc(n)%n113)
       else
          write(LIS_logunit,*) 'Interpolation option not specified for NLDAS1'
          write(LIS_logunit,*) 'Program stopping...'
          call LIS_endrun()
       endif

    !- If Elevation Correction desired, remove existing NLDAS-1 1/8 degree
    !   elevation correction that is currently built into forcing:

       if ( trim(LIS_rc%met_ecor(findex)).ne."none") then     ! ACTIVATE ELEV CORRECTION

       !- Read in original NLDAS-1 EDAS-GTOPO Elevation Difference file::
          allocate(nldas1_struc(n)%orig_ediff(&
                     nldas1_struc(n)%ncold*nldas1_struc(n)%nrold))
          call read_orig_nldas1_elevdiff(n)

       !- Read EDAS elevation file to correct NLDAS-1 EDAS fields:
          call read_edas_elev(n,findex)    ! NEW code -- read EDAS height file only.

       end if
       
       if(nldas1_struc(n)%applymask.eq.1) then 
          allocate(nldas1_struc(n)%conusmask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(rlat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(rlon(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          do r=1,LIS_rc%lnr(n)
             do c=1,LIS_rc%lnc(n)
                call ij_to_latlon(LIS_domain(n)%lisproj,float(c),float(r),&
                     rlat(c,r),rlon(c,r))                
             enddo
          enddo
          
          ftn = LIS_getNextUnitNumber()
          write(LIS_logunit,*) 'Reading the CONUS mask: ',&
               trim(nldas1_struc(n)%conusmaskfile)
          open(ftn,file=trim(nldas1_struc(n)%conusmaskfile),form='unformatted',&
               recl=4,access='direct',status='old')
          nc_dom = nint((nldas1_struc(n)%mask_gd(4)-nldas1_struc(n)%mask_gd(2))/&
               nldas1_struc(n)%mask_gd(5))+1
          do r=1,LIS_rc%lnr(n)
             do c=1,LIS_rc%lnc(n)
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
          call LIS_releaseUnitNumber(ftn)
          
          deallocate(rlat)
          deallocate(rlon)
       endif
    enddo
!    print*, 'writing '
!    open(100,file='conusmask.bin',form='unformatted')
!    write(100) nldas1_struc(1)%conusmask
!    close(100)
!    stop

  end subroutine init_NLDAS1
end module nldas1_forcingMod
