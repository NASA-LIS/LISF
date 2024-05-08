!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE:  AGRMET_makest
! \label{AGRMET_makest}
!
! !REVISION HISTORY:
!
!     1 jul 89  initial version...................mr moore/sddc(agromet)
!    11 feb 91  added 'land' and 'estwt' to arg list of getcld. added   
!               'prob' to AGRMET_makest and calest arg list. eliminated call   
!               to setcoe by adding call to bndslc and adding some  
!               former setcoe code up in makest btwn calls to bndslc
!               and calest. reorganized prolog method discussion. re-   
!               sized pcoef and pcout arrays to reflect new latbnd  
!               and cldtyp array size.  isolated debug stmt. declared   
!               misc loop variables...............mr moore/sddc(agromet)
!    04 may 91  chgd init of 'ra' and 'e' arrays to 9999.0 so rainrates 
!               of 0.0 can be used. chg logic in setting rpcoef. added  
!               'thres' to arg lists so as to pass it down thru getcld  
!               to lodcld to eliminate multiple calls to func 'int'.
!               updtd prolog......................mr moore/sddc(agromet)
!     3 mar 92  made chgs necessary to include cldamt as an index of
!               pcoef and pcout arrays.  fixed error in process to pick 
!               rpcoef. movd read/write of pcoef up to driver and de-   
!               leted pcout array. increased the # of cld types used
!               ...................................mr moore/sys(agromet)
!     1 jun 92  created cac to make cld amt index of pcoef adjustable   
!               fm control file.  added use of rh class and stability   
!               class as indices to pcoef. expanded pcoef array and 
!               reduced prob array.  added use of pwgt..................
!               ...................................mr moore/sys(agromet)
!     1 feb 94  deleted stability class, rh class, cldtyp and latbnd
!               as pcoef indices. added cldtop ht class and tropopause  
!               ht class as pcoef indices.  added calls to geoslc (in   
!               place of bndslc) and zrobog (to zero-fill around bogus  
!               values of zero.  omitted write stmts that produced  
!               smiedr status print.  expanded cldamt array to account  
!               for 3-hrly times. updtd prolog.......................   
!               ....................mr moore, capt bertone, sys(agromet)
!    25 apr 94  added 'files accessed' section as 'id'ed during the 
!               94-01 qaa.................... capt bertone sysm(agromet)
!    15 jun 96  eliminated loading of previous precip estimates in  
!               back end of 'e' array.  reduced array 'e' size from 
!               (64,64,12) to (64,64,4).  eliminated writing of 3-hrly  
!               estd precip amts to file, a function now handled in 
!               phsrel routine.  moved calc of pprob from calest to 
!               to this routine. updtd prolog and brought up to stds.   
!               ..................................mr moore sysm(agromet)
!    07 mar 97  added subroutine crest to create an estimate based on   
!               climatological precip histories, which replaced the 
!               precip coefficient and tropopause height class estimate 
!               routines.  altered the arrangement of the code to run   
!               more efficiently.  brought code up to current afgwc and 
!               sysm standards..........capt andrus, ssgt mccormick/sysm
!    24 mar 97  added the 'source' array which is argued through this   
!               routine into the 'calest' routine.......ssgt miller/sysm
!    19 may 97  redimensioned source array to write each 3-hour sourced 
!               estimate into its own file..............ssgt miller/sysm
!    04 dec 97  changed crest routine to only produce a rtneph estimate 
!               and renamed rtnest. the rtnest routine also now outputs 
!               rtneph precip estimates. modified makest routine to call
!               cliest for a pure climo based estimate and smiest for 
!               an ssmi based estimate. added getcli routine to retrieve
!               monthly climo amounts and interpolate monthly values to 
!               the current day. modified calest to identfy all source
!               type based upon a pure hierarchy. file assign also now
!               assigns the additional and rtneph data files. updated
!               prolog and brought up to standards......................
!               ...............................capt andrus/dnxm(agromet)
!     2 sep 99  added option for geoswch to be 0, 1, or 2.  this allows
!               the rank values computed by geo_precip to either be 
!               used (geoswch=2) or not used (geoswch=1) in the final
!               calculation of the precip estimate.....capt hidalgo/dnxm
!     07 oct 99 ported to ibm sp-2, updated prolog, incorporated
!               FORTRAN 90 features.................capt hidalgo/agrmet
!     09 may 00 modified list of arguments in call to smiest............
!               .....................................capt hidalgo/agrmet
!     01 jun 00 removed estwt from program.  estwt was used to control
!               the way different estimates are ordered for use in the
!               calest sub-routine. changed grnk to integer, added grnk
!               to calest call......................ssgt campbell/agrmet
!     21 feb 01 modified loop and changed some variable names
!               to reflect change to 6-hrly cycles.........mr gayno/dnxm
!     10 jun 02 modified for incorporation of cdfs2 data to replace
!               rtneph.....................................mr gayno/dnxm
! 29Jul2005 Sujay Kumar, Initial Code
!     31 MAR 2010 add multiple resolution handling for ssmis and geoprecip.
!                 ........see in code comments for where..Michael Shaw/WXE
!     10 May 2013 added cmorph and gfs
!                 ................................Ryan Ruhge/16WS/WXE/SEMS
!     11 Aug 2016 Add support to process GALWEM precipitation
!                 ................................James Geiger/NASA/GSFC
!     23 Jan 2017 Return first guess field........Eric Kemp/NASA/GSFC
!     29 May 2017 only convert precipitation from mm/hr to mm/3hr
!                 if that data is used and present.......Chris Franks/SEMS
!
! !INTERFACE:    
subroutine AGRMET_makest(n,findex,j6hr,estpcp,source,cdfs2est,prcpwe, use_twelve, cmorphpixel,back4)
! !USES: 
  use LIS_coreMod,               only : LIS_rc
  use agrmet_forcingMod,         only : agrmet_struc
  use LIS_logMod,                only : LIS_abort, LIS_endrun, &
           LIS_getNextUnitNumber, LIS_releaseUnitNumber 

  implicit none

  integer, intent(in)  :: n
  integer, intent(in)  :: findex
  integer, intent(in)  :: j6hr
  real, intent(out)    :: estpcp(LIS_rc%lnc(n), LIS_rc%lnr(n),4)
  integer, intent(out) :: source(LIS_rc%lnc(n), LIS_rc%lnr(n),4)
  real, intent(out)    :: cdfs2est(LIS_rc%lnc(n), LIS_rc%lnr(n),4)
  real                 :: prcpwe(LIS_rc%lnc(n), LIS_rc%lnr(n),4)
!  real, intent(inout) :: back4(LIS_rc%lnc(n), LIS_rc%lnr(n),4) ! EMK
  real, intent(inout) :: back4(LIS_rc%gnc(n), LIS_rc%gnr(n),4) ! EMK

  logical, intent(in)  :: use_twelve
  real :: gridDesci(50)
  real :: gridDesco(50)
  real :: xmeshl,xmeshl2
  real :: xpnmcaf,xpnmcaf2
  real :: ypnmcaf,ypnmcaf2
  real :: xmesh, xmesh2,orient,orient2,xi1,xi12,xj1,xj12
  real :: alat1,alat12,alon1,alon12
  integer :: ihemi
  integer :: kprs
  logical, intent(out) :: cmorphpixel(LIS_rc%lnc(n), LIS_rc%lnr(n),4)
! declarations for readmask in geoprecip latlon
      character*9                   :: cstat
      character*100                 :: file_name,file_nam
      character*255                 :: message(20)
      integer                       :: rec_length
      integer                       :: istat
      integer                       :: istat1
      integer                       :: ftn
!
! !DESCRIPTION:
!
!    to make a consolidated 'final' precipitation amount estimate based 
!    upon several possible precipitation estimates.   
!
!    \textbf{Method} \newline
!    
!    - initialize output arrays.  \newline
!    - retrieve and interpolate any needed climatological data. \newline
!    - using cliest, make pure climo precip estimate if desired. \newline
!    - loop back through the 6-hour cycle at 3-hour increments \newline
!      - initialize local precip estimate arrays.   \newline
!      - if specified use ssmi data then generate the ssmi estimates
!        using smiest. \newline
!      - if specified use climo precip data and cdfs2 cloud data   
!        to generate a cdfs2 based precip estimate in routine cdfs2\_est.  \newline
!      - using specified hierarchy and data availability determine an   
!        estimate and identify the source of the estimate for every 
!        grid point using subroutine calest. \newline
!
!  The arguments and variables are: 
!  \begin{description}
!  \item[alert\_number]  alert\_number
!  \item[cdfs2est]      cdfs2 cloud based precip estimate (mm/3hr)
!  \item[cdfs2int]      cdfs2 time interval to look for cloud amount
!  \item[cdfs2swch]     cdfs2-based estimate use switch 
!  \item[cest]          pure climatological based precip estimate (mm/3hr)   
!  \item[cmorphpixel]   indicates if each pixel has cmorph data
!  \item[estpcp]        final estimated precip amounts (mm/3hr)  
!  \item[gest]          geo-precip estimated precip values
!  \item[grnk]	   geo-precip rank values for goodness of estimates
!                                 1 = better than ssmi est.
!                                 2 = better than present wx based est.
!                                 3 = better than cdfs2 est.
!                                 4 = default value of geo\_precip est.
!                                 5 = worse than climatological est.
!  \item[hemi]          hemisphere (1=north, 2=south)
!  \item[j6hr]          julian hour of beginning of 12-hour period   
!  \item[j3hr]          loop index,  3-hour julian hour  
!  \item[k]             loop counter (1 through 4)   
!  \item[prcpwe]        bogus precip rates  (mm/3hr)
!  \item[ra]            ssmi precip estimate
!  \item[cmorphdata]    cmorph precip estimate
!  \item[fg\_data]      first guess (GFS or GALWEM) precip estimate
!  \item[fc\_hr]        specify forecast hour when reading gfs files
!  \item[source]        array that contains source-of-the-estimate
!                  information with specified source flags
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[AGRMET\_read\_pcpclimodata] (\ref{AGRMET_read_pcpclimodata}) \newline
!   read precip climatology files
!  \item[AGRMET\_cliest] (\ref{AGRMET_cliest}) \newline
!    create purely climatologically based precip estimate
!  \item[AGRMET\_getcmorph] (\ref{AGRMET_getcmorph}) \newline
!    CMORPH-based estimate
!  \item[AGRMET\_smiest] (\ref{AGRMET_smiest}) \newline
!    create an SSM/I based estimate
!  \item[AGRMET\_cdsf2\_est] (\ref{AGRMET_cdfs2_est}) \newline
!    CDFS2-based estimate
!  \item[AGRMET\_geoest] (\ref{AGRMET_geoest}) \newline
!    GEOPRECIP-based estimtate
!  \item[AGRMET\_calest] (\ref{AGRMET_calest}) \newline
!    calculate a final estimate
!  \end{description}
!EOP
  integer           :: t
  integer           :: mesh 
  integer           :: j3hr
  real              :: cest(LIS_rc%lnc(n), LIS_rc%lnr(n))
  integer           :: i,j,k,count
  real              :: gest(LIS_rc%lnc(n), LIS_rc%lnr(n))
  integer           :: grnk(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real              :: ra(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real              :: cmorphdata(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real              :: fg_data(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real, allocatable :: fg_data_glb(:,:) ! EMK
  integer           :: fc_hr
  integer           :: alert_number
  real, parameter   :: quad9r = 9999.0   
  integer           :: days(12)
! Michael Shaw - Don't need these buffers, but hold over from previous versions
  byte, allocatable :: buffer(:,:,:)
  integer,allocatable  :: buffer2(:,:,:)
  data days    / 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /  
  count = 0    
  
!     ------------------------------------------------------------------
!     executable code begins here ...
!     initialize various arrays, counter and accumulator variables, 
!     statistical variables, and logical variables
!     ------------------------------------------------------------------
  cest     =   LIS_rc%udef
  ra       =   LIS_rc%udef
  cmorphdata = LIS_rc%udef

  allocate(fg_data_glb(LIS_rc%gnc(n), LIS_rc%gnr(n))) ! EMK

!---- for (0, 6Z], k=1, 2. for (6, 12Z], k= 3, 4
  if (use_twelve)  then
    k = 2 
    do t=3,4
       estpcp(:,:,t)   = LIS_rc%udef
       source(:,:,t)   = -9999
       cdfs2est(:,:,t) = LIS_rc%udef
       back4(:,:,t) = LIS_rc%udef ! EMK
    enddo
  else 
    k = 0
    do t=1,2
       estpcp(:,:,t)   = LIS_rc%udef
       source(:,:,t)   = -9999
       cdfs2est(:,:,t) = LIS_rc%udef
       back4(:,:,t) = LIS_rc%udef ! EMK
    enddo
  end if

  if((agrmet_struc(n)%cdfs2swch.eq.1).or.&
       (agrmet_struc(n)%clswch.eq.1) ) then 
     call AGRMET_read_pcpclimodata(n)
!     ------------------------------------------------------------------
!     if cliprc has good values, convert the units to mm/3hr
!     ------------------------------------------------------------------
     do j=1,LIS_rc%lnr(n)
        do i=1,LIS_rc%lnc(n)
           if(agrmet_struc(n)%cliprc(i,j).lt.9990.0) then 
!     ------------------------------------------------------------------
!         cliprc has units of mm/month.  divide by the number of days in
!         the month to get mm/day
!         ...then divide by 8 to get mm/3hr
!     ------------------------------------------------------------------
              agrmet_struc(n)%cliprc(i,j) = &
                   agrmet_struc(n)%cliprc(i,j)&
                   /float(days(LIS_rc%mo))
              agrmet_struc(n)%cliprc(i,j) = &
                   agrmet_struc(n)%cliprc(i,j)/8.0
           endif
        enddo
     enddo
  endif
!     ------------------------------------------------------------------
!       if clippd has good values, convert the units to mm/3hr
!     ------------------------------------------------------------------
                 
  if(agrmet_struc(n)%cdfs2swch.eq.1) then
     do j=1,LIS_rc%lnr(n)
        do i=1,LIS_rc%lnc(n)
           if (agrmet_struc(n)%clippd(i,j) .gt. -9990.0) then
              
!     ------------------------------------------------------------------
!         clippd has units of mm/day.  divide by 8 to get mm/3hr
!     ------------------------------------------------------------------

              agrmet_struc(n)%clippd(i,j) = agrmet_struc(n)%clippd(i,j)/8.0
           endif
        enddo
     enddo
  endif

  if( agrmet_struc(n)%clswch .eq. 1 )then
     call AGRMET_cliest( n, cest, agrmet_struc(n)%cliprc, &
          agrmet_struc(n)%clmult,LIS_rc%udef)
  endif

!     ------------------------------------------------------------------
!     loop backwards thru 6-hr period at 3-hr increments.
!     all satellite based estimates will produce a different 
!     precip est for each 3-hour time.
!     ------------------------------------------------------------------

  do j3hr = j6hr+3, j6hr+6, 3
     count =count+1
     k = k + 1   

!     ------------------------------------------------------------------
!       create a final precip estimate based on specified hierarchy.
!     ------------------------------------------------------------------
!       if the CMORPH switch is set, get CMORPH data
!     ------------------------------------------------------------------
     if(agrmet_struc(n)%cmorswch .eq. 1) then !if using CMORPH estimate
       call AGRMET_getcmorph(n,cmorphdata,j3hr,LIS_rc%udef)
     endif

!     ------------------------------------------------------------------
!       if the ssmi switch is set and CMORPH unavailable, determine the ssmi estimate
!     ------------------------------------------------------------------

! Michael Shaw - handle various polar stereo and latlon ssmis - see commetns
! for the similar geoprecip handling for some more details (below)
     if(agrmet_struc(n)%raswch .eq.1 )then !if using ssmis edr data
        if(agrmet_struc(n)%imaxsmi == 1440)  agrmet_struc(n)%global_or_hemi=1
        if(agrmet_struc(n)%imaxsmi /= agrmet_struc(n)%imax) agrmet_struc(n)%diff_grid=1
        if(agrmet_struc(n)%diff_grid == 1)then
           call AGRMET_smiest( n, j3hr, LIS_rc%udef, ra, &
                agrmet_struc(n)%razero, &
                alert_number,agrmet_struc(n)%imax3,agrmet_Struc(n)%jmax3)
        else
           call AGRMET_smiest( n, j3hr, LIS_rc%udef, ra, &
                agrmet_struc(n)%razero, &
                alert_number,agrmet_struc(n)%imax,agrmet_struc(n)%jmax)
        endif
        if(agrmet_struc(n)%imaxsmi == 1440)  agrmet_struc(n)%global_or_hemi=0
        if(agrmet_struc(n)%imaxsmi /= agrmet_struc(n)%imax) agrmet_struc(n)%diff_grid=0
     endif ! switch for smi

     fg_data = LIS_rc%udef
     fg_data_glb = LIS_rc%udef ! EMK

!     ------------------------------------------------------------------
!       if the GFS precip switch is set, get GFS precip
!     ------------------------------------------------------------------
     if(agrmet_struc(n)%gfsprecswch .eq. 1) then !if using GFS precip
       if ((k .eq. 1) .or. (k .eq. 3)) then
         fc_hr=3
         call AGRMET_fldbld_precip_gfs(n,findex,j6hr,fc_hr,fg_data)
       else if ((k .eq. 2) .or. (k .eq. 4)) then
         fc_hr=6
         call AGRMET_fldbld_precip_gfs(n,findex,j6hr,fc_hr,fg_data)
       end if
     endif

!     ------------------------------------------------------------------
!       if the GALWEM precip switch is set, get GALWEM precip
!     ------------------------------------------------------------------
     if(agrmet_struc(n)%galwemprecswch .eq. 1) then !if using GALWEM precip
       if ((k .eq. 1) .or. (k .eq. 3)) then
         fc_hr=3
         call AGRMET_fldbld_precip_galwem(n,j6hr,fc_hr,fg_data)
       else if ((k .eq. 2) .or. (k .eq. 4)) then
         fc_hr=6
         call AGRMET_fldbld_precip_galwem(n,j6hr,fc_hr,fg_data)
       end if
     endif

!     ------------------------------------------------------------------
!       if the cdfs2 switch is set, retrieve cdfs2 cloud amounts, and
!       create an cdfs2 based precip estimate.
!       Ryan Ruhge - held here for legacy options but will
!        not be included when AGRMET_calest makes the precip blend.
!     ------------------------------------------------------------------

     if( agrmet_struc(n)%cdfs2swch .eq. 1 ) then
        call AGRMET_cdfs2_est( n, k, &
             agrmet_struc(n)%cliprc, agrmet_struc(n)%clippd,&
             agrmet_struc(n)%clirtn, agrmet_struc(n)%mnpr, &
             agrmet_struc(n)%mxpr, &
             agrmet_struc(n)%mnpd, agrmet_struc(n)%mxpd, &
             agrmet_struc(n)%cldth, agrmet_struc(n)%mdmv, &
             agrmet_struc(n)%mdpe, &
             agrmet_struc(n)%ovpe, agrmet_struc(n)%cdfs2int, cdfs2est,&
             alert_number, j3hr )
     endif

!     ------------------------------------------------------------------
!       if geo switch is set and cmorph is unavailable, retrieve geo_precip-based precip estimate.
!       geoswch = 0 indicates do not use the geo-precip estimate.
!       geoswch = 1 indicates use the estimate but do not use the rank
!                   values.
!       geoswch = 2 indicates use the estimate and use the rank values
!     ------------------------------------------------------------------

! Michael Shaw - handle various polar stereo and latlon geoprecip - similar to 
! ssmis handling above
     gest = LIS_rc%udef
     grnk = 4 
     if( ( agrmet_struc(n)%geoswch .eq. 1 ) .or. & !if using geoprecip
             (agrmet_struc(n)%geoswch .eq. 2) )then !if using rank files
       ! if the geoprecip is global or at least different from "native" agrmet 8th mesh
       ! grid, need to indicate this, because interp_agrmetvar and AGRMET_geoest
       ! handle these differently from "same basis as input grids"
       if(agrmet_struc(n)%imaxgp == 3600)                 agrmet_struc(n)%global_or_hemi=1
       if(agrmet_struc(n)%imaxgp /= agrmet_struc(n)%imax) agrmet_struc(n)%diff_grid=1
       if(agrmet_struc(n)%diff_grid == 1)then
          if(agrmet_struc(n)%imaxgp == 4096)then
             if(agrmet_struc(n)%imax == 512)  mesh=8 !Factor difference b/w meshes since using
             if(agrmet_struc(n)%imax == 1024) mesh=4 !16th mesh land mask for 64th mesh for the moment
             call AGRMET_geoest( n, j3hr, agrmet_struc(n)%land, gest, &
                  grnk, LIS_rc%udef, agrmet_struc(n)%imax2, agrmet_struc(n)%jmax2,&
                  agrmet_struc(n)%geoswch,alert_number,mesh)
          else
             call AGRMET_geoest( n, j3hr, agrmet_struc(n)%land2, gest, &
                  grnk, LIS_rc%udef, agrmet_struc(n)%imax2, agrmet_struc(n)%jmax2,&
                  agrmet_struc(n)%geoswch,alert_number,1)
          endif
       else
          call AGRMET_geoest( n, j3hr, agrmet_struc(n)%land, gest, &
               grnk, LIS_rc%udef, agrmet_struc(n)%imax, agrmet_struc(n)%jmax,&
               agrmet_struc(n)%geoswch,alert_number,1) !,shemi,nhemi)
       endif
       ! if the geoprecip is global or at least different from "native" agrmet 8th mesh
       ! grid, need to reset the indicators as other grids (precip, clouds, radiation)
       ! use interp_agrmetvar, too, and so need to be interpolated correctly while
       ! that routine handles things differently depending on global or otherwise
       ! different grid 
       if(agrmet_struc(n)%imaxgp == 3600)                 agrmet_struc(n)%global_or_hemi=0
       if(agrmet_struc(n)%imaxgp /= agrmet_struc(n)%imax) agrmet_struc(n)%diff_grid=0
       endif !switch for gp 

!     ------------------------------------------------------------------
!       create a final precip estimate based on specified hierarchy.
!       for each grid point identify the source of the grid point   
!       estimate, and generate statistics.  
!     ------------------------------------------------------------------
       back4(:,:,k) = fg_data_glb(:,:) ! EMK

        call AGRMET_calest( n, cest, fg_data, estpcp, gest, k,&
             prcpwe, ra, cmorphdata, cmorphpixel, source, &
             LIS_rc%lnc(n), LIS_rc%lnr(n), grnk )

     enddo

     deallocate(fg_data_glb) ! EMK

   end subroutine AGRMET_makest
