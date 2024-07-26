!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
!BOP
! 
! !MODULE: LVT_TSMod 
! \label(LVT_TSMod)
!
! !INTERFACE:
module LVT_TSMod
! 
! !USES: 
  use ESMF
  use LVT_logMod
  use LVT_coreMod
  use LVT_statsDataMod
  use LVT_histDataMod
  use LVT_CIMod
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module provides tools to generate time series output of 
!  model and observational data. The code supports the extraction of 
!  both point data and area averaged data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 2 Oct 2008: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! 
! 
! 
  
  private
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LVT_TSinit
  PUBLIC :: LVT_writeTSinfo
  PUBLIC :: LVT_writeDataBasedStrat
  PUBLIC :: LVT_writeSeasonalCycleInfo
  PUBLIC :: LVT_writeAvgDiurnalCycleInfo
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LVT_TSobj
!EOP

  type, public :: ts_struc

     character*500  :: tslocfile
     character*500  :: tslocname
     real         :: tslat1
     real         :: tslon1
     real         :: tslat2
     real         :: tslon2
     integer      :: npts
     integer      :: ts_cindex1
     integer      :: ts_rindex1
     integer      :: ts_cindex2
     integer      :: ts_rindex2
     real         :: ts_min_pts
     integer, allocatable  :: ts_tindex(:)

  end type ts_struc

  type(ts_struc), allocatable  :: LVT_TSobj(:)
  real, parameter          :: min_param = 1E8
  real, parameter          :: max_param = -1E8

  interface LVT_writeTSinfo
     module procedure writeTSinfo_noensem_noobs
     module procedure writeTSinfo_noensem
     module procedure writeTSinfo3
     module procedure writeTSinfo4
     module procedure writeTSinfo5
  end interface

  interface LVT_writeDataBasedStrat
     module procedure writeDataBasedStrat1
     module procedure writeDataBasedStrat2
  end interface

  interface LVT_writeSeasonalCycleInfo
     module procedure writeSeasonalCycleInfo1
     module procedure writeSeasonalCycleInfo2
  end interface

  interface LVT_writeAvgDiurnalCycleInfo
     module procedure writeAvgDiurnalCycleInfo1
     module procedure writeAvgDiurnalCycleInfo2
  end interface

contains

!BOP
! 
! !ROUTINE: LVT_TSinit
!  \label{LVT_TSinit}
!
! !INTERFACE:   
  subroutine LVT_TSinit()
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine performs the initialization steps required for
!   time series extraction.  The subroutine reads the time series
!   location file, which specifies the locations of interest.
! 
!  The locations can be specified in five different formats:
!  (1) using the lat/lon values; (2) using the column/row indices;
!  (3) using the tile indices; (4) specifying lat/lon values to
!  draw a polygon around a region; and (5) using a categorical
!  from which to define subregions.
! 
!  Note that LVT has been updated so the format of the time series
!  locations file uses a minimum fraction of the domain before the
!  temporal calculations will occur.  Previously, the time series
!  locations file used a minimum number of observations.  A value
!  for this ``min frac'' of 0.1 (for example) implies that at least
!  10 percent of the total number of points in the domain location
!  must be available for the temporal calculations to occur.
! 
!  A sample file with location style 1 is shown below:
! 
!  \begin{verbatim}
!   #Number of locations
!   2
!   #Location style (1-lat/lon, 2-col/row, 3-tile, 4-polygon, 5-map)
!   1
!   #Name (then, next line), SW-lat, SW-lon, NE-lat, NE-lon, min frac
!   WEST_US
!   40.0 -130.0 50.0 -110.0 0.0
!   HIGH_PLAINS_US
!   43.0 -110.0 49.0 -100.0 0.0
!  \end{verbatim}
! 
!  If the location style is 2, the user specifies the column and
!  row indices for the bounding boxes, instead of the corner lat/lon
!  values.  A sample file with location style 2 is shown below:
! 
!  \begin{verbatim}
!   #Number of locations
!   2
!   #Location style (1-lat/lon, 2-col/row, 3-tile, 4-polygon, 5-map)
!   2
!   #Name (then, next line), SW-col, SW-row, NE-col, NE-row, min frac
!   WEST_US
!    1 1 20 30 0.0
!   EAST_US
!   21 1 40 30 0.0
!  \end{verbatim}
! 
!  If the location style is 3, the user specifies the tile indices
!  for specifying the bounds (starting tile index and ending tile
!  index).  A sample file with location style 3 is shown below:
! 
!  \begin{verbatim}
!   #Number of locations
!   2
!   #Location style (1-lat/lon, 2-col/row, 3-tile, 4-polygon, 5-map)
!   3
!   #Name (then, next line), Start index, End index, min frac
!   WEST_US
!    1 20 0.0
!   EAST_US
!   21 40 0.0
!  \end{verbatim}
! 
!  If the location style is 4, the user explicitly specifies the
!  lat/lons of each grid point to be used to specify a region in
!  the shape of a polygon.  Users should be careful with this
!  location style option, as they cannot specify a minimum fraction
!  of the domain that must have valid observations.  A sample file
!  with location style 4 is shown below:
!  
!  \begin{verbatim}
!   #Number of locations
!   2
!   #Location style (1-lat/lon, 2-col/row, 3-tile, 4-polygon, 5-map)
!   4
!   #Number of points followed by lat/lon of each point
!   REGION1
!   3
!   34.4 -103.2
!   33.4 -100.2
!   32.1  -99.3
!   REGION2
!   2
!   40.2 -103.3
!   42.2 -104.2
!  \end{verbatim}
! 
!  If the location style is 5, the user explicitly specifies a
!  categorical map from which to define subregions.  In the map,
!  the categories must be in numerically increasing order from 1.
!  The map must be a binary direct-access file, with point (1,1)
!  in the southwest corner of the domain.  A sample file with
!  location style 5 is shown below:
! 
!  \begin{verbatim}
!   #Number of stations
!   3
!   #Location style (1-lat/lon, 2-col/row, 3-tile, 4-polygon, 5-map)
!   5
!   #Name (then, next line), min frac
!   NEWENGLAND
!   0.0
!   MIDATLANTIC
!   0.0
!   SOUTHATLANTIC
!   0.0
!   #categorical map
!   ../huc02_conus_0.125dg.1gd4r
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer             :: ftn
    integer             :: ts_tindex1(LVT_rc%ntslocs), ts_tindex2(LVT_rc%ntslocs)
    integer             :: i, c,r,k,kk 
    integer, allocatable    :: kk_vals(:)
    real                :: region_cat(LVT_rc%lnc,LVT_rc%lnr)
    character*500       :: region_map
    logical             :: file_exists
    real                :: col, row,lat,lon
    integer             :: rc

!    if(LVT_rc%extractTS.gt.0) then 
       call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%tsspecfile,&
            label="Time series location file:",rc=rc)
       
!       call LVT_verify(rc,'Time series location file: not defined')  
       if(rc.eq.0) then 
       open(10,file=(LVT_rc%tsspecfile),status='old')
       read(10,*) 
       read(10,*) LVT_rc%ntslocs
       read(10,*) 
       read(10,*) LVT_rc%tsspecstyle
       read(10,*) 
       
       allocate(LVT_TSobj(LVT_rc%ntslocs))
       
       write(LVT_logunit,*) '[INFO] Reading location data ....'
       do k=1,LVT_rc%ntslocs
          read(10,fmt='(a40)') LVT_TSobj(k)%tslocname
          write(LVT_logunit,*) '[INFO] Station name: ',trim(LVT_TSobj(k)%tslocname)
          if(LVT_rc%tsspecstyle.eq.1) then 
             read(10,*)  LVT_TSobj(k)%tslat1,LVT_TSobj(k)%tslon1,&
                  LVT_TSobj(k)%tslat2,LVT_TSobj(k)%tslon2, &
                  LVT_TSobj(k)%ts_min_pts
             write(LVT_logunit,*) '[INFO] Station location: ',&
                  LVT_TSobj(k)%tslat1,LVT_TSobj(k)%tslon1,&
                  LVT_TSobj(k)%tslat2,LVT_TSobj(k)%tslon2, &
                  LVT_TSobj(k)%ts_min_pts
          elseif(LVT_rc%tsspecstyle.eq.2) then 
             read(10,*)  LVT_TSobj(k)%ts_cindex1, LVT_TSobj(k)%ts_rindex1, &
                  LVT_TSobj(k)%ts_cindex2, LVT_TSobj(k)%ts_rindex2, &
                  LVT_TSobj(k)%ts_min_pts
             write(LVT_logunit,*) '[INFO] Station location (c,r): ',&
                  LVT_TSobj(k)%ts_cindex1, LVT_TSobj(k)%ts_rindex1, &
                  LVT_TSobj(k)%ts_cindex2, LVT_TSobj(k)%ts_rindex2, &
                  LVT_TSobj(k)%ts_min_pts
          elseif(LVT_rc%tsspecstyle.eq.3) then 
             read(10,*) ts_tindex1(k), ts_tindex2(k), LVT_TSobj(k)%ts_min_pts
             write(LVT_logunit,*) '[INFO] Station location (tile): ',&
                  ts_tindex1(k), ts_tindex2(k),LVT_TSobj(k)%ts_min_pts
          elseif(LVT_rc%tsspecstyle.eq.4) then 
             read(10,*)  LVT_TSobj(k)%npts
             allocate(LVT_TSobj(k)%ts_tindex(LVT_TSObj(k)%npts))
             
             do kk=1, LVT_TSobj(k)%npts
                read(10,*) lat,lon
                
                call latlon_to_ij(LVT_domain%lvtproj, &
                     lat, lon, col,row)
                LVT_TSobj(k)%ts_tindex(kk) = LVT_domain%gindex(&
                     nint(col),nint(row))
             enddo
! Assume with location style 4 that there must be at least a small
! fraction of points in the polygon with valid values.  Previously,
! this was hard-coded to at least one valid value in the polygon.
             LVT_TSobj(k)%ts_min_pts = 0.001
          elseif(LVT_rc%tsspecstyle.eq.5) then
             read(10,*) LVT_TSobj(k)%ts_min_pts
          endif
       enddo
       
       if(LVT_rc%tsspecstyle.eq.1) then 
          do i=1,LVT_rc%ntslocs
             call latlon_to_ij(LVT_domain%lvtproj, &
                  LVT_TSobj(i)%tslat1, LVT_TSobj(i)%tslon1, col,row)

             LVT_TSobj(i)%ts_cindex1 = nint(col)
             LVT_TSobj(i)%ts_rindex1 = nint(row)
             call latlon_to_ij(LVT_domain%lvtproj, &
                  LVT_TSobj(i)%tslat2, LVT_TSobj(i)%tslon2, col,row)

             LVT_TSobj(i)%ts_cindex2 = nint(col)
             LVT_TSobj(i)%ts_rindex2 = nint(row)
             
             !check if the domain is valid compared to the LVT domain. 
             if(LVT_TSobj(i)%ts_cindex1.ge.1.and. &
                  LVT_TSobj(i)%ts_cindex1.le.LVT_rc%lnc.and. &
                  LVT_TSobj(i)%ts_cindex2.ge.1.and. &
                  LVT_TSobj(i)%ts_cindex2.le.LVT_rc%lnc.and. &
                  LVT_TSobj(i)%ts_rindex1.ge.1.and. &
                  LVT_TSobj(i)%ts_rindex1.le.LVT_rc%lnr.and. &
                  LVT_TSobj(i)%ts_rindex2.ge.1.and. &
                  LVT_TSobj(i)%ts_rindex2.le.LVT_rc%lnr) then                 
                
                LVT_TSobj(i)%npts = ((LVT_TSobj(i)%ts_cindex2-&
                     LVT_TSobj(i)%ts_cindex1+1)*&
                     (LVT_TSobj(i)%ts_rindex2-LVT_TSobj(i)%ts_rindex1+1))
                
                allocate(LVT_TSobj(i)%ts_tindex(LVT_TSObj(i)%npts))
                
                kk = 0

                do r=LVT_TSobj(i)%ts_rindex1,LVT_TSobj(i)%ts_rindex2
                   do c=LVT_TSobj(i)%ts_cindex1,LVT_TSobj(i)%ts_cindex2
                      kk = kk+1    
                      if(c.ge.1.and.r.ge.1) then 
                         
                         LVT_TSobj(i)%ts_tindex(kk) = LVT_domain%gindex(c,r)
                      else
                         LVT_TSobj(i)%ts_tindex(kk) = -1
                      endif
                   enddo
                enddo
             else
                write(LVT_logunit,*) '[ERR] The subdomain '//trim(LVT_TSobj(i)%tslocname)
                write(LVT_logunit,*) '[ERR] is outside the LVT domain. Please modify '
                write(LVT_logunit,*) '[ERR] or remove it from the time series locations file'
                call LVT_endrun()
             endif
          enddo
       elseif(LVT_rc%tsspecstyle.eq.2) then 
          do i=1,LVT_rc%ntslocs             
             LVT_TSobj(i)%npts = &
                  ((LVT_TSobj(i)%ts_cindex2-LVT_TSobj(i)%ts_cindex1+1)*&
                  (LVT_TSobj(i)%ts_rindex2-LVT_TSobj(i)%ts_rindex1+1))
             
             allocate(LVT_TSobj(i)%ts_tindex(LVT_TSObj(i)%npts))
             kk = 0
             do r=LVT_TSobj(i)%ts_rindex1,LVT_TSobj(i)%ts_rindex2
                do c=LVT_TSobj(i)%ts_cindex1,LVT_TSobj(i)%ts_cindex2
                   kk = kk+1                                     
                   LVT_TSobj(i)%ts_tindex(kk) = LVT_domain%gindex(c,r)
                enddo
             enddo
          enddo
       elseif(LVT_rc%tsspecstyle.eq.3) then 
          do i=1,LVT_rc%ntslocs
             LVT_TSobj(i)%npts = (ts_tindex2(i)-ts_tindex1(i)+1)
             
             do kk=1,LVT_TSobj(i)%npts
                LVT_TSobj(i)%ts_tindex(kk) = ts_tindex1(i)+kk
             enddo
          enddo
       elseif(LVT_rc%tsspecstyle.eq.4) then !a series of points
          !processed above
          
       elseif(LVT_rc%tsspecstyle.eq.5) then 
!read the categorical data, considering each category as a separate 
          read(10,*)
          read(10,fmt='(a100)') region_map
          ftn = LVT_getNextUnitNumber()
          write(LVT_logunit,*) '[INFO] Reading categorical map: ',trim(region_map)
          inquire(file=trim(region_map),exist=file_exists) 
          if(file_exists) then 
             open(ftn,file=region_map,form='unformatted',access='direct',&
                  recl=LVT_rc%lnc*LVT_rc%lnr*4)
             read(ftn,rec=1) region_cat
             call LVT_releaseUnitNumber(ftn)
          else
             write(LVT_logunit,*) '[INFO] Categorical map '//trim(region_map)//' not found'
             call LVT_endrun()
          endif
          
          do k=1,LVT_rc%ntslocs
             LVT_TSobj(k)%npts = 0 
          enddo
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(region_cat(c,r).ne.LVT_rc%udef) then 
                   i = nint(region_cat(c,r))
                   ! The first n categories will be chosen
                   if(i.le.LVT_rc%ntslocs.and.i.gt.0) then 
                      LVT_TSobj(i)%npts = LVT_TSobj(i)%npts + 1
                   endif
                endif
             enddo
          enddo
          
          do k=1,LVT_rc%ntslocs
             allocate(LVT_TSobj(k)%ts_tindex(LVT_TSobj(k)%npts))
          enddo

          allocate(kk_vals(LVT_rc%ntslocs)) 
          kk_vals = 0
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(region_cat(c,r).ne.LVT_rc%udef) then 
                   i = nint(region_cat(c,r))
                   ! The first n categories will be chosen
                   if(i.gt.0.and.i.le.LVT_rc%ntslocs) then 
                      kk_vals(i) = kk_vals(i) +1
                      LVT_TSobj(i)%ts_tindex(kk_vals(i)) = &
                           LVT_domain%gindex(c,r)
                   endif
                endif
             enddo
          enddo
          deallocate(kk_vals)
       endif

    endif

    close(10)
  end subroutine LVT_TSinit

!BOP
! 
! !ROUTINE: writeTSinfo_noensem_noobs
! \label{writeTSinfo_noensem_noobs}
!
! !INTERFACE: 
  subroutine writeTSinfo_noensem_noobs(ftn_ts_loc,model,nsize,metric_ts,count_metric_ts)
! 
! !USES:     

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: ftn_ts_loc(LVT_rc%ntslocs)
    type(LVT_metaDataEntry) :: model
    integer                 :: nsize
    real                    :: metric_ts(nsize,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: count_metric_ts(nsize,model%selectNlevs,&
         LVT_rc%strat_nlevels)


    real     :: sum_v
    integer  :: nsum_v
    real     :: wsum_v
    integer  :: wnsum_v
    real     :: ens_std
    real     :: mu_m(1)
    integer  :: nmu_m(1)
    real     :: mu_2
    integer  :: nmu
    real     :: mean_v, med_v, sstd_v
    real     :: ci_val
    integer  :: i,k,l,kk,tid
    integer  :: stid, m,t
    real     :: metric_tsdom(nsize)
    real     :: maxv, minv
    integer  :: nensem

    nensem =1
    do i=1, LVT_rc%ntslocs
       do k=1,model%selectNlevs
          do l=1,LVT_rc%strat_nlevels
             sum_v = 0 
             nsum_v = 0 
             wsum_v = 0 
             wnsum_v = 0 
             maxv = max_param
             minv = min_param
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   t = tid
                   if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                      sum_v = sum_v + metric_ts(t,k,l)
                      nsum_v = nsum_v + 1
                      metric_tsdom(nsum_v) = metric_ts(t,k,l)

                      if(metric_ts(t,k,l).lt.minv) then 
                         minv = metric_ts(t,k,l)
                      endif
                      if(metric_ts(t,k,l).gt.maxv) then 
                         maxv = metric_ts(t,k,l)
                      endif
                   endif
                endif
             enddo
             if(maxv.eq.max_param) maxv = LVT_rc%udef 
             if(minv.eq.min_param) minv = LVT_rc%udef

!             if(nsum_v.ge.LVT_TSobj(i)%ts_min_pts.and.nsum_v.ne.0) then 
             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts).and.&
                  nsum_v.ne.0) then 
                mean_v = sum_v/nsum_v
             else
                mean_v = LVT_rc%udef
             endif

             sstd_v = 0 
             nsum_v = 0 

             if(mean_v.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1) then 
                      t = tid
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         sstd_v = sstd_v + &
                              (metric_ts(t,k,l)-mean_v)**2
                         nsum_v = nsum_v + 1
                      endif
                   endif
                enddo

                if(nsum_v.gt.0) then
                   sstd_v = sqrt(sstd_v/nsum_v)
                else
                   sstd_v = LVT_rc%udef
                endif
             else
                sstd_v = LVT_rc%udef
             endif
             if(nsum_v.le.1) then 
                sstd_v = LVT_rc%udef
             endif

             mu_m = 0 
             nmu_m = 0 
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)*nensem
                   do m=1,nensem
                      t = stid+m
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         mu_m(m) = mu_m(m) + metric_ts(t,k,l)
                         nmu_m(m) = nmu_m(m) + 1
                      endif
                   enddo
                end if
             enddo
             mu_2 = 0 
             nmu = 0 
             do m=1,nensem
                if(nmu_m(m).gt.0) then 
                   mu_m(m) = mu_m(m)/nmu_m(m)
                   mu_2 = mu_2 + mu_m(m)
                   nmu = nmu + 1
                endif
             enddo
             if(nmu.gt.0) then 
                mu_2 = mu_2/nmu
             else
                mu_2 = LVT_rc%udef
             endif
             
             ens_std = 0 
             if(mu_2.ne.LVT_rc%udef) then 
                do m=1,nensem
                   ens_std = ens_std + (mu_m(m) - mu_2)**2
                enddo
                ens_std = sqrt(ens_std/nensem)
             else
                ens_std = LVT_rc%udef
             endif             
             
             ci_val = LVT_rc%udef
             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                call LVT_computeCI(metric_tsdom(1:nsum_v),nsum_v,&
                     LVT_rc%pval_CI,ci_val)
             endif
             write(ftn_ts_loc(i),201,advance='no') &
                  mean_v, sstd_v, minv, maxv, ens_std, ci_val
          enddo
       enddo
    enddo

201 format(6E14.6)
  end subroutine writeTSinfo_noensem_noobs


!BOP
! 
! !ROUTINE: writeTSinfo_noensem
! \label{writeTSinfo_noensem}
!
! !INTERFACE: 
  subroutine writeTSinfo_noensem(ftn_ts_loc,model,nsize_m,&
       metric_ts,count_metric_ts,&
       nsize_o,obs_ts, count_obs_ts)
! 
! !USES:     

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: ftn_ts_loc(LVT_rc%ntslocs)
    type(LVT_metaDataEntry) :: model
    integer                 :: nsize_m
    real                    :: metric_ts(nsize_m,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: count_metric_ts(nsize_m,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: nsize_o
    real                    :: obs_ts(nsize_o,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: count_obs_ts(nsize_o,model%selectNlevs,&
         LVT_rc%strat_nlevels)


    real     :: sum_v,sum_v1
    integer  :: nsum_v,nsum_v1
    real     :: mean_v, sstd_v, mean_v1,sstd_v1
    real     :: mu_m(1)
    integer  :: nmu_m(1)
    real     :: mu_2
    integer  :: nmu
    real     :: ens_std
    real     :: ci_val,ci_val1
    integer  :: i,k,l,kk,tid
    integer  :: stid, m,t
    real     :: metric_tsdom(nsize_m)
    real     :: maxv, minv
    real     :: maxv1, minv1
    integer  :: nensem

    nensem = 1

    do i=1, LVT_rc%ntslocs
       do k=1,model%selectNlevs
          do l=1,LVT_rc%strat_nlevels
             sum_v = 0 
             nsum_v = 0 
             
             maxv = max_param
             minv = min_param

             do kk=1, LVT_TSobj(i)%npts 
                tid = LVT_Tsobj(i)%ts_tindex(kk)

                if(tid.ne.-1) then 
                   t = tid
                   if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                      sum_v = sum_v + metric_ts(t,k,l)
                      nsum_v = nsum_v + 1
                      metric_tsdom(nsum_v) = metric_ts(t,k,l)
                      if(metric_ts(t,k,l).lt.minv) then 
                         minv = metric_ts(t,k,l)
                      endif
                      if(metric_ts(t,k,l).gt.maxv) then 
                         maxv = metric_ts(t,k,l)
                      endif
                   endif
                endif
             enddo
             if(maxv.eq.max_param) maxv = LVT_rc%udef
             if(minv.eq.min_param) minv = LVT_rc%udef

             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                mean_v = sum_v/nsum_v
             else
                mean_v = LVT_rc%udef
             endif
             
             sstd_v = 0 
             nsum_v = 0 
             
             if(mean_v.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1) then 
                      t = tid
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         sstd_v = sstd_v + &
                              (metric_ts(t,k,l)-mean_v)**2
                         nsum_v = nsum_v + 1                            
                      endif
                   endif
                enddo
                if(nsum_v.gt.0) then 
                   sstd_v = sqrt(sstd_v/nsum_v)
                else
                   sstd_v = LVT_rc%udef
                endif
             else
                sstd_v = LVT_rc%udef
             endif
             if(nsum_v.le.1) sstd_v = LVT_rc%udef

             mu_m = 0 
             nmu_m = 0 
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   t = tid
                   m = 1
                   if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                      mu_m(m) = mu_m(m) + metric_ts(t,k,l)
                      nmu_m(m) = nmu_m(m) + 1
                   endif
                end if
             enddo
             mu_2 = 0 
             nmu = 0 
             do m=1,nensem
                if(nmu_m(m).gt.0) then 
                   mu_m(m) = mu_m(m)/nmu_m(m)
                   mu_2 = mu_2 + mu_m(m)
                   nmu = nmu + 1
                endif
             enddo
             if(nmu.gt.0) then 
                mu_2 = mu_2/nmu
             else
                mu_2 = LVT_rc%udef
             endif
             
             ens_std = 0 
             if(mu_2.ne.LVT_rc%udef) then 
                do m=1,nensem
                   ens_std = ens_std + (mu_m(m) - mu_2)**2
                enddo
                ens_std = sqrt(ens_std/nensem)
             else
                ens_std = LVT_rc%udef
             endif

             ci_val = LVT_rc%udef
             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                call LVT_computeCI(metric_tsdom(1:nsum_v),nsum_v,&
                     LVT_rc%pval_CI,ci_val)
             endif
             sum_v1 = 0 
             nsum_v1 = 0 
             maxv1 = max_param
             minv1 = min_param

             metric_tsdom = LVT_rc%udef
             
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   if(obs_ts(tid,k,l).ne.LVT_rc%udef) then 
                      sum_v1 = sum_v1 + obs_ts(tid,k,l)
                      nsum_v1 = nsum_v1 + 1
                      metric_tsdom(nsum_v1) = obs_ts(tid,k,l)
                      if(obs_ts(tid,k,l).lt.minv1) then 
                         minv1 = obs_ts(tid,k,l)
                      endif
                      if(obs_ts(tid,k,l).gt.maxv1) then 
                         maxv1 = obs_ts(tid,k,l) 
                      endif
                   endif
                endif
             enddo

             if(maxv1.eq.max_param) maxv1 = LVT_rc%udef
             if(minv1.eq.min_param) minv1 = LVT_rc%udef

             if(nsum_v1.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts).and.&
                  nsum_v1.gt.0) then 
                mean_v1 = sum_v1/nsum_v1
             else
                mean_v1 = LVT_rc%udef
             endif
                   
             sstd_v1 = 0 
             nsum_v1 = 0                       
             
             sstd_v1 = 0 
             nsum_v1 = 0 
             if(mean_v1.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1.) then 
                      if(obs_ts(tid,k,l).ne.LVT_rc%udef) then 
                         sstd_v1 = sstd_v1 + (obs_ts(tid,k,l)-mean_v1)**2
                         nsum_v1 = nsum_v1 + 1
                      endif
                   endif
                enddo
                sstd_v1 = sqrt(sstd_v1/nsum_v1)
             else
                sstd_v1 = LVT_rc%udef
             endif           
             if(nsum_v1.le.1) sstd_v1 = LVT_rc%udef
             ci_val1 = LVT_rc%udef  
             if(nsum_v1.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                call LVT_computeCI(metric_tsdom(1:nsum_v1),nsum_v1,&
                     LVT_rc%pval_CI,ci_val1)
             endif
!observation ensemble standard deviation is assumed to be undefined. 
             write(ftn_ts_loc(i),203,advance='no') &
                  mean_v, sstd_v, minv, maxv, ens_std,ci_val, &
                  mean_v1, sstd_v1, minv1, maxv1, LVT_rc%udef, ci_val1
          enddo
       enddo
    enddo

203 format(12E14.6)
  end subroutine writeTSinfo_noensem

!BOP
! 
! !ROUTINE: writeTSinfo3
! \label{writeTSinfo3}
!
! !INTERFACE: 
  subroutine writeTSinfo3(ftn_ts_loc,model,nsize,metric_ts,count_metric_ts,&
       dummy)
! 
! !USES:     

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: ftn_ts_loc(LVT_rc%ntslocs)
    type(LVT_metaDataEntry) :: model
    integer                 :: nsize
    real                    :: metric_ts(nsize,model%selectNlevs,&
         LVT_rc%nparam)
    integer                 :: count_metric_ts(nsize,model%selectNlevs,&
         LVT_rc%nparam)
    integer                 :: dummy


    real     :: sum_v
    integer  :: nsum_v
    real     :: wsum_v
    integer  :: wnsum_v
    real     :: mu_m(LVT_LIS_rc(1)%nensem)
    integer  :: nmu_m(LVT_LIS_rc(1)%nensem)
    real     :: mu_2
    integer  :: nmu
    real     :: ens_std
    real     :: mean_v, sstd_v
    real     :: ci_val
    integer  :: i,k,l,kk,tid
    integer  :: stid, m,t
    real     :: metric_tsdom(LVT_LIS_rc(1)%ntiles)
    real     :: maxv, minv
    integer  :: nensem

    nensem = LVT_LIS_rc(1)%nensem
    if(nsize.eq.LVT_rc%ngrid) then 
       nensem = 1
    endif

    do i=1, LVT_rc%ntslocs
       do k=1,model%selectNlevs
          do l=1,LVT_rc%nparam
             sum_v = 0 
             nsum_v = 0 
             wsum_v = 0 
             wnsum_v = 0 
             maxv = max_param
             minv = min_param

!works only if the tilespace info is passed down
#if 0              
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)*nensem
                   do m=1,nensem
                      t = stid+m
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         sum_v = sum_v + metric_ts(t,k,l)
                         nsum_v = nsum_v + 1
                         metric_tsdom(nsum_v) = metric_ts(t,k,l)
                         if(metric_ts(t,k,l).lt.minv) then 
                            minv = metric_ts(t,k,l)
                         endif
                         if(metric_ts(t,k,l).gt.maxv) then 
                            maxv = metric_ts(t,k,l)
                         endif
                      endif
                   enddo
                endif
             enddo
#endif
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   if(metric_ts(tid,k,l).ne.LVT_rc%udef) then 
                      sum_v = sum_v + metric_ts(tid,k,l)
                      nsum_v = nsum_v + 1
                      metric_tsdom(nsum_v) = metric_ts(tid,k,l)
                      if(metric_ts(tid,k,l).lt.minv) then 
                         minv = metric_ts(tid,k,l)
                      endif
                      if(metric_ts(tid,k,l).gt.maxv) then 
                         maxv = metric_ts(tid,k,l)
                      endif
                   endif
                endif
             enddo
             if(maxv.eq.max_param) maxv = LVT_rc%udef
             if(minv.eq.min_param) minv = LVT_rc%udef

             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                mean_v = sum_v/nsum_v
             else
                mean_v = LVT_rc%udef
             endif
             sstd_v = 0 
             nsum_v = 0 
#if 0 
             if(mean_v.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1) then 
                      stid = (tid-1)*nensem
                      do m=1,nensem
                         t = stid+m
                         if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                            sstd_v = sstd_v + &
                                 (metric_ts(t,k,l)-mean_v)**2
                            nsum_v = nsum_v + 1
                         endif
                      enddo
                   endif
                enddo
                sstd_v = sqrt(sstd_v/nsum_v)
             else
                sstd_v = LVT_rc%udef
             endif
#endif
             if(mean_v.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1) then 
                      if(metric_ts(tid,k,l).ne.LVT_rc%udef) then 
                         sstd_v = sstd_v + &
                              (metric_ts(tid,k,l)-mean_v)**2
                         nsum_v = nsum_v + 1
                      endif
                   endif
                enddo
                sstd_v = sqrt(sstd_v/nsum_v)
             else
                sstd_v = LVT_rc%udef
             endif
             if(nsum_v.le.1) then 
                sstd_v = LVT_rc%udef
             endif

             mu_m = 0 
             nmu_m = 0 
#if 0 
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)*nensem
                   do m=1,nensem
                      t = stid+m
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         mu_m(m) = mu_m(m) + metric_ts(t,k,l)
                         nmu_m(m) = nmu_m(m) + 1
                      endif
                   enddo
                end if
             enddo
#endif
             mu_2 = 0 
             nmu = 0 
             do m=1,nensem
                if(nmu_m(m).gt.0) then 
                   mu_m(m) = mu_m(m)/nmu_m(m)
                   mu_2 = mu_2 + mu_m(m)
                   nmu = nmu + 1
                endif
             enddo
             if(nmu.gt.0) then 
                mu_2 = mu_2/nmu
             else
                mu_2 = LVT_rc%udef
             endif
             
             ens_std = 0 
             if(mu_2.ne.LVT_rc%udef) then 
                do m=1,nensem
                   ens_std = ens_std + (mu_m(m) - mu_2)**2
                enddo
                ens_std = sqrt(ens_std/nensem)
             else
                ens_std = LVT_rc%udef
             endif

             ci_val = LVT_rc%udef
             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                call LVT_computeCI(metric_tsdom(1:nsum_v),nsum_v,&
                     LVT_rc%pval_CI,ci_val)
             endif
             write(ftn_ts_loc(i),204,advance='no') &
                  mean_v, sstd_v, minv, maxv, ens_std,ci_val
          enddo
       enddo
    enddo

204 format(6E14.6)
  end subroutine writeTSinfo3

!BOP
! 
! !ROUTINE: writeTSinfo4
! \label{writeTSinfo4}
!
! !INTERFACE: 
  subroutine writeTSinfo4(ftn_ts_loc,model,nsize_m,metric_ts,count_metric_ts,&
       nsize_o,obs_ts, count_obs_ts, obs_std_ts, count_obs_std_ts)
! 
! !USES:     

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This interface is used when both model and observation stats are to be
!  written and when observation standard deviation is specified. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: ftn_ts_loc(LVT_rc%ntslocs)
    type(LVT_metaDataEntry) :: model
    integer                 :: nsize_m
    real                    :: metric_ts(nsize_m,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: count_metric_ts(nsize_m,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: nsize_o
    real                    :: obs_ts(nsize_o,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: count_obs_ts(nsize_o,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    real                    :: obs_std_ts(nsize_o,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: count_obs_std_ts(nsize_o,model%selectNlevs,&
         LVT_rc%strat_nlevels)


    real     :: sum_v,sum_v1
    integer  :: nsum_v,nsum_v1
    real     :: mean_v, sstd_v, mean_v1,sstd_v1
    real     :: ens_std
    real     :: mu_m(LVT_LIS_rc(1)%nensem)
    integer  :: nmu_m(LVT_LIS_rc(1)%nensem)
    real     :: mu_2
    integer  :: nmu
    real     :: ci_val,ci_val1
    integer  :: i,k,l,kk,tid
    integer  :: stid, m,t
    real     :: metric_tsdom(LVT_LIS_rc(1)%ntiles)
    real     :: maxv, minv
    real     :: maxv1, minv1
    integer  :: nensem

    nensem = LVT_LIS_rc(1)%nensem
    if(nsize_m.eq.LVT_rc%ngrid) then 
       nensem = 1
    endif

    do i=1, LVT_rc%ntslocs
       do k=1,model%selectNlevs
          do l=1,LVT_rc%strat_nlevels
             sum_v = 0 
             nsum_v = 0 
             maxv = max_param
             minv = min_param

             do kk=1, LVT_TSobj(i)%npts             
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)*nensem
                   do m=1,nensem
                      t = stid+m
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         sum_v = sum_v + metric_ts(t,k,l)
                         nsum_v = nsum_v + 1
                         metric_tsdom(nsum_v) = metric_ts(t,k,l)
                         if(metric_ts(t,k,l).lt.minv) then 
                            minv = metric_ts(t,k,l)
                         endif
                         if(metric_ts(t,k,l).gt.maxv) then 
                            maxv = metric_ts(t,k,l)
                         endif
                      endif
                   enddo                   
                endif
             enddo

             if(maxv.eq.max_param) maxv = LVT_rc%udef
             if(minv.eq.min_param) minv = LVT_rc%udef

             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                mean_v = sum_v/nsum_v
             else
                mean_v = LVT_rc%udef
             endif
             
             sstd_v = 0 
             nsum_v = 0 

             if(mean_v.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1) then 
                      stid = (tid-1)*nensem
                      do m=1,nensem
                         t = stid+m
                         if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                            sstd_v = sstd_v + &
                                 (metric_ts(t,k,l)-mean_v)**2
                            nsum_v = nsum_v + 1                            
                         endif
                      enddo
                   endif
                enddo
                if(nsum_v.gt.0) then 
                   sstd_v = sqrt(sstd_v/nsum_v)
                else
                   sstd_v = LVT_rc%udef
                endif
             else
                sstd_v = LVT_rc%udef
             endif
             if(nsum_v.le.1) sstd_v = LVT_rc%udef

             mu_m = 0 
             nmu_m = 0 
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)*nensem
                   do m=1,nensem
                      t = stid+m
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         mu_m(m) = mu_m(m) + metric_ts(t,k,l)
                         nmu_m(m) = nmu_m(m) + 1
                      endif
                   enddo
                end if
             enddo
             mu_2 = 0 
             nmu = 0 
             do m=1,nensem
                if(nmu_m(m).gt.0) then 
                   mu_m(m) = mu_m(m)/nmu_m(m)
                   mu_2 = mu_2 + mu_m(m)
                   nmu = nmu + 1
                endif
             enddo
             if(nmu.gt.0) then 
                mu_2 = mu_2/nmu
             else
                mu_2 = LVT_rc%udef
             endif

             ens_std = 0 
             if(mu_2.ne.LVT_rc%udef) then 
                do m=1,nensem
                   ens_std = ens_std + (mu_m(m) - mu_2)**2
                enddo                
                ens_std = sqrt(ens_std/nensem)
             else
                ens_std = LVT_rc%udef
             endif

             ci_val = LVT_rc%udef
             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                call LVT_computeCI(metric_tsdom(1:nsum_v),nsum_v,&
                     LVT_rc%pval_CI,ci_val)
             endif
             sum_v1 = 0 
             nsum_v1 = 0 
             maxv1 = max_param
             minv1 = min_param
             metric_tsdom = LVT_rc%udef
             
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   if(obs_ts(tid,k,l).ne.LVT_rc%udef) then 
                      sum_v1 = sum_v1 + obs_ts(tid,k,l)
                      nsum_v1 = nsum_v1 + 1
                      metric_tsdom(nsum_v1) = obs_ts(tid,k,l)
                      if(obs_ts(tid,k,l).lt.minv1) then 
                         minv1 = obs_ts(tid,k,l)
                      endif
                      if(obs_ts(tid,k,l).gt.maxv1) then 
                         maxv1 = obs_ts(tid,k,l)
                      endif
                   endif
                endif
             enddo
            
             if(maxv1.eq.max_param) maxv1 = LVT_rc%udef
             if(minv1.eq.min_param) minv1 = LVT_rc%udef

             if(nsum_v1.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                mean_v1 = sum_v1/nsum_v1
             else
                mean_v1 = LVT_rc%udef
             endif
                   
             sstd_v1 = 0 
             nsum_v1 = 0                       
             
             sstd_v1 = 0 
             nsum_v1 = 0 
             if(mean_v1.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1.) then 
                      if(obs_std_ts(tid,k,l).ne.LVT_rc%udef) then 
                         sstd_v1 = sstd_v1 + obs_std_ts(tid,k,l)
                         nsum_v1 = nsum_v1 + 1
                      endif
                   endif
                enddo
                if(nsum_v1.gt.0) then
                   sstd_v1 = sstd_v1/nsum_v1
                else
                   sstd_v1 = LVT_rc%udef
                endif
             else
                sstd_v1 = LVT_rc%udef
             endif           
             ci_val1 = LVT_rc%udef  
             if(nsum_v1.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                call LVT_computeCI(metric_tsdom(1:nsum_v1),nsum_v1,&
                     LVT_rc%pval_CI,ci_val1)
             endif
!observation ensemble standard deviation is assumed to be undefined. 
             write(ftn_ts_loc(i),205,advance='no') &
                  mean_v, sstd_v, minv, maxv, ens_std, ci_val, &
                  mean_v1, sstd_v1, minv1,maxv1,LVT_rc%udef, ci_val1
          enddo
       enddo
    enddo

205 format(12E14.6)
  end subroutine writeTSinfo4

!BOP
! 
! !ROUTINE: writeTSinfo5
! \label{writeTSinfo5}
!
! !INTERFACE: 
  subroutine writeTSinfo5(ftn_ts_loc,model,nsize,metric_ts,count_metric_ts,&
       dummy1,dummy2)
! 
! !USES:     

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This is used to write TS info when the grid space variable is passed down
!  while working with an ensemble output
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: ftn_ts_loc(LVT_rc%ntslocs)
    type(LVT_metaDataEntry) :: model
    integer                 :: nsize
    real                    :: metric_ts(nsize,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: count_metric_ts(nsize,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: dummy1,dummy2

    real     :: sum_v
    integer  :: nsum_v
    real     :: wsum_v
    integer  :: wnsum_v
    real     :: ens_std
    real     :: mu_m(1)
    integer  :: nmu_m(1)
    real     :: mu_2
    integer  :: nmu
    real     :: mean_v, sstd_v
    real     :: ci_val
    integer  :: i,k,l,kk,tid
    integer  :: stid, m,t
    real     :: metric_tsdom(LVT_LIS_rc(1)%ntiles)
    real     :: maxv, minv
    integer  :: nensem

    nensem = LVT_LIS_rc(1)%nensem
    if(nsize.eq.LVT_rc%ngrid) then 
       nensem = 1
    endif

    do i=1, LVT_rc%ntslocs
       do k=1,model%selectNlevs
!       do k=1,model%selectStats
          do l=1,LVT_rc%strat_nlevels
             sum_v = 0 
             nsum_v = 0 
             wsum_v = 0 
             wnsum_v = 0 
             maxv = max_param
             minv = min_param

             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)
                   do m=1,1
                      t = stid+m
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         sum_v = sum_v + metric_ts(t,k,l)
                         nsum_v = nsum_v + 1
                         metric_tsdom(nsum_v) = metric_ts(t,k,l)
                         if(metric_ts(t,k,l).lt.minv) then 
                            minv = metric_ts(t,k,l)
                         endif
                         if(metric_ts(t,k,l).gt.maxv) then 
                            maxv = metric_ts(t,k,l)
                         endif
                      endif
                   enddo
                endif
             enddo

             if(maxv.eq.max_param) maxv = LVT_rc%udef
             if(minv.eq.min_param) minv = LVT_rc%udef

             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                mean_v = sum_v/nsum_v
             else
                mean_v = LVT_rc%udef
             endif

             sstd_v = 0 
             nsum_v = 0 

             if(mean_v.ne.LVT_rc%udef) then 
                do kk=1, LVT_TSobj(i)%npts
                   tid = LVT_Tsobj(i)%ts_tindex(kk)
                   if(tid.ne.-1) then 
                      stid = (tid-1)
                      do m=1,1 !nensem
                         t = stid+m
                         if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                            sstd_v = sstd_v + &
                                 (metric_ts(t,k,l)-mean_v)**2
                            nsum_v = nsum_v + 1
                         endif
                      enddo
                   endif
                enddo

                if(nsum_v.gt.0) then
                   sstd_v = sqrt(sstd_v/nsum_v)
                else
                   sstd_v = LVT_rc%udef
                endif
             else
                sstd_v = LVT_rc%udef
             endif
             if(nsum_v.le.1) then 
                sstd_v = LVT_rc%udef
             endif

             mu_m = 0 
             nmu_m = 0 
             do kk=1, LVT_TSobj(i)%npts
                tid = LVT_Tsobj(i)%ts_tindex(kk)
                if(tid.ne.-1) then 
                   stid = (tid-1)
                   do m=1,1
                      t = stid+m
                      if(metric_ts(t,k,l).ne.LVT_rc%udef) then 
                         mu_m(m) = mu_m(m) + metric_ts(t,k,l)
                         nmu_m(m) = nmu_m(m) + 1
                      endif
                   enddo
                end if
             enddo
             mu_2 = 0 
             nmu = 0 
             do m=1,1
                if(nmu_m(m).gt.0) then 
                   mu_m(m) = mu_m(m)/nmu_m(m)
                   mu_2 = mu_2 + mu_m(m)
                   nmu = nmu + 1
                endif
             enddo
             if(nmu.gt.0) then 
                mu_2 = mu_2/nmu
             else
                mu_2 = LVT_rc%udef
             endif
             
             ens_std = 0 
             if(mu_2.ne.LVT_rc%udef) then 
                do m=1,1
                   ens_std = ens_std + (mu_m(m) - mu_2)**2
                enddo
                ens_std = sqrt(ens_std)
             else
                ens_std = LVT_rc%udef
             endif             
             
             ci_val = LVT_rc%udef
             if(nsum_v.ge.(LVT_TSobj(i)%ts_min_pts*LVT_TSobj(i)%npts)) then 
                call LVT_computeCI(metric_tsdom(1:nsum_v),nsum_v,&
                     LVT_rc%pval_CI,ci_val)
             endif
             write(ftn_ts_loc(i),201,advance='no') &
                  mean_v, sstd_v, minv, maxv, ens_std,ci_val
          enddo
       enddo
    enddo

201 format(6E14.6)
  end subroutine writeTSinfo5

  
!BOP
! 
! !ROUTINE: writeDataBasedStrat1
! \label{writeDataBasedStrat1}
!
! !INTERFACE: 
  subroutine writeDataBasedStrat1(model,obs,stats,metric,nsize, &       
       metric_total)
! 
! !USES:     

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 

    type(LVT_metaDataEntry)  :: model
    type(LVT_metaDataEntry)  :: obs
    type(LVT_statsEntry)     :: stats
    type(LVT_metricEntry)     :: metric
    integer                 :: nsize
    real                    :: metric_total(nsize,model%selectNlevs,&
         LVT_rc%strat_nlevels)


    integer, allocatable :: nbin(:)
    real, allocatable :: bin(:)
    
    integer        :: ftn,i,t,k,kk,binval
    character*500  :: filename

    if(LVT_rc%data_based_strat.eq.1) then       
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then  
          do i=1,LVT_rc%data_based_nstrats
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(metric%short_name)//'_by_'//&
                  trim(LVT_rc%data_based_strat_var(i))//'_'//&
                  trim(stats%short_name)//'.dat'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=(filename),form='formatted')
             
             allocate(bin(LVT_rc%data_based_strat_nbins(i)))
             allocate(nbin(LVT_rc%data_based_strat_nbins(i)))
             bin = 0 
             nbin = 0 
             
             do t=1,nsize
                do k=1,model%selectNlevs
                   if(metric_total(t,k,1).ne.LVT_rc%udef) then 
                      binval = nint((LVT_rc%strat_data(i,t) & 
                           - LVT_rc%data_based_strat_min(i))/&
                           LVT_rc%data_based_strat_delta(i)) + 1
                      if(binval.le.0) binval = 1
                      if(binval.gt.LVT_rc%data_based_strat_nbins(i)) &
                           binval = LVT_rc%data_based_strat_nbins(i)
                      bin(binval) = bin(binval) + metric_total(t,k,1)
                      nbin(binval) = nbin(binval) + 1
                   endif
                enddo
             enddo
             
             do kk=1,LVT_rc%data_based_strat_nbins(i)
                if(nbin(kk).ne.0) then 
                   write(ftn,204) &
                        LVT_rc%data_based_strat_min(i)+(kk-1)*&
                        LVT_rc%data_based_strat_delta(i),bin(kk)/nbin(kk), &
                        nbin(kk)
                else
                   write(ftn,204) &
                        LVT_rc%data_based_strat_min(i)+(kk-1)*&
                        LVT_rc%data_based_strat_delta(i),LVT_rc%udef,&
                        nbin(kk)
                endif
             enddo
             
             deallocate(bin)
             deallocate(nbin)
             call LVT_releaseUnitNumber(ftn)
             
          enddo
       endif
    endif
204 format(2E14.6,I14.3)
  end subroutine WriteDataBasedStrat1

!BOP
! 
! !ROUTINE: writeDataBasedStrat2
! \label{writeDataBasedStrat2}
!
! !INTERFACE: 
  subroutine writeDataBasedStrat2(model,obs,stats,metric, &       
       nsize_m,metric_total,nsize_o,obs_total)
! 
! !USES:     

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 

    type(LVT_metaDataEntry)  :: model
    type(LVT_metaDataEntry)  :: obs
    type(LVT_statsEntry)     :: stats
    type(LVT_metricEntry)     :: metric
    integer                 :: nsize_m
    real                    :: metric_total(nsize_m,model%selectNlevs,&
         LVT_rc%strat_nlevels)
    integer                 :: nsize_o
    real                    :: obs_total(nsize_o,model%selectNlevs,&
         LVT_rc%strat_nlevels)


    integer, allocatable :: nbin(:)
    real, allocatable :: bin(:)
    integer, allocatable :: nobin(:)
    real, allocatable :: obin(:)
    
    integer        :: ftn,i,t,k,kk,binval
    character*500  :: filename

    if(LVT_rc%data_based_strat.eq.1) then       
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then  
          do i=1,LVT_rc%data_based_nstrats
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(metric%short_name)//'_by_'//&
                  trim(LVT_rc%data_based_strat_var(i))//'_'//&
                  trim(stats%short_name)//'.dat'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=(filename),form='formatted')
             
             allocate(bin(LVT_rc%data_based_strat_nbins(i)))
             allocate(nbin(LVT_rc%data_based_strat_nbins(i)))
             bin = 0 
             nbin = 0 
             allocate(obin(LVT_rc%data_based_strat_nbins(i)))
             allocate(nobin(LVT_rc%data_based_strat_nbins(i)))
             obin = 0 
             nobin = 0 
             
             do t=1,nsize_m
                do k=1,model%selectNlevs
                   if(metric_total(t,k,1).ne.LVT_rc%udef) then 
                      binval = nint((LVT_rc%strat_data(i,t) & 
                           - LVT_rc%data_based_strat_min(i))/&
                           LVT_rc%data_based_strat_delta(i)) + 1
                      if(binval.le.0) binval = 1
                      if(binval.gt.LVT_rc%data_based_strat_nbins(i)) &
                           binval = LVT_rc%data_based_strat_nbins(i)
                      bin(binval) = bin(binval) + metric_total(t,k,1)
                      nbin(binval) = nbin(binval) + 1
                   endif
                enddo
             enddo
             do t=1,nsize_o
                do k=1,obs%selectNlevs
                   if(obs_total(t,k,1).ne.LVT_rc%udef) then 
                      binval = nint((LVT_rc%strat_data(i,t) & 
                           - LVT_rc%data_based_strat_min(i))/&
                           LVT_rc%data_based_strat_delta(i)) + 1
                      if(binval.le.0) binval = 1
                      if(binval.gt.LVT_rc%data_based_strat_nbins(i)) &
                           binval = LVT_rc%data_based_strat_nbins(i)
                      obin(binval) = obin(binval) + obs_total(t,k,1)
                      nobin(binval) = nobin(binval) + 1
                   endif                   
                enddo
             enddo

             do kk=1,LVT_rc%data_based_strat_nbins(i)
                if(nbin(kk).ne.0) then 
                   bin(kk) = bin(kk)/nbin(kk)
                else
                   bin(kk) = LVT_rc%udef
                endif
                if(nobin(kk).ne.0) then 
                   obin(kk) = obin(kk)/nobin(kk)
                else
                   obin(kk) = LVT_rc%udef
                endif
                write(ftn,205) &
                     LVT_rc%data_based_strat_min(i)+(kk-1)*&
                     LVT_rc%data_based_strat_delta(i),bin(kk), &
                     nbin(kk),obin(kk),nobin(kk)
             enddo
             
             deallocate(bin)
             deallocate(nbin)
             deallocate(obin)
             deallocate(nobin)
             call LVT_releaseUnitNumber(ftn)
             
          enddo
       endif
    endif
205 format(2E14.6,I14.3,E14.6,I14.3)
  end subroutine WriteDataBasedStrat2


!BOP
! 
! !ROUTINE: writeSeasonalCycleInfo1
! \label(writeSeasonalCycleInfo1)
!
! !INTERFACE:
  subroutine writeSeasonalCycleInfo1(model,obs,stats,metric,nsize,&
       metric_asc,&
       count_metric_asc)
! 
! !USES:   

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    
    type(LVT_metaDataEntry)   :: model
    type(LVT_metaDataEntry)   :: obs
    type(LVT_statsEntry)      :: stats
    type(LVT_metricEntry)     :: metric
    integer                   :: nsize
    real                      :: metric_asc(nsize,model%selectNlevs,&
         LVT_rc%nasc)
    integer                   :: count_metric_asc(nsize,model%selectNlevs,&
         LVT_rc%nasc)

    real     :: sum_v_sc(LVT_rc%nasc)
    integer  :: nsum_v_sc(LVT_rc%nasc)
    real     :: wsum_v_sc(LVT_rc%nasc)
    integer  :: wnsum_v_sc(LVT_rc%nasc)
    integer                   :: tid,k,l,i,kk,ftn
    character*500             :: filename

    if(metric%computeSC.eq.1.and.metric%timeOpt.eq.1) then
       if(stats%selectOpt.eq.1.and.model%selectNlevs.ge.1) then 
                    
          do i=1, LVT_rc%ntslocs
             ftn = LVT_getNextUnitNumber()
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(metric%short_name)//'_ASC_'//&
                  trim(LVT_TSobj(i)%tslocname)//'_'//&
                  trim(stats%short_name)//'.dat'
             open(ftn,file=(filename),form='formatted')
             
             do k=1,model%selectNlevs
                sum_v_sc = 0 
                nsum_v_sc = 0 
                wsum_v_sc = 0 
                wnsum_v_sc = 0 
                
                do l=1,LVT_rc%nasc                   
                   do kk=1,LVT_TSobj(i)%npts
                      tid = LVT_Tsobj(i)%ts_tindex(kk)
                      if(tid.ne.-1) then 
                         if(metric_asc(tid,k,l).ne.LVT_rc%udef) then 
                            sum_v_sc(l) = sum_v_sc(l) + metric_asc(tid,k,l)
                            nsum_v_sc(l) = nsum_v_sc(l)+ 1
!                            wnsum_v_sc(l) = wnsum_v_sc(l) + &
!                                 metric_asc(tid,k,l)*&
!                                 count_metric_asc(tid,k,l)
!                            wnsum_v_sc(l) = wnsum_v_sc(l) + &
!                                 count_metric_asc(tid,k,l)
                            
                         endif
                      endif
                   enddo
                enddo
                
                do l=1,LVT_rc%nasc
                   if(nsum_v_sc(l).gt.0.and.nsum_v_sc(l).gt.&
                        LVT_rc%scCountThreshold) then 
                      sum_v_sc(l) = sum_v_sc(l)/nsum_v_sc(l)
!                      wsum_v_sc(l) = wsum_v_sc(l)/wnsum_v_sc(l)
                   else
                      sum_v_sc(l) = LVT_rc%udef
!                      wsum_v_sc(l) = LVT_rc%udef
                   endif
                   write(ftn,222) l, sum_v_sc(l), nsum_v_sc(l)
                enddo
             enddo
             call LVT_releaseUnitNumber(ftn)
          enddo
       endif
    endif
222 format(I2.2,E14.6,I14.6)

  end subroutine WriteSeasonalCycleInfo1

!BOP
! 
! !ROUTINE: writeSeasonalCycleInfo2
! \label(writeSeasonalCycleInfo2)
!
! !INTERFACE:
  subroutine writeSeasonalCycleInfo2(model,obs,stats,metric,nsize_m,&
       metric_asc,&
       count_metric_asc,nsize_o,obs_asc,count_obs_asc)
! 
! !USES:   

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    
    type(LVT_metaDataEntry)   :: model
    type(LVT_metaDataEntry)   :: obs
    type(LVT_statsEntry)      :: stats
    type(LVT_metricEntry)     :: metric
    integer                   :: nsize_m
    real                      :: metric_asc(nsize_m,model%selectNlevs,&
         LVT_rc%nasc)
    integer                   :: count_metric_asc(nsize_m,model%selectNlevs,&
         LVT_rc%nasc)
    integer                   :: nsize_o
    real                      :: obs_asc(nsize_o,model%selectNlevs,&
         LVT_rc%nasc)
    integer                   :: count_obs_asc(nsize_o,model%selectNlevs,&
         LVT_rc%nasc)

    real     :: sum_v_sc(LVT_rc%nasc)
    integer  :: nsum_v_sc(LVT_rc%nasc)
    real     :: sum_v_sc1(LVT_rc%nasc)
    integer  :: nsum_v_sc1(LVT_rc%nasc)
    real     :: wsum_v_sc(LVT_rc%nasc)
    integer  :: wnsum_v_sc(LVT_rc%nasc)
    integer                   :: tid,k,l,i,kk,ftn
    character*500             :: filename

    if(metric%computeSC.eq.1.and.metric%timeOpt.eq.1) then
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
                    
          do i=1, LVT_rc%ntslocs
             ftn = LVT_getNextUnitNumber()
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(metric%short_name)//'_ASC_'//&
                  trim(LVT_TSobj(i)%tslocname)//'_'//&
                  trim(stats%short_name)//'.dat'
             open(ftn,file=(filename),form='formatted')
             
             do k=1,model%selectNlevs
                sum_v_sc = 0 
                nsum_v_sc = 0 
                wsum_v_sc = 0 
                wnsum_v_sc = 0 
                
                sum_v_sc1 = 0
                nsum_v_sc1 = 0 
                
                do l=1,LVT_rc%nasc                   
                   do kk=1,LVT_TSobj(i)%npts
                      tid = LVT_Tsobj(i)%ts_tindex(kk)
                      if(tid.ne.-1) then 
                         if(metric_asc(tid,k,l).ne.LVT_rc%udef) then 
                            sum_v_sc(l) = sum_v_sc(l) + metric_asc(tid,k,l)
                            nsum_v_sc(l) = nsum_v_sc(l)+ 1
!                            wnsum_v_sc(l) = wnsum_v_sc(l) + &
!                                 metric_asc(tid,k,l)*&
!                                 count_metric_asc(tid,k,l)
!                            wnsum_v_sc(l) = wnsum_v_sc(l) + &
!                                 count_metric_asc(tid,k,l)
                            
                         endif
                      endif
                   enddo
                enddo

                do l=1,LVT_rc%nasc                   
                   do kk=1,LVT_TSobj(i)%npts
                      tid = LVT_Tsobj(i)%ts_tindex(kk)
                      if(tid.ne.-1) then 
                         if(obs_asc(tid,k,l).ne.LVT_rc%udef) then 
                            sum_v_sc1(l) = sum_v_sc1(l) + obs_asc(tid,k,l)
                            nsum_v_sc1(l) = nsum_v_sc1(l)+ 1
                         endif
                      endif
                   enddo
                enddo
                
                do l=1,LVT_rc%nasc
                   if(nsum_v_sc(l).gt.0) then !.and.nsum_v_sc(l).gt.&
!                        LVT_rc%scCountThreshold) then 
                      sum_v_sc(l) = sum_v_sc(l)/nsum_v_sc(l)
!                      wsum_v_sc(l) = wsum_v_sc(l)/wnsum_v_sc(l)
                   else
                      sum_v_sc(l) = LVT_rc%udef
!                      wsum_v_sc(l) = LVT_rc%udef
                   endif
                   if(nsum_v_sc1(l).gt.0.and.nsum_v_sc1(l).gt.&
                        LVT_rc%scCountThreshold) then 
                      sum_v_sc1(l) = sum_v_sc1(l)/nsum_v_sc1(l)
                   else
                      sum_v_sc1(l) = LVT_rc%udef
                   endif
                   write(ftn,222) l, sum_v_sc(l), nsum_v_sc(l), &
                        sum_v_sc1(l),nsum_v_sc1(l)
                enddo
             enddo
             call LVT_releaseUnitNumber(ftn)
          enddo
       endif
    endif
222 format(I2.2,E14.6,I14.3,E14.6,I14.3)

  end subroutine WriteSeasonalCycleInfo2

!BOP
! 
! !ROUTINE: writeAvgDiurnalCycleInfo1
! \label(writeAvgDiurnalCycleInfo1)
!
! !INTERFACE:
  subroutine writeAvgDiurnalCycleInfo1(model,obs,stats,metric,&
       nsize,metric_adc,count_metric_adc)
! 
! !USES:   

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    
    type(LVT_metaDataEntry)  :: model
    type(LVT_metaDataEntry)  :: obs
    type(LVT_statsEntry)     :: stats
    type(LVT_metricEntry)    :: metric
    integer                  :: nsize
    real                     :: metric_adc(nsize,model%selectNlevs,&
         LVT_rc%nadc)
    integer                  :: count_metric_adc(nsize,model%selectNlevs,&
         LVT_rc%nadc)

    real     :: sum_v_adc(LVT_rc%nadc)
    integer  :: nsum_v_adc(LVT_rc%nadc)
    real     :: wsum_v_adc(LVT_rc%nadc)
    integer  :: wnsum_v_adc(LVT_rc%nadc)
    integer                   :: tid,k,l,i,kk,ftn
    character*500             :: filename

    if(metric%computeADC.eq.1.and.metric%timeOpt.eq.1) then
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
                    
          do i=1, LVT_rc%ntslocs
             ftn = LVT_getNextUnitNumber()
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(metric%short_name)//'_ADC_'//&
                  trim(LVT_TSobj(i)%tslocname)//'_'//&
                  trim(stats%short_name)//'.dat'
             open(ftn,file=(filename),form='formatted')
             
             do k=1,model%selectNlevs
                sum_v_adc = 0 
                nsum_v_adc = 0 
                wsum_v_adc = 0 
                wnsum_v_adc = 0 
                
                do l=1,LVT_rc%nadc                   
                   do kk=1,LVT_TSobj(i)%npts
                      tid = LVT_Tsobj(i)%ts_tindex(kk)
                      if(tid.ne.-1) then 
                         if(metric_adc(tid,k,l).ne.LVT_rc%udef) then 
                            sum_v_adc(l) = sum_v_adc(l) + metric_adc(tid,k,l)
                            nsum_v_adc(l) = nsum_v_adc(l)+ 1
!                            wnsum_v_adc(l) = wnsum_v_adc(l) + &
!                                 metric_adc(tid,k,l)*&
!                                 count_metric_adc(tid,k,l)
!                            wnsum_v_adc(l) = wnsum_v_adc(l) + &
!                                 count_metric_adc(tid,k,l)
                            
                         endif
                      endif
                   enddo
                enddo
                
                do l=1,LVT_rc%nadc
                   if(nsum_v_adc(l).gt.0.and.nsum_v_adc(l).gt.&
                        LVT_rc%adcCountThreshold) then 
                      sum_v_adc(l) = sum_v_adc(l)/nsum_v_adc(l)
!                      wsum_v_adc(l) = wsum_v_adc(l)/wnsum_v_adc(l)
                   else
                      sum_v_adc(l) = LVT_rc%udef
!                      wsum_v_adc(l) = LVT_rc%udef
                   endif
                   write(ftn,333) l*LVT_rc%statswriteint/3600,&
                        sum_v_adc(l), nsum_v_adc(l)
                enddo
             enddo
             call LVT_releaseUnitNumber(ftn)
          enddo
       endif
    endif
333 format(I2.2,E14.6,I14.3)

  end subroutine WriteAvgDiurnalCycleInfo1

!BOP
! 
! !ROUTINE: writeAvgDiurnalCycleInfo2
! \label(writeAvgDiurnalCycleInfo2)
!
! !INTERFACE:
  subroutine writeAvgDiurnalCycleInfo2(model,obs,stats,metric,&
       nsize_m,metric_adc,count_metric_adc,nsize_o,obs_adc,count_obs_adc)
! 
! !USES:   

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    
    type(LVT_metaDataEntry)  :: model
    type(LVT_metaDataEntry)  :: obs
    type(LVT_statsEntry)     :: stats
    type(LVT_metricEntry)    :: metric
    integer                  :: nsize_m
    real                     :: metric_adc(nsize_m,model%selectNlevs,&
         LVT_rc%nadc)
    integer                  :: count_metric_adc(nsize_m,model%selectNlevs,&
         LVT_rc%nadc)
    integer                  :: nsize_o
    real                     :: obs_adc(nsize_o,model%selectNlevs,&
         LVT_rc%nadc)
    integer                  :: count_obs_adc(nsize_o,model%selectNlevs,&
         LVT_rc%nadc)

    real     :: sum_v_adc(LVT_rc%nadc)
    integer  :: nsum_v_adc(LVT_rc%nadc)
    real     :: sum_v_adc1(LVT_rc%nadc)
    integer  :: nsum_v_adc1(LVT_rc%nadc)
    real     :: wsum_v_adc(LVT_rc%nadc)
    integer  :: wnsum_v_adc(LVT_rc%nadc)
    integer                   :: tid,k,l,i,kk,ftn
    character*500             :: filename

    if(metric%computeADC.eq.1.and.metric%timeOpt.eq.1) then
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
                    
          do i=1, LVT_rc%ntslocs
             ftn = LVT_getNextUnitNumber()
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(metric%short_name)//'_ADC_'//&
                  trim(LVT_TSobj(i)%tslocname)//'_'//&
                  trim(stats%short_name)//'.dat'
             open(ftn,file=(filename),form='formatted')
             
             do k=1,model%selectNlevs
                sum_v_adc = 0 
                nsum_v_adc = 0 
                sum_v_adc1 = 0 
                nsum_v_adc1 = 0 
                wsum_v_adc = 0 
                wnsum_v_adc = 0 
                
                do l=1,LVT_rc%nadc                   
                   do kk=1,LVT_TSobj(i)%npts
                      tid = LVT_Tsobj(i)%ts_tindex(kk)
                      if(tid.ne.-1) then 
                         if(metric_adc(tid,k,l).ne.LVT_rc%udef) then 
                            sum_v_adc(l) = sum_v_adc(l) + metric_adc(tid,k,l)
                            nsum_v_adc(l) = nsum_v_adc(l)+ 1
!                            wnsum_v_adc(l) = wnsum_v_adc(l) + &
!                                 metric_adc(tid,k,l)*&
!                                 count_metric_adc(tid,k,l)
!                            wnsum_v_adc(l) = wnsum_v_adc(l) + &
!                                 count_metric_adc(tid,k,l)
                            
                         endif
                         if(obs_adc(tid,k,l).ne.LVT_rc%udef) then 
                            sum_v_adc1(l) = sum_v_adc1(l) + obs_adc(tid,k,l)
                            nsum_v_adc1(l) = nsum_v_adc1(l)+ 1
                         endif
                      endif
                   enddo
                enddo
                
                do l=1,LVT_rc%nadc
                   if(nsum_v_adc(l).gt.0.and.nsum_v_adc(l).gt.&
                        LVT_rc%adcCountThreshold) then 
                      sum_v_adc(l) = sum_v_adc(l)/nsum_v_adc(l)
!                      wsum_v_adc(l) = wsum_v_adc(l)/wnsum_v_adc(l)
                   else
                      sum_v_adc(l) = LVT_rc%udef
!                      wsum_v_adc(l) = LVT_rc%udef
                   endif
                   if(nsum_v_adc1(l).gt.0.and.nsum_v_adc1(l).gt.&
                        LVT_rc%adcCountThreshold) then 
                      sum_v_adc1(l) = sum_v_adc1(l)/nsum_v_adc1(l)
                   else
                      sum_v_adc1(l) = LVT_rc%udef
                   endif
                   write(ftn,334) l*LVT_rc%statswriteint/3600,&
                        sum_v_adc(l), nsum_v_adc(l), &
                        sum_v_adc1(l), nsum_v_adc1(l)
                enddo
             enddo
             call LVT_releaseUnitNumber(ftn)
          enddo
       endif
    endif
334 format(I2.2,E14.6,I14.3,E14.6,I14.3)

  end subroutine WriteAvgDiurnalCycleInfo2


end module LVT_TSMod
