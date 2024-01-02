!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah39_dynsetup
! \label{noah39_dynsetup}
!
! !REVISION HISTORY:
!  28 Jan 2002: Sujay Kumar, Initial Specification
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  14 Jan 2014: David Mocko, reconfirmed Noah3.3 in LIS7.0
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah39_dynsetup(n)
! !USES: 
  use LIS_coreMod, only : LIS_rc,LIS_domain, LIS_surface
  use LIS_snowMod, only : LIS_snow_struc
  use LIS_vegDataMod,    only : LIS_gfrac, LIS_lai, LIS_roughness
  use LIS_emissMod,      only : LIS_emiss
  use LIS_albedoMod,     only : LIS_alb
  use noah39_lsmMod, only : noah39_struc
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use module_sf_noah36lsm, only : month_d, calc_localtime
! 
! !DESCRIPTION: 
!  This routine sets up the time-dependent variables in Noah-3.9
! 
!EOP     
  implicit none
  integer, intent(in) :: n 

  integer   :: tid,i
  integer   :: t,gid,change,local_hour
  integer   :: locdoy,locyr,locmo,locda,lochr,locmn,locss
  real*8    :: loctime
  real      :: interp_fraction
  real      :: locgmt
  integer   :: col,row,sftype
  integer   :: ncount(LIS_rc%ngrid(n))
  logical   :: Bondvillecheck

  Bondvillecheck = .false.
  do i=1,LIS_rc%nmetforc
     if (trim(LIS_rc%metforc(i)).eq."Bondville") Bondvillecheck = .true.
  enddo

  if(LIS_rc%snowsrc(n).gt.0) then
     ncount = 0 
     LIS_snow_struc(n)%sneqv = 0 
     LIS_snow_struc(n)%snowdepth = 0 

#if 0 
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
        row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
        gid = LIS_domain(n)%gindex(col,row)
        LIS_snow_struc(n)%sneqv(gid)     = LIS_snow_struc(n)%sneqv(gid) + & 
             noah39_struc(n)%noah(t)%sneqv
        LIS_snow_struc(n)%snowdepth(gid) = LIS_snow_struc(n)%snowdepth(gid) + & 
             noah39_struc(n)%noah(t)%snowh
        ncount(gid) = ncount(gid)+1 
     enddo

     do t=1,LIS_rc%ngrid(n)
        LIS_snow_struc(n)%sneqv(t) = LIS_snow_struc(n)%sneqv(t)/ncount(t)
        LIS_snow_struc(n)%snowdepth(t) = &
             LIS_snow_struc(n)%snowdepth(t)/ncount(t)
     enddo
#endif
     LIS_snow_struc(n)%sneqv(:)     = 0.0
     LIS_snow_struc(n)%snowdepth(:) = 0.0
     
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
        LIS_snow_struc(n)%sneqv(tid)     = LIS_snow_struc(n)%sneqv(tid) + & 
             noah39_struc(n)%noah(t)%sneqv
     enddo
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        LIS_snow_struc(n)%snowdepth(gid) = LIS_snow_struc(n)%snowdepth(gid) + & 
                noah39_struc(n)%noah(t)%snowh
        ncount(gid) = ncount(gid) + 1
     enddo

     do t=1,LIS_rc%ngrid(n)
        if(ncount(t).gt.0) then 
           LIS_snow_struc(n)%snowdepth(t) = &
                LIS_snow_struc(n)%snowdepth(t)/ncount(t)
        else
           LIS_snow_struc(n)%snowdepth(t) = 0.0
        endif
     enddo

  endif

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (LIS_rc%usegreennessmap(n).eq."none") then
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        call calc_localtime(LIS_rc%hr,LIS_domain(n)%grid(gid)%lon,     &
                            local_hour,change)
! For a true benchmark against the Noah-3.9 testcase from NCAR,
! set "change = 0".  This line allows the code to _incorrectly_
! run in the same way as the Bondville testcase, which runs on
! local time instead of on UTC time. - dmm
        if (Bondvillecheck) change = 0
        change = change * 3600
        locyr = LIS_rc%yr
        locmo = LIS_rc%mo
        locda = LIS_rc%da
        lochr = LIS_rc%hr
        locmn = LIS_rc%mn
        locss = LIS_rc%ss
        call LIS_date2time(loctime,locdoy,locgmt,locyr,locmo,          &
                           locda,lochr,locmn,locss)
        call LIS_tick(loctime,locdoy,locgmt,locyr,locmo,               &
                      locda,lochr,locmn,locss,real(change))
        noah39_struc(n)%noah(t)%shdfac =                               &
                            month_d(noah39_struc(n)%shdfac_monthly,    &
                            locyr,locmo,locda)
     else
        tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
        noah39_struc(n)%noah(t)%shdfac =  LIS_gfrac(n)%greenness(tid)
     endif

     if (noah39_struc(n)%noah(t)%vegt.eq.LIS_rc%bareclass)             &
              noah39_struc(n)%noah(t)%shdfac = 0.0
     if (noah39_struc(n)%noah(t)%vegt.eq.LIS_rc%snowclass)             &
              noah39_struc(n)%noah(t)%shdfac = 0.0
     if (noah39_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass)            &
              noah39_struc(n)%noah(t)%shdfac = 0.0

     noah39_struc(n)%noah(t)%z0brd_old = noah39_struc(n)%noah(t)%z0brd

     if(noah39_struc(n)%noah(t)%shdfac.ge.noah39_struc(n)%noah(t)%shdmax) then
        if(LIS_rc%useemissmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah39_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
        else
           if(noah39_struc(n)%embrd_upd.eq.0) then
              noah39_struc(n)%noah(t)%embrd =                             &
                   noah39_struc(n)%emissmax(noah39_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah39_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if(noah39_struc(n)%z0brd_upd.eq.0) then
              noah39_struc(n)%noah(t)%z0brd =                             &
                   noah39_struc(n)%z0max(noah39_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
           if(noah39_struc(n)%alb_upd.eq.0) then
              noah39_struc(n)%noah(t)%alb =                            &
                    noah39_struc(n)%albmin(noah39_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
           if(noah39_struc(n)%lai_upd.eq.0) then
              noah39_struc(n)%noah(t)%lai =                            &
                    noah39_struc(n)%laimax(noah39_struc(n)%noah(t)%vegt)          
           else
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
              noah39_struc(n)%noah(t)%lai = &
                   LIS_lai(n)%tlai(tid)
           endif
        endif

     elseif(noah39_struc(n)%noah(t)%shdfac.le.noah39_struc(n)%noah(t)%shdmin) then
        if(LIS_rc%useemissmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah39_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
        else
           if(noah39_struc(n)%embrd_upd.eq.0) then
              noah39_struc(n)%noah(t)%embrd =                             &
                   noah39_struc(n)%emissmin(noah39_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah39_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if(noah39_struc(n)%z0brd_upd.eq.0) then
              noah39_struc(n)%noah(t)%z0brd =                             &
                   noah39_struc(n)%z0min(noah39_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
           if(noah39_struc(n)%alb_upd.eq.0) then
              noah39_struc(n)%noah(t)%alb =                            &
                    noah39_struc(n)%albmax(noah39_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
           if(noah39_struc(n)%lai_upd.eq.0) then
              noah39_struc(n)%noah(t)%lai =                            &
                    noah39_struc(n)%laimin(noah39_struc(n)%noah(t)%vegt)
           else
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
              noah39_struc(n)%noah(t)%lai = &
                   LIS_lai(n)%tlai(tid)
           endif
        endif

     else
        if(noah39_struc(n)%noah(t)%shdmax.gt.noah39_struc(n)%noah(t)%shdmin) then
           interp_fraction =                                           &
                (noah39_struc(n)%noah(t)%shdfac-noah39_struc(n)%noah(t)%shdmin) / &
                (noah39_struc(n)%noah(t)%shdmax-noah39_struc(n)%noah(t)%shdmin)
           interp_fraction = min(interp_fraction,1.0)
           interp_fraction = max(interp_fraction,0.0)

           if(LIS_rc%useemissmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah39_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              if(noah39_struc(n)%embrd_upd.eq.0) then
                 noah39_struc(n)%noah(t)%embrd = ((1.0 - interp_fraction) * &
                      noah39_struc(n)%emissmin(noah39_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah39_struc(n)%emissmax(noah39_struc(n)%noah(t)%vegt))
              endif
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah39_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              if(noah39_struc(n)%z0brd_upd.eq.0) then
                 noah39_struc(n)%noah(t)%z0brd = ((1.0 - interp_fraction) * &
                      noah39_struc(n)%z0min(noah39_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah39_struc(n)%z0max(noah39_struc(n)%noah(t)%vegt))
              endif
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              if(noah39_struc(n)%lai_upd.eq.0) then
                 noah39_struc(n)%noah(t)%lai = ((1.0-interp_fraction) * &
                      noah39_struc(n)%laimin(noah39_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah39_struc(n)%laimax(noah39_struc(n)%noah(t)%vegt))
              else
                 tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
                 noah39_struc(n)%noah(t)%lai = &
                      LIS_lai(n)%tlai(tid)
              endif
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              if(noah39_struc(n)%alb_upd.eq.0) then
                 noah39_struc(n)%noah(t)%alb = ((1.0-interp_fraction) * &
                      noah39_struc(n)%albmax(noah39_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah39_struc(n)%albmin(noah39_struc(n)%noah(t)%vegt))
              endif
           endif
        else
           if(LIS_rc%useemissmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah39_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              if(noah39_struc(n)%embrd_upd.eq.0) then
                 noah39_struc(n)%noah(t)%embrd =                          &
                      0.5 * noah39_struc(n)%emissmin(noah39_struc(n)%noah(t)%vegt) + &
                      0.5 * noah39_struc(n)%emissmax(noah39_struc(n)%noah(t)%vegt)
              endif
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah39_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              if(noah39_struc(n)%z0brd_upd.eq.0) then
                 noah39_struc(n)%noah(t)%z0brd =                          &
                      0.5 * noah39_struc(n)%z0min(noah39_struc(n)%noah(t)%vegt) + &
                      0.5 * noah39_struc(n)%z0max(noah39_struc(n)%noah(t)%vegt)
              endif
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              if(noah39_struc(n)%lai_upd.eq.0) then
                 noah39_struc(n)%noah(t)%lai =                         &
                      0.5 * noah39_struc(n)%laimin(noah39_struc(n)%noah(t)%vegt) + &
                      0.5 * noah39_struc(n)%laimax(noah39_struc(n)%noah(t)%vegt)
              else
                 tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
                 noah39_struc(n)%noah(t)%lai = &
                      LIS_lai(n)%tlai(tid)
              endif
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah39_struc(n)%noah(t)%alb =                            &
                   0.5 * noah39_struc(n)%albmin(noah39_struc(n)%noah(t)%vegt) + &
                   0.5 * noah39_struc(n)%albmax(noah39_struc(n)%noah(t)%vegt)
           endif
        endif
     endif

     if (LIS_rc%usealbedomap(n).ne."none") then
        if(noah39_struc(n)%alb_upd.eq.0) then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
           noah39_struc(n)%noah(t)%alb = LIS_alb(n)%albsf(tid)
        endif
     endif
     if (LIS_rc%uselaimap(n).ne."none") then
        if(noah39_struc(n)%lai_upd.eq.0) then 
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
           noah39_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
        endif
     endif
  enddo

  return
end subroutine noah39_dynsetup
