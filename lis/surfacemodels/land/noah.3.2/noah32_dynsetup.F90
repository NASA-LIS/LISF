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
! !ROUTINE: noah32_dynsetup
! \label{noah32_dynsetup}
!
! !REVISION HISTORY:
!  28 Jan 2002: Sujay Kumar, Initial Specification
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!
! !INTERFACE:
subroutine noah32_dynsetup(n)
! !USES: 
  use LIS_coreMod, only : LIS_rc,LIS_domain, LIS_surface
  use LIS_snowMod, only : LIS_snow_struc
  use LIS_vegDataMod,    only : LIS_gfrac, LIS_lai, LIS_roughness
  use LIS_emissMod,      only : LIS_emiss
  use LIS_albedoMod,     only : LIS_alb
  use noah32_lsmMod, only : noah32_struc
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use module_sf_noah32lsm, only : month_d, calc_localtime
! 
! !DESCRIPTION: 
!  This routine sets up the time-dependent variables in Noah3.2
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
  integer   :: col,row,sftype,ncount(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical   :: Bondvillecheck

  Bondvillecheck = .false.
  do i=1,LIS_rc%nmetforc
     if (trim(LIS_rc%metforc(i)).eq."Bondville") Bondvillecheck = .true.
  enddo

  if(LIS_rc%snowsrc(n).gt.0) then
     ncount = 0 
     LIS_snow_struc(n)%sneqv = 0 
     LIS_snow_struc(n)%snowdepth = 0 

     LIS_snow_struc(n)%sneqv(:)     = 0.0
     LIS_snow_struc(n)%snowdepth(:) = 0.0
     
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
        LIS_snow_struc(n)%sneqv(tid)     = LIS_snow_struc(n)%sneqv(tid) + & 
             noah32_struc(n)%noah(t)%sneqv
        LIS_snow_struc(n)%snowdepth(tid) = LIS_snow_struc(n)%snowdepth(tid) + & 
                noah32_struc(n)%noah(t)%snowh
     enddo
  endif

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (LIS_rc%usegreennessmap(n).eq."none") then
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        call calc_localtime(LIS_rc%hr,LIS_domain(n)%grid(gid)%lon,     &
                            local_hour,change)
! For a true benchmark against the Noah3.2 testcase from NCAR,
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
        noah32_struc(n)%noah(t)%shdfac =                               &
                            month_d(noah32_struc(n)%shdfac_monthly,    &
                            locyr,locmo,locda)
     else
        tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
        noah32_struc(n)%noah(t)%shdfac =  LIS_gfrac(n)%greenness(tid)
     endif

     if (noah32_struc(n)%noah(t)%vegt.eq.LIS_rc%bareclass)             &
              noah32_struc(n)%noah(t)%shdfac = 0.0
     if (noah32_struc(n)%noah(t)%vegt.eq.LIS_rc%snowclass)             &
              noah32_struc(n)%noah(t)%shdfac = 0.0
     if (noah32_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass)            &
              noah32_struc(n)%noah(t)%shdfac = 0.0

     noah32_struc(n)%noah(t)%z0brd_old = noah32_struc(n)%noah(t)%z0brd

     if(noah32_struc(n)%noah(t)%shdfac.ge.noah32_struc(n)%noah(t)%shdmax) then
        if(LIS_rc%useemissmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
        else
           if(noah32_struc(n)%embrd_upd.eq.0) then
              noah32_struc(n)%noah(t)%embrd =                             &
                   noah32_struc(n)%emissmax(noah32_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if(noah32_struc(n)%z0brd_upd.eq.0) then
              noah32_struc(n)%noah(t)%z0brd =                             &
                   noah32_struc(n)%z0max(noah32_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
           if(noah32_struc(n)%alb_upd.eq.0) then
              noah32_struc(n)%noah(t)%alb =                            &
                    noah32_struc(n)%albmin(noah32_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
           if(noah32_struc(n)%lai_upd.eq.0) then
              noah32_struc(n)%noah(t)%lai =                            &
                    noah32_struc(n)%laimax(noah32_struc(n)%noah(t)%vegt)          
           else
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%lai = &
                   LIS_lai(n)%tlai(tid)
           endif
        endif

     elseif(noah32_struc(n)%noah(t)%shdfac.le.noah32_struc(n)%noah(t)%shdmin) then
        if(LIS_rc%useemissmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
        else
           if(noah32_struc(n)%embrd_upd.eq.0) then
              noah32_struc(n)%noah(t)%embrd =                             &
                   noah32_struc(n)%emissmin(noah32_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
           noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if(noah32_struc(n)%z0brd_upd.eq.0) then
              noah32_struc(n)%noah(t)%z0brd =                             &
                   noah32_struc(n)%z0min(noah32_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
           if(noah32_struc(n)%alb_upd.eq.0) then
              noah32_struc(n)%noah(t)%alb =                            &
                    noah32_struc(n)%albmax(noah32_struc(n)%noah(t)%vegt)
           endif
        endif
        if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
           if(noah32_struc(n)%lai_upd.eq.0) then
              noah32_struc(n)%noah(t)%lai =                            &
                    noah32_struc(n)%laimin(noah32_struc(n)%noah(t)%vegt)
           else
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%lai = &
                   LIS_lai(n)%tlai(tid)
           endif
        endif

     else
        if(noah32_struc(n)%noah(t)%shdmax.gt.noah32_struc(n)%noah(t)%shdmin) then
           interp_fraction =                                           &
                (noah32_struc(n)%noah(t)%shdfac-noah32_struc(n)%noah(t)%shdmin) / &
                (noah32_struc(n)%noah(t)%shdmax-noah32_struc(n)%noah(t)%shdmin)
           interp_fraction = min(interp_fraction,1.0)
           interp_fraction = max(interp_fraction,0.0)

           if(LIS_rc%useemissmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              if(noah32_struc(n)%embrd_upd.eq.0) then
                 noah32_struc(n)%noah(t)%embrd = ((1.0 - interp_fraction) * &
                      noah32_struc(n)%emissmin(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%emissmax(noah32_struc(n)%noah(t)%vegt))
              endif
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              if(noah32_struc(n)%z0brd_upd.eq.0) then
                 noah32_struc(n)%noah(t)%z0brd = ((1.0 - interp_fraction) * &
                      noah32_struc(n)%z0min(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%z0max(noah32_struc(n)%noah(t)%vegt))
              endif
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              if(noah32_struc(n)%lai_upd.eq.0) then
                 noah32_struc(n)%noah(t)%lai = ((1.0-interp_fraction) * &
                      noah32_struc(n)%laimin(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%laimax(noah32_struc(n)%noah(t)%vegt))
              else
                 tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%lai = &
                      LIS_lai(n)%tlai(tid)
              endif
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              if(noah32_struc(n)%alb_upd.eq.0) then
                 noah32_struc(n)%noah(t)%alb = ((1.0-interp_fraction) * &
                      noah32_struc(n)%albmax(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%albmin(noah32_struc(n)%noah(t)%vegt))
              endif
           endif
        else
           if(LIS_rc%useemissmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              if(noah32_struc(n)%embrd_upd.eq.0) then
                 noah32_struc(n)%noah(t)%embrd =                          &
                      0.5 * noah32_struc(n)%emissmin(noah32_struc(n)%noah(t)%vegt) + &
                      0.5 * noah32_struc(n)%emissmax(noah32_struc(n)%noah(t)%vegt)
              endif
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then
              tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id 
              noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              if(noah32_struc(n)%z0brd_upd.eq.0) then
                 noah32_struc(n)%noah(t)%z0brd =                          &
                      0.5 * noah32_struc(n)%z0min(noah32_struc(n)%noah(t)%vegt) + &
                      0.5 * noah32_struc(n)%z0max(noah32_struc(n)%noah(t)%vegt)
              endif
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              if(noah32_struc(n)%lai_upd.eq.0) then
                 noah32_struc(n)%noah(t)%lai =                         &
                      0.5 * noah32_struc(n)%laimin(noah32_struc(n)%noah(t)%vegt) + &
                      0.5 * noah32_struc(n)%laimax(noah32_struc(n)%noah(t)%vegt)
              else
                 tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%lai = &
                      LIS_lai(n)%tlai(tid)
              endif
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah32_struc(n)%noah(t)%alb =                            &
                   0.5 * noah32_struc(n)%albmin(noah32_struc(n)%noah(t)%vegt) + &
                   0.5 * noah32_struc(n)%albmax(noah32_struc(n)%noah(t)%vegt)
           endif
        endif
     endif

     if (LIS_rc%usealbedomap(n).ne."none") then
        if(noah32_struc(n)%alb_upd.eq.0) then
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
           noah32_struc(n)%noah(t)%alb = LIS_alb(n)%albsf(tid)
        endif
     endif
     if (LIS_rc%uselaimap(n).ne."none") then
        if(noah32_struc(n)%lai_upd.eq.0) then 
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
           noah32_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
        endif
     endif
  enddo

  return
end subroutine noah32_dynsetup
