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
!
! !ROUTINE: noah36_setvegparms
! \label{noah36_setvegparms}
!
! !REVISION HISTORY:
!  07 May 2009: Sujay Kumar; Initial Implementation
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!  20 Jan 2011: David Mocko, added max/min greenness
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  14 Jan 2014: David Mocko, reconfirmed Noah3.3 in LIS7.0
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah36_setvegparms(mtype)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_vegDataMod,    only : LIS_gfrac, LIS_read_shdmax, &
       LIS_read_shdmin, LIS_lai, LIS_roughness
  use LIS_emissMod,      only : LIS_emiss
  use LIS_albedoMod,     only : LIS_alb
  use LIS_logMod,  only : LIS_logunit, &
                          LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use module_sf_noah36lsm, only : month_d, calc_localtime
  use noah36_lsmMod      

! !DESCRIPTION:
!  This subroutine retrieves Noah-3.6 vegetation parameters.
!  The current implementation uses a table-based lookup based
!  on vegetation classes to initialize the following parameters
!  \begin{verbatim}
!   nroot  - number of root zones
!   rsmin  - minimum canopy resistance
!   rgl    - solar radiation term in canopy resistance function
!   hs     - vapor pressure deficit term
!   snup   - threshold snow depth
!   shdfac - greenness vegetation fraction
!   shdmax - maximum greenness vegetation fraction during year
!   shdmin - minimum greenness vegetation fraction during year
!  The following variables are also a function of shdfac:
!   albedo - background (snow-free) surface albedo
!   embrd  - background (snow-free) emissivity
!   z0     - roughness length
!   lai    - leaf area index
!  \end{verbatim}
!EOP

  implicit none
  integer       :: mtype,ftn
  integer       :: n,j,t,m,gid,tid,i
  character*4   :: lutype,mminlu
  integer       :: iindex,lumatch
  real          :: interp_fraction
  real          :: topt_data
  real          :: cmcmax_data
  real          :: cfactr_data
  real          :: rsmax_data
  real          :: bare
  real          :: natural 
  integer       :: change,local_hour,locdoy,locyr,locmo,locda,lochr,locmn,locss
  real*8        :: loctime
  real          :: locgmt
  real, allocatable :: placeshdmax(:,:), placeshdmin(:,:)
  logical       :: Bondvillecheck

  Bondvillecheck = .false.
  do i=1,LIS_rc%nmetforc
     if (trim(LIS_rc%metforc(i)).eq."Bondville") Bondvillecheck = .true.
  enddo

  do n=1,LIS_rc%nnest

! Set the section of the vegetation parameter table to read for Noah-3.6
     if (LIS_rc%lcscheme.eq."UMD")      mminlu = 'UMD '
     if (LIS_rc%lcscheme.eq."USGS")     mminlu = 'USGS'
     if (LIS_rc%lcscheme.eq."IGBPNCEP") mminlu = 'MODI'
     if (LIS_rc%lcscheme.eq."MODIS")    mminlu = 'MODI'
     if (LIS_rc%lcscheme.eq."ECOCLIMAP2") mminlu = 'ECM2'

!-----------------------------------------------------------------------
! Set Noah-3.6 vegetation type at tile from the LIS domain
!-----------------------------------------------------------------------
     if (noah36_struc(n)%fixedvegtype.ne.0) then
        write(LIS_logunit,*) '[INFO] Fixing Noah-3.6 vegetation to type: ',   &
                             noah36_struc(n)%fixedvegtype
     endif
     do t=1,LIS_rc%npatch(n,mtype)
        if (noah36_struc(n)%fixedvegtype.eq.0) then
           noah36_struc(n)%noah(t)%vegt = LIS_surface(n,mtype)%tile(t)%vegt
        else
           noah36_struc(n)%noah(t)%vegt = noah36_struc(n)%fixedvegtype
        endif
     enddo
     
!-----------------------------------------------------------------------
! Get Vegetation Parameters for Noah-3.6 Model in Tile Space
! Read in the Noah-3.6 Static Vegetation Parameter Files
!-----------------------------------------------------------------------
     write(LIS_logunit,*) '[INFO] Reading Noah-3.6 vegetation parameter file: ',&
          trim(noah36_struc(n)%vfile)
     ftn = LIS_getNextUnitNumber()
     open(unit=ftn,file=noah36_struc(n)%vfile,status='old')
     lumatch = 0
     do while (lumatch.eq.0)
        read(ftn,*,end=2002) lutype
        if (lutype.eq.mminlu) then
           read(ftn,*) noah36_struc(n)%lucats,iindex
           write(LIS_logunit,*) '[INFO] Noah-3.6 Landuse type ',mminlu,       &
                                ' - found ',noah36_struc(n)%lucats,    &
                                ' categories'
           lumatch = 1
        endif
     enddo

 2002   continue

     allocate(noah36_struc(n)%emissmin(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%emissmax(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%laimin(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%laimax(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%albmin(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%albmax(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%z0min(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%z0max(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%shdtbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%nrotbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%rstbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%rgltbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%hstbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%snuptbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%maxalb(noah36_struc(n)%lucats))
!     allocate(noah36_struc(n)%cziltbl(noah36_struc(n)%lucats))
!     allocate(noah36_struc(n)%cmcmaxtbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%ztopvtbl(noah36_struc(n)%lucats))
     allocate(noah36_struc(n)%zbotvtbl(noah36_struc(n)%lucats))

     do j=1,NOAH36_STRUC(N)%LUCATS
        read(ftn,*) IINDEX, noah36_struc(n)%shdtbl(j), &
             noah36_struc(n)%nrotbl(j), noah36_struc(n)%rstbl(j), &
             noah36_struc(n)%rgltbl(j), &
             noah36_struc(n)%hstbl(j), noah36_struc(n)%snuptbl(j),    &
             noah36_struc(n)%maxalb(j), &
             noah36_struc(n)%laimin(j),  noah36_struc(n)%laimax(j),   &
             noah36_struc(n)%emissmin(j),noah36_struc(n)%emissmax(j), &
             noah36_struc(n)%albmin(j),  noah36_struc(n)%albmax(j),   &
             noah36_struc(n)%z0min(j),   noah36_struc(n)%z0max(j),    &
!             noah36_struc(n)%cziltbl(j), noah36_struc(n)%cmcmaxtbl(j)
!             noah36_struc(n)%cziltbl(j), noah36_struc(n)%cmcmaxtbl(j),&
             noah36_struc(n)%ztopvtbl(j),noah36_struc(n)%zbotvtbl(j)

!        read(ftn,*) IINDEX, albtbl(j), z0tbl(j), shdtbl(j), &
!             nrotbl(j), rstbl(j), rgltbl(j), hstbl(j), &
!             snuptbl(j), laitbl(j), maxalb(j)
     enddo
!    optimum transpiration air temperature (K)
     READ (ftn,*)
     READ (ftn,*)TOPT_DATA
     READ (ftn,*)
!    maximum canopy water capacity
     READ (ftn,*)CMCMAX_DATA
     READ (ftn,*)
!   parameter used in the canopy interception calculation
     READ (ftn,*)CFACTR_DATA
     READ (ftn,*)
!  maximum stomatal resistence 
     READ (ftn,*)RSMAX_DATA
     READ (ftn,*)
     READ (ftn,*)BARE
! This is for UCM 
     READ (ftn,*)
     READ (ftn,*)NATURAL 
     close(ftn)
     call LIS_releaseUnitNumber(ftn)
     
     if (LIS_rc%usegreennessmap(n).eq."none") then
        write(LIS_logunit,*) '[INFO] Reading Noah-3.6 greenness fraction ',   &
                             'from lis.config file'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah-3.6 greenness fraction ',   &
                             'from greenness maps'
        allocate(placeshdmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        allocate(placeshdmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        call LIS_read_shdmax(n,placeshdmax)
        call LIS_read_shdmin(n,placeshdmin)
        do t=1,LIS_rc%npatch(n,mtype)
           if (placeshdmax(LIS_surface(n,mtype)%tile(t)%col,                  &
                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
              noah36_struc(n)%noah(t)%shdmax =                         &
                            placeshdmax(LIS_surface(n,mtype)%tile(t)%col,     &
                                        LIS_surface(n,mtype)%tile(t)%row)
           else
              noah36_struc(n)%noah(t)%shdmax = 1.0
              write(LIS_logunit,*) '[WARN] Noah-3.6 maximum greenness -- ',   &
                                   'not defined for point ',t
              write(LIS_logunit,*) '[INFO] Noah-3.6 maximum greenness -- ',   &
                                   'set to default value of 1.0'
           endif
           if (placeshdmin(LIS_surface(n,mtype)%tile(t)%col,                  &
                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
              noah36_struc(n)%noah(t)%shdmin =                         &
                            placeshdmin(LIS_surface(n,mtype)%tile(t)%col,     &
                                        LIS_surface(n,mtype)%tile(t)%row)
           else
              noah36_struc(n)%noah(t)%shdmin = 0.0
              write(LIS_logunit,*) '[WARN] Noah-3.6 minimum greenness -- ',   &
                                   'not defined for point ',t
              write(LIS_logunit,*) '[WARN] Noah-3.6 minimum greenness -- ',   &
                                   'set to default value of 0.0'
           endif
        enddo
        deallocate(placeshdmax)
        deallocate(placeshdmin)
     endif
     if (LIS_rc%usealbedomap(n).eq."none") then
        write(LIS_logunit,*) '[INFO] Calculating Noah-3.6 albedo from ',      &
                             'greenness and albedo min/max in table'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah-3.6 albedo from ',          &
                             'albedo maps'
     endif
     if (LIS_rc%uselaimap(n).eq."none") then
        write(LIS_logunit,*) '[INFO] Calculating Noah-3.6 LAI from ',         &
                             'greenness and LAI min/max in table'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah-3.6 LAI from ',             &
                             'LAI maps'
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        if (LIS_rc%startcode.eq."coldstart") then
           noah36_struc(n)%noah(t)%albedo = month_d(noah36_struc(n)%albedo_monthly, &
                                                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%z0 = LIS_roughness(n)%roughness(tid)
           else
              noah36_struc(n)%noah(t)%z0 = &
                   month_d(noah36_struc(n)%z0brd_monthly,      &
                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           endif

           noah36_struc(n)%noah(t)%z0_old = noah36_struc(n)%noah(t)%z0
        endif

        if (LIS_rc%usegreennessmap(n).eq."none") then
           gid = LIS_surface(n,mtype)%tile(t)%index
           call calc_localtime(LIS_rc%hr,LIS_domain(n)%grid(gid)%lon,  &
                               local_hour,change)
! For a true benchmark against the Noah-3.6 testcase from NCAR,
! set "change = 0".  This line allows the code to _incorrectly_
! run in the same way as the Bondville testcase, which runs on
! local time instead of on UTC time. - dmm
           if (Bondvillecheck) then
              change = 0
              if (t.eq.1) &
              write(LIS_logunit,*) '[INFO] Performing benchmark of Noah-3.6 ',&
                                   'against Bondville testcase'
           endif
           change = change * 3600
           locyr = LIS_rc%yr
           locmo = LIS_rc%mo
           locda = LIS_rc%da
           lochr = LIS_rc%hr
           locmn = LIS_rc%mn
           locss = LIS_rc%ss
           call LIS_date2time(loctime,locdoy,locgmt,locyr,locmo,       &
                              locda,lochr,locmn,locss)
           call LIS_tick(loctime,locdoy,locgmt,locyr,locmo,            &
                         locda,lochr,locmn,locss,real(change))
           noah36_struc(n)%noah(t)%shdfac =                            &
                               month_d(noah36_struc(n)%shdfac_monthly, &
                               locyr,locmo,locda)
           noah36_struc(n)%noah(t)%shdmax = 0.0
           noah36_struc(n)%noah(t)%shdmin = 1.0
           do m = 1,12
              noah36_struc(n)%noah(t)%shdmax =                         &
                                  max(noah36_struc(n)%noah(t)%shdmax,  &
                                      noah36_struc(n)%shdfac_monthly(m))
              noah36_struc(n)%noah(t)%shdmin =                         &
                                  min(noah36_struc(n)%noah(t)%shdmin,  &
                                      noah36_struc(n)%shdfac_monthly(m))
           enddo
        else
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%shdfac =  LIS_gfrac(n)%greenness(tid)
        endif

        if (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%bareclass)          &
                 noah36_struc(n)%noah(t)%shdfac = 0.0
        if (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%snowclass)          &
                 noah36_struc(n)%noah(t)%shdfac = 0.0
        if (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass)         &
                 noah36_struc(n)%noah(t)%shdfac = 0.0

        if(noah36_struc(n)%noah(t)%shdfac.ge.noah36_struc(n)%noah(t)%shdmax) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah36_struc(n)%noah(t)%embrd =                             &
                   noah36_struc(n)%emissmax(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah36_struc(n)%noah(t)%z0brd =                             &
                   noah36_struc(n)%z0max(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%alb =                            &
                   noah36_struc(n)%albmin(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%lai =                            &
                   noah36_struc(n)%laimax(noah36_struc(n)%noah(t)%vegt) 
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        elseif(noah36_struc(n)%noah(t)%shdfac.le.noah36_struc(n)%noah(t)%shdmin) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah36_struc(n)%noah(t)%embrd =                             &
                   noah36_struc(n)%emissmin(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah36_struc(n)%noah(t)%z0brd =                             &
                   noah36_struc(n)%z0min(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%alb =                            &
                    noah36_struc(n)%albmax(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%lai =                            &
                    noah36_struc(n)%laimin(noah36_struc(n)%noah(t)%vegt)
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        else
           if(noah36_struc(n)%noah(t)%shdmax.gt.noah36_struc(n)%noah(t)%shdmin) then
              interp_fraction =                                        &
                   (noah36_struc(n)%noah(t)%shdfac-noah36_struc(n)%noah(t)%shdmin) / &
                   (noah36_struc(n)%noah(t)%shdmax-noah36_struc(n)%noah(t)%shdmin)
              interp_fraction = min(interp_fraction,1.0)
              interp_fraction = max(interp_fraction,0.0)

              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah36_struc(n)%noah(t)%embrd = ((1.0 - interp_fraction) * &
                      noah36_struc(n)%emissmin(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%emissmax(noah36_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah36_struc(n)%noah(t)%z0brd = ((1.0 - interp_fraction) * &
                      noah36_struc(n)%z0min(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%z0max(noah36_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%lai = ((1.0-interp_fraction) * &
                      noah36_struc(n)%laimin(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%laimax(noah36_struc(n)%noah(t)%vegt))
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%alb = ((1.0-interp_fraction) * &
                      noah36_struc(n)%albmax(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%albmin(noah36_struc(n)%noah(t)%vegt))
              endif
           else
              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah36_struc(n)%noah(t)%embrd =                          &
                      0.5 * noah36_struc(n)%emissmin(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%emissmax(noah36_struc(n)%noah(t)%vegt)
              endif

              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah36_struc(n)%noah(t)%z0brd =                          &
                      0.5 * noah36_struc(n)%z0min(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%z0max(noah36_struc(n)%noah(t)%vegt)
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%lai =                         &
                      0.5 * noah36_struc(n)%laimin(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%laimax(noah36_struc(n)%noah(t)%vegt)
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%alb =                         &
                      0.5 * noah36_struc(n)%albmin(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%albmax(noah36_struc(n)%noah(t)%vegt)
              endif
           endif
        endif

        if ((LIS_rc%usealbedomap(n).ne."none").and.(LIS_rc%startcode.eq."coldstart")) then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%alb = LIS_alb(n)%albsf(tid)
           noah36_struc(n)%noah(t)%albedo = noah36_struc(n)%noah(t)%alb
        endif
        if (LIS_rc%uselaimap(n).ne."none") then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
        endif

        noah36_struc(n)%noah(t)%topt = topt_data
        noah36_struc(n)%noah(t)%cmcmax = cmcmax_data
        noah36_struc(n)%noah(t)%rsmax = rsmax_data
        noah36_struc(n)%noah(t)%cfactr = cfactr_data
        
     enddo

     if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
        write(LIS_logunit,*) '[INFO] Fixing Noah-3.6 max snow albedo: ',      &
                              noah36_struc(n)%fixedmxsnalb
        if (noah36_struc(n)%fixedmxsnalb.eq.0.0) then
           write(LIS_logunit,*) '[WARN] Noah-3.6 max snow albedo ',  &
                                'set to zero!'
        endif
     else
        if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
           write(LIS_logunit,*) '[INFO] Reading Noah-3.6 max snow albedo ',   &
                                'from max snow albedo maps'
        else
           write(LIS_logunit,*) '[INFO] Reading Noah-3.6 max snow albedo ',   &
                                'from vegetation parameter file'
        endif
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        noah36_struc(n)%noah(t)%nroot = &
             noah36_struc(n)%nrotbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%rsmin = &
             noah36_struc(n)%rstbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%rgl = &
             noah36_struc(n)%rgltbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%hs = &
             noah36_struc(n)%hstbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%snup = &
             noah36_struc(n)%snuptbl(noah36_struc(n)%noah(t)%vegt)
!        noah36_struc(n)%noah(t)%czil = &
!             noah36_struc(n)%cziltbl(noah36_struc(n)%noah(t)%vegt)
!        noah36_struc(n)%noah(t)%cmcmax = &
!             noah36_struc(n)%cmcmaxtbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%ztopv = &
             noah36_struc(n)%ztopvtbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%zbotv = &
             noah36_struc(n)%zbotvtbl(noah36_struc(n)%noah(t)%vegt)

        if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
           noah36_struc(n)%noah(t)%mxsnalb = noah36_struc(n)%fixedmxsnalb
        else
           if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%mxsnalb = LIS_alb(n)%mxsnalb(tid)
           else
              noah36_struc(n)%noah(t)%mxsnalb = noah36_struc(n)%maxalb(noah36_struc(n)%noah(t)%vegt)/100.0
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then 
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if (LIS_rc%startcode.eq."coldstart") then
              noah36_struc(n)%noah(t)%z0brd = noah36_struc(n)%noah(t)%z0
           endif
        endif
     enddo

! If not coupled to WRF, set any open water points to (ice = 1),
! to prevent the Noah LSM from running in noah36_main.F90 - dmm
#if (!defined COUPLED)
     do t=1,LIS_rc%npatch(n,mtype)
        if (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%glacierclass) then
           noah36_struc(n)%noah(t)%xice = -1
           noah36_struc(n)%noah(t)%alb =                               &
                    noah36_struc(n)%albmax(noah36_struc(n)%noah(t)%vegt)
           noah36_struc(n)%noah(t)%z0brd =                             &
                     noah36_struc(n)%z0min(noah36_struc(n)%noah(t)%vegt)
        elseif (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass) then
           noah36_struc(n)%noah(t)%xice = 1
        else
           noah36_struc(n)%noah(t)%xice = 0
        endif
     enddo
#endif
  enddo

end subroutine noah36_setvegparms

!BOP
!
! !ROUTINE: noah36_resetvegparms
! \label{noah36_resetvegparms}
!
! !REVISION HISTORY:
!  07 May 2009: Sujay Kumar; Initial Implementation
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!  20 Jan 2011: David Mocko, added max/min greenness
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah36_resetvegparms(mtype)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_vegDataMod,    only : LIS_gfrac, LIS_read_shdmax, &
       LIS_read_shdmin, LIS_lai, LIS_roughness
  use LIS_emissMod,      only : LIS_emiss
  use LIS_albedoMod,     only : LIS_alb
  use LIS_logMod,  only : LIS_logunit
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use module_sf_noah36lsm, only : month_d, calc_localtime
  use noah36_lsmMod      

! !DESCRIPTION:
!  This subroutine resets Noah-3.6 vegetation parameters.
!  The current implementation uses a table-based lookup based
!  on vegetation classes to initialize the following parameters
!  \begin{verbatim}
!   nroot  - number of root zones
!   rsmin  - minimum canopy resistance
!   rgl    - solar radiation term in canopy resistance function
!   hs     - vapor pressure deficit term
!   snup   - threshold snow depth
!   shdfac - greenness vegetation fraction
!   shdmax - maximum greenness vegetation fraction during year
!   shdmin - minimum greenness vegetation fraction during year
!  The following variables are also a function of shdfac:
!   albedo - background (snow-free) surface albedo
!   embrd  - background (snow-free) emissivity
!   z0     - roughness length
!   lai    - leaf area index
!  \end{verbatim}
!EOP

  implicit none
  integer :: mtype
  integer :: n,j,t,m,gid,tid
!!!!!  character*4      :: lutype,mminlu
!!!!!  integer          :: iindex,lumatch
  real             :: interp_fraction
!!!!!  real             :: topt_data
!!!!!  real             :: cmcmax_data
!!!!!  real             :: cfactr_data
!!!!!  real             :: rsmax_data
!!!!!  real             :: bare
!!!!!  real             :: natural 
  integer   :: change,local_hour,locdoy,locyr,locmo,locda,lochr,locmn,locss
  real*8    :: loctime
  real      :: locgmt
!!!!!  real, allocatable :: placeshdmax(:,:), placeshdmin(:,:)

  do n=1,LIS_rc%nnest
!!!!!
!!!!!       write(LIS_logunit,*)                        &
!!!!!            'Noah-3.6 resetting vegparms'
!!!!!
!!!!!
!!!!!! Reset the section of the vegetation parameter table to read for Noah-3.6
!!!!!     if (LIS_rc%lcscheme.eq."AVHRR") mminlu = 'UMD '
!!!!!     if (LIS_rc%lcscheme.eq."USGS") mminlu = 'USGS'
!!!!!     if (LIS_rc%lcscheme.eq."MODIS") mminlu = 'MODI'
!!!!!     if (LIS_rc%lcscheme.eq."ECOCLIMAP2") mminlu = 'ECM2'
!!!!!!-----------------------------------------------------------------------
!!!!!! Reset Noah-3.6 vegetation type at tile from the LIS domain
!!!!!!-----------------------------------------------------------------------
!!!!!     if (noah36_struc(n)%fixedvegtype.ne.0) then
!!!!!        write(LIS_logunit,*) 'Fixing Noah-3.6 vegetation to type: ',    &
!!!!!                             noah36_struc(n)%fixedvegtype
!!!!!     endif
!!!!!     do t=1,LIS_rc%npatch(n,mtype)
!!!!!        if (noah36_struc(n)%fixedvegtype.eq.0) then
!!!!!           noah36_struc(n)%noah(t)%vegt = LIS_surface(n,mtype)%tile(t)%vegt
!!!!!        else
!!!!!           noah36_struc(n)%noah(t)%vegt = noah36_struc(n)%fixedvegtype
!!!!!        endif
!!!!!     enddo
!!!!!     
!!!!!!-----------------------------------------------------------------------
!!!!!! Get Vegetation Parameters for Noah-3.6 Model in Tile Space
!!!!!! Read in the Noah-3.6 Static Vegetation Parameter Files
!!!!!!-----------------------------------------------------------------------
!!!!!     write(LIS_logunit,*) 'Reading Noah-3.6 vegetation parameter file: ',&
!!!!!          noah36_struc(n)%vfile
!!!!!     open(unit=11,file=noah36_struc(n)%vfile,status='old')
!!!!!     lumatch = 0
!!!!!     do while (lumatch.eq.0)
!!!!!        read(11,*,end=2002) lutype
!!!!!        if (lutype.eq.mminlu) then
!!!!!           read(11,*) noah36_struc(n)%lucats,iindex
!!!!!           write(LIS_logunit,*) 'Noah-3.6 Landuse type ',mminlu,       &
!!!!!                                ' - found ',noah36_struc(n)%lucats,    &
!!!!!                                ' categories'
!!!!!           lumatch = 1
!!!!!        endif
!!!!!     enddo
!!!!!
!!!!! 2002   continue
!!!!!
!!!!!!!!!     allocate(noah36_struc(n)%emissmin(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%emissmax(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%laimin(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%laimax(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%albmin(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%albmax(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%z0min(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%z0max(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%shdtbl(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%nrotbl(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%rstbl(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%rgltbl(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%hstbl(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%snuptbl(noah36_struc(n)%lucats))
!!!!!!!!!     allocate(noah36_struc(n)%maxalb(noah36_struc(n)%lucats))
!!!!!!!!!
!!!!!     do j=1,NOAH36_STRUC(N)%LUCATS
!!!!!        read(11,*) IINDEX, noah36_struc(n)%shdtbl(j), &
!!!!!             noah36_struc(n)%nrotbl(j), noah36_struc(n)%rstbl(j),&
!!!!!             noah36_struc(n)%rgltbl(j), &
!!!!!             noah36_struc(n)%hstbl(j), noah36_struc(n)%snuptbl(j), &
!!!!!             noah36_struc(n)%maxalb(j), &
!!!!!             noah36_struc(n)%laimin(j), noah36_struc(n)%laimax(j),&
!!!!!             noah36_struc(n)%emissmin(j), noah36_struc(n)%emissmax(j), &
!!!!!             noah36_struc(n)%albmin(j), noah36_struc(n)%albmax(j), &
!!!!!             noah36_struc(n)%z0min(j), noah36_struc(n)%z0max(j)
!!!!!
!!!!!     enddo
!!!!!!    optimum transpiration air temperature (K)
!!!!!     READ (11,*)
!!!!!     READ (11,*)TOPT_DATA
!!!!!     READ (11,*)
!!!!!!    maximum canopy water capacity
!!!!!     READ (11,*)CMCMAX_DATA
!!!!!     READ (11,*)
!!!!!!   parameter used in the canopy interception calculation
!!!!!     READ (11,*)CFACTR_DATA
!!!!!     READ (11,*)
!!!!!!  maximum stomatal resistence 
!!!!!     READ (11,*)RSMAX_DATA
!!!!!     READ (11,*)
!!!!!     READ (11,*)BARE
!!!!!! This is for UCM 
!!!!!     READ (11,*)
!!!!!     READ (11,*)NATURAL 
!!!!!     close(11)
!!!!!  
!!!!!     if (LIS_rc%usegreennessmap(n).eq."none") then
!!!!!        write(LIS_logunit,*) 'Reading Noah-3.6 greenness fraction ',   &
!!!!!                             'from lis.config file'
!!!!!     else
!!!!!        write(LIS_logunit,*) 'Reading Noah-3.6 greenness fraction ',   &
!!!!!                             'from greenness maps'
!!!!!        allocate(placeshdmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!!!        allocate(placeshdmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!!!        call LIS_read_shdmax(n,placeshdmax)
!!!!!        call LIS_read_shdmin(n,placeshdmin)
!!!!!        do t=1,LIS_rc%npatch(n,mtype)
!!!!!           if (placeshdmax(LIS_surface(n,mtype)%tile(t)%col,                  &
!!!!!                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
!!!!!              noah36_struc(n)%noah(t)%shdmax =                         &
!!!!!                            placeshdmax(LIS_surface(n,mtype)%tile(t)%col,     &
!!!!!                                        LIS_surface(n,mtype)%tile(t)%row)
!!!!!           else
!!!!!              noah36_struc(n)%noah(t)%shdmax = 1.0
!!!!!              write(LIS_logunit,*) 'Noah-3.6 maximum greenness -- ',   &
!!!!!                                   'not defined for point ',t
!!!!!              write(LIS_logunit,*) 'Noah-3.6 maximum greenness -- ',   &
!!!!!                                   'reset to default value of 1.0'
!!!!!           endif
!!!!!           if (placeshdmin(LIS_surface(n,mtype)%tile(t)%col,                  &
!!!!!                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
!!!!!              noah36_struc(n)%noah(t)%shdmin =                         &
!!!!!                            placeshdmin(LIS_surface(n,mtype)%tile(t)%col,     &
!!!!!                                        LIS_surface(n,mtype)%tile(t)%row)
!!!!!           else
!!!!!              noah36_struc(n)%noah(t)%shdmin = 0.0
!!!!!              write(LIS_logunit,*) 'Noah-3.6 minimum greenness -- ',   &
!!!!!                                   'not defined for point ',t
!!!!!              write(LIS_logunit,*) 'Noah-3.6 minimum greenness -- ',   &
!!!!!                                   'reset to default value of 0.0'
!!!!!           endif
!!!!!        enddo
!!!!!        deallocate(placeshdmax)
!!!!!        deallocate(placeshdmin)
!!!!!     endif
  if (LIS_rc%usealbedomap(n).eq."none") then
        write(LIS_logunit,*) 'Calculating Noah-3.6 albedo from ',      &
                             'greenness and albedo min/max in table'
     else
        write(LIS_logunit,*) 'Reading Noah-3.6 albedo from ',          &
                             'albedo maps'
     endif
     if (LIS_rc%uselaimap(n).eq."none") then
        write(LIS_logunit,*) 'Calculating Noah-3.6 LAI from ',         &
                             'greenness and LAI min/max in table'
     else
        write(LIS_logunit,*) 'Reading Noah-3.6 LAI from ',             &
                             'LAI maps'
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        if (LIS_rc%startcode.eq."coldstart") then
           noah36_struc(n)%noah(t)%albedo = month_d(noah36_struc(n)%albedo_monthly, &
                                                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%z0 = LIS_roughness(n)%roughness(tid)
           else
              noah36_struc(n)%noah(t)%z0 = &
                   month_d(noah36_struc(n)%z0brd_monthly,      &
                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           endif

           noah36_struc(n)%noah(t)%z0_old = noah36_struc(n)%noah(t)%z0
        endif

        if (LIS_rc%usegreennessmap(n).eq."none") then
           gid = LIS_surface(n,mtype)%tile(t)%index
           call calc_localtime(LIS_rc%hr,LIS_domain(n)%grid(gid)%lon,  &
                               local_hour,change)
           change = change * 3600
           locyr = LIS_rc%yr
           locmo = LIS_rc%mo
           locda = LIS_rc%da
           lochr = LIS_rc%hr
           locmn = LIS_rc%mn
           locss = LIS_rc%ss
           call LIS_date2time(loctime,locdoy,locgmt,locyr,locmo,       &
                locda,lochr,locmn,locss)
           call LIS_tick(loctime,locdoy,locgmt,locyr,locmo,            &
                locda,lochr,locmn,locss,real(change))
           noah36_struc(n)%noah(t)%shdfac =                            &
                month_d(noah36_struc(n)%shdfac_monthly, &
                locyr,locmo,locda)
           noah36_struc(n)%noah(t)%shdmax = 0.0
           noah36_struc(n)%noah(t)%shdmin = 1.0
           do m = 1,12
              noah36_struc(n)%noah(t)%shdmax =                         &
                   max(noah36_struc(n)%noah(t)%shdmax,  &
                   noah36_struc(n)%shdfac_monthly(m))
              noah36_struc(n)%noah(t)%shdmin =                         &
                   min(noah36_struc(n)%noah(t)%shdmin,  &
                   noah36_struc(n)%shdfac_monthly(m))
           enddo
        else
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%shdfac =  LIS_gfrac(n)%greenness(tid)
        endif

        if (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%bareclass)          &
                 noah36_struc(n)%noah(t)%shdfac = 0.0
        if (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%snowclass)          &
                 noah36_struc(n)%noah(t)%shdfac = 0.0
        if (noah36_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass)         &
                 noah36_struc(n)%noah(t)%shdfac = 0.0

        if(noah36_struc(n)%noah(t)%shdfac.ge.noah36_struc(n)%noah(t)%shdmax) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah36_struc(n)%noah(t)%embrd =                             &
                   noah36_struc(n)%emissmax(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah36_struc(n)%noah(t)%z0brd =                             &
                   noah36_struc(n)%z0max(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%alb =                            &
                   noah36_struc(n)%albmin(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%lai =                            &
                   noah36_struc(n)%laimax(noah36_struc(n)%noah(t)%vegt) 
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        elseif(noah36_struc(n)%noah(t)%shdfac.le.noah36_struc(n)%noah(t)%shdmin) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah36_struc(n)%noah(t)%embrd =                             &
                   noah36_struc(n)%emissmin(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah36_struc(n)%noah(t)%z0brd =                             &
                   noah36_struc(n)%z0min(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%alb =                            &
                    noah36_struc(n)%albmax(noah36_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah36_struc(n)%noah(t)%lai =                            &
                    noah36_struc(n)%laimin(noah36_struc(n)%noah(t)%vegt)
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        else
           if(noah36_struc(n)%noah(t)%shdmax.gt.noah36_struc(n)%noah(t)%shdmin) then
              interp_fraction =                                        &
                   (noah36_struc(n)%noah(t)%shdfac-noah36_struc(n)%noah(t)%shdmin) / &
                   (noah36_struc(n)%noah(t)%shdmax-noah36_struc(n)%noah(t)%shdmin)
              interp_fraction = min(interp_fraction,1.0)
              interp_fraction = max(interp_fraction,0.0)

              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah36_struc(n)%noah(t)%embrd = ((1.0 - interp_fraction) * &
                      noah36_struc(n)%emissmin(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%emissmax(noah36_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah36_struc(n)%noah(t)%z0brd = ((1.0 - interp_fraction) * &
                      noah36_struc(n)%z0min(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%z0max(noah36_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%lai = ((1.0-interp_fraction) * &
                      noah36_struc(n)%laimin(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%laimax(noah36_struc(n)%noah(t)%vegt))
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%alb = ((1.0-interp_fraction) * &
                      noah36_struc(n)%albmax(noah36_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah36_struc(n)%albmin(noah36_struc(n)%noah(t)%vegt))
              endif
           else
              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah36_struc(n)%noah(t)%embrd =                          &
                      0.5 * noah36_struc(n)%emissmin(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%emissmax(noah36_struc(n)%noah(t)%vegt)
              endif

              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah36_struc(n)%noah(t)%z0brd =                          &
                      0.5 * noah36_struc(n)%z0min(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%z0max(noah36_struc(n)%noah(t)%vegt)
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%lai =                         &
                      0.5 * noah36_struc(n)%laimin(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%laimax(noah36_struc(n)%noah(t)%vegt)
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah36_struc(n)%noah(t)%alb =                         &
                      0.5 * noah36_struc(n)%albmin(noah36_struc(n)%noah(t)%vegt) + &
                      0.5 * noah36_struc(n)%albmax(noah36_struc(n)%noah(t)%vegt)
              endif
           endif
        endif

        if ((LIS_rc%usealbedomap(n).ne."none").and.(LIS_rc%startcode.eq."coldstart")) then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%alb = LIS_alb(n)%albsf(tid)
           noah36_struc(n)%noah(t)%albedo = noah36_struc(n)%noah(t)%alb
        endif
        if (LIS_rc%uselaimap(n).ne."none") then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
        endif

!!!!!        noah36_struc(n)%noah(t)%topt = topt_data
!!!!!        noah36_struc(n)%noah(t)%cmcmax = cmcmax_data
!!!!!        noah36_struc(n)%noah(t)%rsmax = rsmax_data
!!!!!        noah36_struc(n)%noah(t)%cfactr = cfactr_data
        
     enddo

     if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
        write(LIS_logunit,*) 'Fixing Noah-3.6 max snow albedo: ',      &
                              noah36_struc(n)%fixedmxsnalb
        if (noah36_struc(n)%fixedmxsnalb.eq.0.0) then
           write(LIS_logunit,*) 'WARNING: Noah-3.6 max snow albedo ',  &
                                'set to zero!'
        endif
     else
        if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
           write(LIS_logunit,*) 'Reading Noah-3.6 max snow albedo ',   &
                                'from max snow albedo maps'
        else
           write(LIS_logunit,*) 'Reading Noah-3.6 max snow albedo ',   &
                                'from vegetation parameter file'
        endif
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        noah36_struc(n)%noah(t)%nroot = noah36_struc(n)%nrotbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%rsmin = noah36_struc(n)%rstbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%rgl = noah36_struc(n)%rgltbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%hs = noah36_struc(n)%hstbl(noah36_struc(n)%noah(t)%vegt)
        noah36_struc(n)%noah(t)%snup = noah36_struc(n)%snuptbl(noah36_struc(n)%noah(t)%vegt)
        if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
           noah36_struc(n)%noah(t)%mxsnalb = noah36_struc(n)%fixedmxsnalb
        else
           if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah36_struc(n)%noah(t)%mxsnalb = LIS_alb(n)%mxsnalb(tid)
           else
              noah36_struc(n)%noah(t)%mxsnalb = noah36_struc(n)%maxalb(noah36_struc(n)%noah(t)%vegt)/100.0
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then 
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah36_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if (LIS_rc%startcode.eq."coldstart") then
              noah36_struc(n)%noah(t)%z0brd = noah36_struc(n)%noah(t)%z0
           endif
        endif
     enddo         
  enddo

end subroutine noah36_resetvegparms
