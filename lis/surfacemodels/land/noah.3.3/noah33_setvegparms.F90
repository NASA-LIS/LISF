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
! !ROUTINE: noah33_setvegparms
! \label{noah33_setvegparms}
!
! !REVISION HISTORY:
!  07 May 2009: Sujay Kumar; Initial Implementation
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!  20 Jan 2011: David Mocko, added max/min greenness
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  14 Jan 2014: David Mocko, reconfirmed Noah3.3 in LIS7.0
!  24 Oct 2014: David Mocko, ensured max albedo maps used
!
! !INTERFACE:
subroutine noah33_setvegparms(mtype)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_vegDataMod,    only : LIS_gfrac, LIS_read_shdmax, &
       LIS_read_shdmin, LIS_lai, LIS_roughness
  use LIS_emissMod,      only : LIS_emiss
  use LIS_albedoMod,     only : LIS_alb
  use LIS_logMod,  only : LIS_logunit, LIS_endrun, &
                          LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use module_sf_noah33lsm, only : month_d, calc_localtime
  use noah33_lsmMod      

! !DESCRIPTION:
!  This subroutine retrieves Noah3.3 vegetation parameters.
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
  logical       :: file_exists

  Bondvillecheck = .false.
  do i=1,LIS_rc%nmetforc
     if (trim(LIS_rc%metforc(i)).eq."Bondville") Bondvillecheck = .true.
  enddo

  do n=1,LIS_rc%nnest

! Set the section of the vegetation parameter table to read for Noah3.3
     if (LIS_rc%lcscheme.eq."UMD")      mminlu = 'UMD '
     if (LIS_rc%lcscheme.eq."USGS")     mminlu = 'USGS'
     if (LIS_rc%lcscheme.eq."IGBPNCEP") mminlu = 'MODI'
     if (LIS_rc%lcscheme.eq."MODIS")    mminlu = 'MODI'
     if (LIS_rc%lcscheme.eq."ECOCLIMAP2") mminlu = 'ECM2'

!-----------------------------------------------------------------------
! Set Noah3.3 vegetation type at tile from the LIS domain
!-----------------------------------------------------------------------
     if (noah33_struc(n)%fixedvegtype.ne.0) then
        write(LIS_logunit,*) '[INFO] Fixing Noah3.3 vegetation to type: ',    &
                             noah33_struc(n)%fixedvegtype
     endif
     do t=1,LIS_rc%npatch(n,mtype)
        if (noah33_struc(n)%fixedvegtype.eq.0) then
           noah33_struc(n)%noah(t)%vegt = LIS_surface(n,mtype)%tile(t)%vegt
        else
           noah33_struc(n)%noah(t)%vegt = noah33_struc(n)%fixedvegtype
        endif
     enddo
     
!-----------------------------------------------------------------------
! Get Vegetation Parameters for Noah3.3 Model in Tile Space
! Read in the Noah3.3 Static Vegetation Parameter Files
!-----------------------------------------------------------------------
     inquire(file=noah33_struc(n)%vfile,exist=file_exists)
     if(.not.file_exists) then
        write(LIS_logunit,*) '[ERR] Noah3.3 veg parameter table file, ',  &
             trim(noah33_struc(n)%vfile),', does not exist. '
        write(LIS_logunit,*) '[ERR] Program stopping ...'
       call LIS_endrun()
     endif

     write(LIS_logunit,*) '[INFO] Reading Noah3.3 vegetation parameter file: ',&
          trim(noah33_struc(n)%vfile)
     ftn = LIS_getNextUnitNumber()
     open(unit=ftn,file=noah33_struc(n)%vfile,status='old')
     lumatch = 0
     do while (lumatch.eq.0)
        read(ftn,*,end=2002) lutype
        if (lutype.eq.mminlu) then
           read(ftn,*) noah33_struc(n)%lucats,iindex
           write(LIS_logunit,*) '[INFO] Noah3.3 Landuse type ',mminlu,        &
                                ' - found ',noah33_struc(n)%lucats,    &
                                ' categories'
           lumatch = 1
        endif
     enddo

 2002   continue

     allocate(noah33_struc(n)%emissmin(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%emissmax(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%laimin(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%laimax(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%albmin(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%albmax(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%z0min(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%z0max(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%shdtbl(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%nrotbl(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%rstbl(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%rgltbl(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%hstbl(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%snuptbl(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%maxalb(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%cziltbl(noah33_struc(n)%lucats))
     allocate(noah33_struc(n)%cmcmaxtbl(noah33_struc(n)%lucats))

     do j=1,NOAH33_STRUC(N)%LUCATS
        read(ftn,*) IINDEX, noah33_struc(n)%shdtbl(j), &
             noah33_struc(n)%nrotbl(j), noah33_struc(n)%rstbl(j), &
             noah33_struc(n)%rgltbl(j), &
             noah33_struc(n)%hstbl(j), noah33_struc(n)%snuptbl(j),    &
             noah33_struc(n)%maxalb(j), &
             noah33_struc(n)%laimin(j),  noah33_struc(n)%laimax(j),   &
             noah33_struc(n)%emissmin(j),noah33_struc(n)%emissmax(j), &
             noah33_struc(n)%albmin(j),  noah33_struc(n)%albmax(j),   &
             noah33_struc(n)%z0min(j),   noah33_struc(n)%z0max(j),    &
             noah33_struc(n)%cziltbl(j), noah33_struc(n)%cmcmaxtbl(j)

!        read(ftn,*) IINDEX, albtbl(j), z0tbl(j), shdtbl(j), &
!             nrotbl(j), rstbl(j), rgltbl(j), hstbl(j), &
!             snuptbl(j), laitbl(j), maxalb(j)
     enddo
!    optimum transpiration air temperature (K)
     READ (ftn,*)
     READ (ftn,*)TOPT_DATA
     READ (ftn,*)
!    maximum canopy water capacity
!     READ (ftn,*)CMCMAX_DATA
!     READ (ftn,*)
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
        write(LIS_logunit,*) '[INFO] Reading Noah3.3 greenness fraction ',    &
                             'from lis.config file'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah3.3 greenness fraction ',    &
                             'from greenness maps'
        allocate(placeshdmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        allocate(placeshdmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        call LIS_read_shdmax(n,placeshdmax)
        call LIS_read_shdmin(n,placeshdmin)
        do t=1,LIS_rc%npatch(n,mtype)
           if (placeshdmax(LIS_surface(n,mtype)%tile(t)%col,                  &
                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
              noah33_struc(n)%noah(t)%shdmax =                         &
                            placeshdmax(LIS_surface(n,mtype)%tile(t)%col,     &
                                        LIS_surface(n,mtype)%tile(t)%row)
           else
              noah33_struc(n)%noah(t)%shdmax = 1.0
              write(LIS_logunit,*) '[INFO] Noah3.3 maximum greenness -- ',    &
                                   'not defined for point ',t
              write(LIS_logunit,*) '[INFO] Noah3.3 maximum greenness -- ',    &
                                   'set to default value of 1.0'
           endif
           if (placeshdmin(LIS_surface(n,mtype)%tile(t)%col,                  &
                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
              noah33_struc(n)%noah(t)%shdmin =                         &
                            placeshdmin(LIS_surface(n,mtype)%tile(t)%col,     &
                                        LIS_surface(n,mtype)%tile(t)%row)
           else
              noah33_struc(n)%noah(t)%shdmin = 0.0
              write(LIS_logunit,*) '[INFO] Noah3.3 minimum greenness -- ',    &
                                   'not defined for point ',t
              write(LIS_logunit,*) '[INFO] Noah3.3 minimum greenness -- ',    &
                                   'set to default value of 0.0'
           endif
        enddo
        deallocate(placeshdmax)
        deallocate(placeshdmin)
     endif
     if (LIS_rc%usealbedomap(n).eq."none") then
        write(LIS_logunit,*) '[INFO] Calculating Noah3.3 albedo from ',       &
                             'greenness and albedo min/max in table'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah3.3 albedo from ',           &
                             'albedo maps'
     endif
     if (LIS_rc%uselaimap(n).eq."none") then
        write(LIS_logunit,*) '[INFO] Calculating Noah3.3 LAI from ',          &
                             'greenness and LAI min/max in table'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah3.3 LAI from ',              &
                             'LAI maps'
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        if (LIS_rc%startcode.eq."coldstart") then
           noah33_struc(n)%noah(t)%albedo = &
                month_d(noah33_struc(n)%albedo_monthly, &
                LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%z0 = LIS_roughness(n)%roughness(tid)
           else
              noah33_struc(n)%noah(t)%z0 = &
                   month_d(noah33_struc(n)%z0brd_monthly,      &
                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           endif

           noah33_struc(n)%noah(t)%z0_old = noah33_struc(n)%noah(t)%z0
        endif

        if (LIS_rc%usegreennessmap(n).eq."none") then
           gid = LIS_surface(n,mtype)%tile(t)%index
           call calc_localtime(LIS_rc%hr,LIS_domain(n)%grid(gid)%lon,  &
                               local_hour,change)
! For a true benchmark against the Noah3.3 testcase from NCAR,
! set "change = 0".  This line allows the code to _incorrectly_
! run in the same way as the Bondville testcase, which runs on
! local time instead of on UTC time. - dmm
           if (Bondvillecheck) then
              change = 0
              if (t.eq.1) &
              write(LIS_logunit,*) '[INFO] Performing benchmark of Noah3.3 ',&
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
           noah33_struc(n)%noah(t)%shdfac =                            &
                               month_d(noah33_struc(n)%shdfac_monthly, &
                               locyr,locmo,locda)
           noah33_struc(n)%noah(t)%shdmax = 0.0
           noah33_struc(n)%noah(t)%shdmin = 1.0
           do m = 1,12
              noah33_struc(n)%noah(t)%shdmax =                         &
                                  max(noah33_struc(n)%noah(t)%shdmax,  &
                                      noah33_struc(n)%shdfac_monthly(m))
              noah33_struc(n)%noah(t)%shdmin =                         &
                                  min(noah33_struc(n)%noah(t)%shdmin,  &
                                      noah33_struc(n)%shdfac_monthly(m))
           enddo
        else
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%shdfac =  LIS_gfrac(n)%greenness(tid)
        endif

        if (noah33_struc(n)%noah(t)%vegt.eq.LIS_rc%bareclass)          &
                 noah33_struc(n)%noah(t)%shdfac = 0.0
        if (noah33_struc(n)%noah(t)%vegt.eq.LIS_rc%snowclass)          &
                 noah33_struc(n)%noah(t)%shdfac = 0.0
        if (noah33_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass)         &
                 noah33_struc(n)%noah(t)%shdfac = 0.0

        if(noah33_struc(n)%noah(t)%shdfac.ge.noah33_struc(n)%noah(t)%shdmax) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah33_struc(n)%noah(t)%embrd =                             &
                   noah33_struc(n)%emissmax(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah33_struc(n)%noah(t)%z0brd =                             &
                   noah33_struc(n)%z0max(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%alb =                            &
                   noah33_struc(n)%albmin(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%lai =                            &
                   noah33_struc(n)%laimax(noah33_struc(n)%noah(t)%vegt) 
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        elseif(noah33_struc(n)%noah(t)%shdfac.le.noah33_struc(n)%noah(t)%shdmin) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah33_struc(n)%noah(t)%embrd =                             &
                   noah33_struc(n)%emissmin(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah33_struc(n)%noah(t)%z0brd =                             &
                   noah33_struc(n)%z0min(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%alb =                            &
                    noah33_struc(n)%albmax(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%lai =                            &
                    noah33_struc(n)%laimin(noah33_struc(n)%noah(t)%vegt)
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        else
           if(noah33_struc(n)%noah(t)%shdmax.gt.noah33_struc(n)%noah(t)%shdmin) then
              interp_fraction =                                        &
                   (noah33_struc(n)%noah(t)%shdfac-noah33_struc(n)%noah(t)%shdmin) / &
                   (noah33_struc(n)%noah(t)%shdmax-noah33_struc(n)%noah(t)%shdmin)
              interp_fraction = min(interp_fraction,1.0)
              interp_fraction = max(interp_fraction,0.0)

              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah33_struc(n)%noah(t)%embrd = ((1.0 - interp_fraction) * &
                      noah33_struc(n)%emissmin(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%emissmax(noah33_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah33_struc(n)%noah(t)%z0brd = ((1.0 - interp_fraction) * &
                      noah33_struc(n)%z0min(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%z0max(noah33_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%lai = ((1.0-interp_fraction) * &
                      noah33_struc(n)%laimin(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%laimax(noah33_struc(n)%noah(t)%vegt))
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%alb = ((1.0-interp_fraction) * &
                      noah33_struc(n)%albmax(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%albmin(noah33_struc(n)%noah(t)%vegt))
              endif
           else
              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah33_struc(n)%noah(t)%embrd =                          &
                      0.5 * noah33_struc(n)%emissmin(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%emissmax(noah33_struc(n)%noah(t)%vegt)
              endif

              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah33_struc(n)%noah(t)%z0brd =                          &
                      0.5 * noah33_struc(n)%z0min(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%z0max(noah33_struc(n)%noah(t)%vegt)
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%lai =                         &
                      0.5 * noah33_struc(n)%laimin(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%laimax(noah33_struc(n)%noah(t)%vegt)
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%alb =                         &
                      0.5 * noah33_struc(n)%albmin(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%albmax(noah33_struc(n)%noah(t)%vegt)
              endif
           endif
        endif

        if ((LIS_rc%usealbedomap(n).ne."none").and.(LIS_rc%startcode.eq."coldstart")) then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%alb = LIS_alb(n)%albsf(tid)
           noah33_struc(n)%noah(t)%albedo = noah33_struc(n)%noah(t)%alb
        endif
        if (LIS_rc%uselaimap(n).ne."none") then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
        endif

        noah33_struc(n)%noah(t)%topt = topt_data
!        noah33_struc(n)%noah(t)%cmcmax = cmcmax_data
        noah33_struc(n)%noah(t)%rsmax = rsmax_data
        noah33_struc(n)%noah(t)%cfactr = cfactr_data
        
     enddo

     if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
        write(LIS_logunit,*) '[INFO] Fixing Noah3.3 max snow albedo: ',       &
                              noah33_struc(n)%fixedmxsnalb
        if (noah33_struc(n)%fixedmxsnalb.eq.0.0) then
           write(LIS_logunit,*) '[WARN] Noah3.3 max snow albedo ',   &
                                'set to zero!'
        endif
     else
        if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
           write(LIS_logunit,*) '[INFO] Reading Noah3.3 max snow albedo ',    &
                                'from max snow albedo maps'
        else
           write(LIS_logunit,*) '[INFO] Reading Noah3.3 max snow albedo ',    &
                                'from vegetation parameter file'
        endif
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(t)%nroot = &
             noah33_struc(n)%nrotbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%rsmin = &
             noah33_struc(n)%rstbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%rgl = &
             noah33_struc(n)%rgltbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%hs = &
             noah33_struc(n)%hstbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%snup = &
             noah33_struc(n)%snuptbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%czil = &
             noah33_struc(n)%cziltbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%cmcmax = &
             noah33_struc(n)%cmcmaxtbl(noah33_struc(n)%noah(t)%vegt)

        if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
           noah33_struc(n)%noah(t)%mxsnalb = noah33_struc(n)%fixedmxsnalb
        else
           if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%mxsnalb = LIS_alb(n)%mxsnalb(tid)
           else
              noah33_struc(n)%noah(t)%mxsnalb = noah33_struc(n)%maxalb(noah33_struc(n)%noah(t)%vegt)/100.0
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then 
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if (LIS_rc%startcode.eq."coldstart") then
              noah33_struc(n)%noah(t)%z0brd = noah33_struc(n)%noah(t)%z0
           endif
        endif
     enddo         
  enddo

end subroutine noah33_setvegparms

!BOP
!
! !ROUTINE: noah33_resetvegparms
! \label{noah33_resetvegparms}
!
! !REVISION HISTORY:
!  07 May 2009: Sujay Kumar; Initial Implementation
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!  20 Jan 2011: David Mocko, added max/min greenness
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!
! !INTERFACE:
subroutine noah33_resetvegparms(mtype)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_vegDataMod,    only : LIS_gfrac, LIS_read_shdmax, &
       LIS_read_shdmin, LIS_lai, LIS_roughness
  use LIS_emissMod,      only : LIS_emiss
  use LIS_albedoMod,     only : LIS_alb
  use LIS_logMod,  only : LIS_logunit
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use module_sf_noah33lsm, only : month_d, calc_localtime
  use noah33_lsmMod      

! !DESCRIPTION:
!  This subroutine resets Noah3.3 vegetation parameters.
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
!!!!!            'Noah3.3 resetting vegparms'
!!!!!
!!!!!
!!!!!! Reset the section of the vegetation parameter table to read for Noah3.3
!!!!!     if (LIS_rc%lcscheme.eq."AVHRR") mminlu = 'UMD '
!!!!!     if (LIS_rc%lcscheme.eq."USGS") mminlu = 'USGS'
!!!!!     if (LIS_rc%lcscheme.eq."MODIS") mminlu = 'MODI'
!!!!!     if (LIS_rc%lcscheme.eq."ECOCLIMAP2") mminlu = 'ECM2'
!!!!!!-----------------------------------------------------------------------
!!!!!! Reset Noah3.3 vegetation type at tile from the LIS domain
!!!!!!-----------------------------------------------------------------------
!!!!!     if (noah33_struc(n)%fixedvegtype.ne.0) then
!!!!!        write(LIS_logunit,*) 'Fixing Noah3.3 vegetation to type: ',    &
!!!!!                             noah33_struc(n)%fixedvegtype
!!!!!     endif
!!!!!     do t=1,LIS_rc%npatch(n,mtype)
!!!!!        if (noah33_struc(n)%fixedvegtype.eq.0) then
!!!!!           noah33_struc(n)%noah(t)%vegt = LIS_surface(n,mtype)%tile(t)%vegt
!!!!!        else
!!!!!           noah33_struc(n)%noah(t)%vegt = noah33_struc(n)%fixedvegtype
!!!!!        endif
!!!!!     enddo
!!!!!     
!!!!!!-----------------------------------------------------------------------
!!!!!! Get Vegetation Parameters for Noah3.3 Model in Tile Space
!!!!!! Read in the Noah3.3 Static Vegetation Parameter Files
!!!!!!-----------------------------------------------------------------------
!!!!!     write(LIS_logunit,*) 'Reading Noah3.3 vegetation parameter file: ',&
!!!!!          noah33_struc(n)%vfile
!!!!!     open(unit=11,file=noah33_struc(n)%vfile,status='old')
!!!!!     lumatch = 0
!!!!!     do while (lumatch.eq.0)
!!!!!        read(11,*,end=2002) lutype
!!!!!        if (lutype.eq.mminlu) then
!!!!!           read(11,*) noah33_struc(n)%lucats,iindex
!!!!!           write(LIS_logunit,*) 'Noah3.3 Landuse type ',mminlu,        &
!!!!!                                ' - found ',noah33_struc(n)%lucats,    &
!!!!!                                ' categories'
!!!!!           lumatch = 1
!!!!!        endif
!!!!!     enddo
!!!!!
!!!!! 2002   continue
!!!!!
!!!!!!!!!     allocate(noah33_struc(n)%emissmin(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%emissmax(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%laimin(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%laimax(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%albmin(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%albmax(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%z0min(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%z0max(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%shdtbl(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%nrotbl(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%rstbl(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%rgltbl(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%hstbl(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%snuptbl(noah33_struc(n)%lucats))
!!!!!!!!!     allocate(noah33_struc(n)%maxalb(noah33_struc(n)%lucats))
!!!!!!!!!
!!!!!     do j=1,NOAH33_STRUC(N)%LUCATS
!!!!!        read(11,*) IINDEX, noah33_struc(n)%shdtbl(j), &
!!!!!             noah33_struc(n)%nrotbl(j), noah33_struc(n)%rstbl(j),&
!!!!!             noah33_struc(n)%rgltbl(j), &
!!!!!             noah33_struc(n)%hstbl(j), noah33_struc(n)%snuptbl(j), &
!!!!!             noah33_struc(n)%maxalb(j), &
!!!!!             noah33_struc(n)%laimin(j), noah33_struc(n)%laimax(j),&
!!!!!             noah33_struc(n)%emissmin(j), noah33_struc(n)%emissmax(j), &
!!!!!             noah33_struc(n)%albmin(j), noah33_struc(n)%albmax(j), &
!!!!!             noah33_struc(n)%z0min(j), noah33_struc(n)%z0max(j)
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
!!!!!        write(LIS_logunit,*) 'Reading Noah3.3 greenness fraction ',    &
!!!!!                             'from lis.config file'
!!!!!     else
!!!!!        write(LIS_logunit,*) 'Reading Noah3.3 greenness fraction ',    &
!!!!!                             'from greenness maps'
!!!!!        allocate(placeshdmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!!!        allocate(placeshdmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!!!        call LIS_read_shdmax(n,placeshdmax)
!!!!!        call LIS_read_shdmin(n,placeshdmin)
!!!!!        do t=1,LIS_rc%npatch(n,mtype)
!!!!!           if (placeshdmax(LIS_surface(n,mtype)%tile(t)%col,                  &
!!!!!                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
!!!!!              noah33_struc(n)%noah(t)%shdmax =                         &
!!!!!                            placeshdmax(LIS_surface(n,mtype)%tile(t)%col,     &
!!!!!                                        LIS_surface(n,mtype)%tile(t)%row)
!!!!!           else
!!!!!              noah33_struc(n)%noah(t)%shdmax = 1.0
!!!!!              write(LIS_logunit,*) 'Noah3.3 maximum greenness -- ',    &
!!!!!                                   'not defined for point ',t
!!!!!              write(LIS_logunit,*) 'Noah3.3 maximum greenness -- ',    &
!!!!!                                   'reset to default value of 1.0'
!!!!!           endif
!!!!!           if (placeshdmin(LIS_surface(n,mtype)%tile(t)%col,                  &
!!!!!                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
!!!!!              noah33_struc(n)%noah(t)%shdmin =                         &
!!!!!                            placeshdmin(LIS_surface(n,mtype)%tile(t)%col,     &
!!!!!                                        LIS_surface(n,mtype)%tile(t)%row)
!!!!!           else
!!!!!              noah33_struc(n)%noah(t)%shdmin = 0.0
!!!!!              write(LIS_logunit,*) 'Noah3.3 minimum greenness -- ',    &
!!!!!                                   'not defined for point ',t
!!!!!              write(LIS_logunit,*) 'Noah3.3 minimum greenness -- ',    &
!!!!!                                   'reset to default value of 0.0'
!!!!!           endif
!!!!!        enddo
!!!!!        deallocate(placeshdmax)
!!!!!        deallocate(placeshdmin)
!!!!!     endif
  if (LIS_rc%usealbedomap(n).eq."none") then
        write(LIS_logunit,*) '[INFO] Calculating Noah3.3 albedo from ',       &
                             'greenness and albedo min/max in table'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah3.3 albedo from ',           &
                             'albedo maps'
     endif
     if (LIS_rc%uselaimap(n).eq."none") then
        write(LIS_logunit,*) '[INFO] Calculating Noah3.3 LAI from ',          &
                             'greenness and LAI min/max in table'
     else
        write(LIS_logunit,*) '[INFO] Reading Noah3.3 LAI from ',              &
                             'LAI maps'
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        if (LIS_rc%startcode.eq."coldstart") then
           noah33_struc(n)%noah(t)%albedo = month_d(noah33_struc(n)%albedo_monthly, &
                                                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%z0 = LIS_roughness(n)%roughness(tid)
           else
              noah33_struc(n)%noah(t)%z0 = &
                   month_d(noah33_struc(n)%z0brd_monthly,      &
                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           endif

           noah33_struc(n)%noah(t)%z0_old = noah33_struc(n)%noah(t)%z0
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
           noah33_struc(n)%noah(t)%shdfac =                            &
                month_d(noah33_struc(n)%shdfac_monthly, &
                locyr,locmo,locda)
           noah33_struc(n)%noah(t)%shdmax = 0.0
           noah33_struc(n)%noah(t)%shdmin = 1.0
           do m = 1,12
              noah33_struc(n)%noah(t)%shdmax =                         &
                   max(noah33_struc(n)%noah(t)%shdmax,  &
                   noah33_struc(n)%shdfac_monthly(m))
              noah33_struc(n)%noah(t)%shdmin =                         &
                   min(noah33_struc(n)%noah(t)%shdmin,  &
                   noah33_struc(n)%shdfac_monthly(m))
           enddo
        else
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%shdfac =  LIS_gfrac(n)%greenness(tid)
        endif

        if (noah33_struc(n)%noah(t)%vegt.eq.LIS_rc%bareclass)          &
                 noah33_struc(n)%noah(t)%shdfac = 0.0
        if (noah33_struc(n)%noah(t)%vegt.eq.LIS_rc%snowclass)          &
                 noah33_struc(n)%noah(t)%shdfac = 0.0
        if (noah33_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass)         &
                 noah33_struc(n)%noah(t)%shdfac = 0.0

        if(noah33_struc(n)%noah(t)%shdfac.ge.noah33_struc(n)%noah(t)%shdmax) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah33_struc(n)%noah(t)%embrd =                             &
                   noah33_struc(n)%emissmax(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah33_struc(n)%noah(t)%z0brd =                             &
                   noah33_struc(n)%z0max(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%alb =                            &
                   noah33_struc(n)%albmin(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%lai =                            &
                   noah33_struc(n)%laimax(noah33_struc(n)%noah(t)%vegt) 
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        elseif(noah33_struc(n)%noah(t)%shdfac.le.noah33_struc(n)%noah(t)%shdmin) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah33_struc(n)%noah(t)%embrd =                             &
                   noah33_struc(n)%emissmin(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah33_struc(n)%noah(t)%z0brd =                             &
                   noah33_struc(n)%z0min(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%alb =                            &
                    noah33_struc(n)%albmax(noah33_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah33_struc(n)%noah(t)%lai =                            &
                    noah33_struc(n)%laimin(noah33_struc(n)%noah(t)%vegt)
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        else
           if(noah33_struc(n)%noah(t)%shdmax.gt.noah33_struc(n)%noah(t)%shdmin) then
              interp_fraction =                                        &
                   (noah33_struc(n)%noah(t)%shdfac-noah33_struc(n)%noah(t)%shdmin) / &
                   (noah33_struc(n)%noah(t)%shdmax-noah33_struc(n)%noah(t)%shdmin)
              interp_fraction = min(interp_fraction,1.0)
              interp_fraction = max(interp_fraction,0.0)

              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah33_struc(n)%noah(t)%embrd = ((1.0 - interp_fraction) * &
                      noah33_struc(n)%emissmin(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%emissmax(noah33_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah33_struc(n)%noah(t)%z0brd = ((1.0 - interp_fraction) * &
                      noah33_struc(n)%z0min(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%z0max(noah33_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%lai = ((1.0-interp_fraction) * &
                      noah33_struc(n)%laimin(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%laimax(noah33_struc(n)%noah(t)%vegt))
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%alb = ((1.0-interp_fraction) * &
                      noah33_struc(n)%albmax(noah33_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah33_struc(n)%albmin(noah33_struc(n)%noah(t)%vegt))
              endif
           else
              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah33_struc(n)%noah(t)%embrd =                          &
                      0.5 * noah33_struc(n)%emissmin(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%emissmax(noah33_struc(n)%noah(t)%vegt)
              endif

              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah33_struc(n)%noah(t)%z0brd =                          &
                      0.5 * noah33_struc(n)%z0min(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%z0max(noah33_struc(n)%noah(t)%vegt)
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%lai =                         &
                      0.5 * noah33_struc(n)%laimin(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%laimax(noah33_struc(n)%noah(t)%vegt)
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah33_struc(n)%noah(t)%alb =                         &
                      0.5 * noah33_struc(n)%albmin(noah33_struc(n)%noah(t)%vegt) + &
                      0.5 * noah33_struc(n)%albmax(noah33_struc(n)%noah(t)%vegt)
              endif
           endif
        endif

        if ((LIS_rc%usealbedomap(n).ne."none").and.(LIS_rc%startcode.eq."coldstart")) then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%alb = LIS_alb(n)%albsf(tid)
           noah33_struc(n)%noah(t)%albedo = noah33_struc(n)%noah(t)%alb
        endif
        if (LIS_rc%uselaimap(n).ne."none") then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
        endif

!!!!!        noah33_struc(n)%noah(t)%topt = topt_data
!!!!!        noah33_struc(n)%noah(t)%cmcmax = cmcmax_data
!!!!!        noah33_struc(n)%noah(t)%rsmax = rsmax_data
!!!!!        noah33_struc(n)%noah(t)%cfactr = cfactr_data
        
     enddo

     if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
        write(LIS_logunit,*) '[INFO] Fixing Noah3.3 max snow albedo: ',       &
                              noah33_struc(n)%fixedmxsnalb
        if (noah33_struc(n)%fixedmxsnalb.eq.0.0) then
           write(LIS_logunit,*) '[WARN] Noah3.3 max snow albedo ',   &
                                'set to zero!'
        endif
     else
        if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
           write(LIS_logunit,*) '[INFO] Reading Noah3.3 max snow albedo ',    &
                                'from max snow albedo maps'
        else
           write(LIS_logunit,*) '[INFO] Reading Noah3.3 max snow albedo ',    &
                                'from vegetation parameter file'
        endif
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(t)%nroot = noah33_struc(n)%nrotbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%rsmin = noah33_struc(n)%rstbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%rgl = noah33_struc(n)%rgltbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%hs = noah33_struc(n)%hstbl(noah33_struc(n)%noah(t)%vegt)
        noah33_struc(n)%noah(t)%snup = noah33_struc(n)%snuptbl(noah33_struc(n)%noah(t)%vegt)
        if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
           noah33_struc(n)%noah(t)%mxsnalb = noah33_struc(n)%fixedmxsnalb
        else
           if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah33_struc(n)%noah(t)%mxsnalb = LIS_alb(n)%mxsnalb(tid)
           else
              noah33_struc(n)%noah(t)%mxsnalb = noah33_struc(n)%maxalb(noah33_struc(n)%noah(t)%vegt)/100.0
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then 
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah33_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if (LIS_rc%startcode.eq."coldstart") then
              noah33_struc(n)%noah(t)%z0brd = noah33_struc(n)%noah(t)%z0
           endif
        endif
     enddo         
  enddo

end subroutine noah33_resetvegparms
