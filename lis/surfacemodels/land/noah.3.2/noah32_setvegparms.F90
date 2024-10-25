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
! !ROUTINE: noah32_setvegparms
! \label{noah32_setvegparms}
!
! !REVISION HISTORY:
!  07 May 2009: Sujay Kumar; Initial Implementation
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!  20 Jan 2011: David Mocko, added max/min greenness
!  24 Oct 2014: David Mocko, ensured max albedo maps used
!
! !INTERFACE:
subroutine noah32_setvegparms(mtype)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_vegDataMod,    only : LIS_gfrac, LIS_read_shdmax, &
       LIS_read_shdmin, LIS_lai, LIS_roughness
  use LIS_emissMod,      only : LIS_emiss
  use LIS_albedoMod,     only : LIS_alb
  use LIS_logMod,  only : LIS_logunit, &
                          LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_tick
  use module_sf_noah32lsm, only : month_d, calc_localtime
  use noah32_lsmMod      

! !DESCRIPTION:
!  This subroutine retrieves Noah3.2 vegetation parameters.
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

! Set the section of the vegetation parameter table to read for Noah3.2
     if (LIS_rc%lcscheme.eq."UMD")      mminlu = 'UMD '
     if (LIS_rc%lcscheme.eq."USGS")     mminlu = 'USGS'
     if (LIS_rc%lcscheme.eq."IGBPNCEP") mminlu = 'MODI'
     if (LIS_rc%lcscheme.eq."MODIS")    mminlu = 'MODI'
     if (LIS_rc%lcscheme.eq."ECOCLIMAP2") mminlu = 'ECM2'

!-----------------------------------------------------------------------
! Set Noah3.2 vegetation type at tile from the LIS domain
!-----------------------------------------------------------------------
     if (noah32_struc(n)%fixedvegtype.ne.0) then
        write(LIS_logunit,*) 'Fixing Noah3.2 vegetation to type: ',    &
                             noah32_struc(n)%fixedvegtype
     endif
     do t=1,LIS_rc%npatch(n,mtype)
        if (noah32_struc(n)%fixedvegtype.eq.0) then
           noah32_struc(n)%noah(t)%vegt = LIS_surface(n,mtype)%tile(t)%vegt
        else
           noah32_struc(n)%noah(t)%vegt = noah32_struc(n)%fixedvegtype
        endif
     enddo
     
!-----------------------------------------------------------------------
! Get Vegetation Parameters for Noah3.2 Model in Tile Space
! Read in the Noah3.2 Static Vegetation Parameter Files
!-----------------------------------------------------------------------
     write(LIS_logunit,*) 'Reading Noah3.2 vegetation parameter file: ',&
          trim(noah32_struc(n)%vfile)
     ftn = LIS_getNextUnitNumber()
     open(unit=ftn,file=noah32_struc(n)%vfile,status='old')
     lumatch = 0
     do while (lumatch.eq.0)
        read(ftn,*,end=2002) lutype
        if (lutype.eq.mminlu) then
           read(ftn,*) noah32_struc(n)%lucats,iindex
           write(LIS_logunit,*) 'Noah3.2 Landuse type ',mminlu,        &
                                ' - found ',noah32_struc(n)%lucats,    &
                                ' categories'
           lumatch = 1
        endif
     enddo

 2002   continue

     allocate(noah32_struc(n)%emissmin(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%emissmax(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%laimin(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%laimax(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%albmin(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%albmax(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%z0min(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%z0max(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%shdtbl(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%nrotbl(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%rstbl(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%rgltbl(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%hstbl(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%snuptbl(noah32_struc(n)%lucats))
     allocate(noah32_struc(n)%maxalb(noah32_struc(n)%lucats))

     do j=1,NOAH32_STRUC(N)%LUCATS
        read(ftn,*) IINDEX, noah32_struc(n)%shdtbl(j), &
             noah32_struc(n)%nrotbl(j), noah32_struc(n)%rstbl(j), &
             noah32_struc(n)%rgltbl(j), &
             noah32_struc(n)%hstbl(j), noah32_struc(n)%snuptbl(j), &
             noah32_struc(n)%maxalb(j), &
             noah32_struc(n)%laimin(j), noah32_struc(n)%laimax(j), &
             noah32_struc(n)%emissmin(j), noah32_struc(n)%emissmax(j), &
             noah32_struc(n)%albmin(j), noah32_struc(n)%albmax(j), &
             noah32_struc(n)%z0min(j), noah32_struc(n)%z0max(j)

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
        write(LIS_logunit,*) 'Reading Noah3.2 greenness fraction ',    &
                             'from lis.config file'
     else
        write(LIS_logunit,*) 'Reading Noah3.2 greenness fraction ',    &
                             'from greenness maps'
        allocate(placeshdmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        allocate(placeshdmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        call LIS_read_shdmax(n,placeshdmax)
        call LIS_read_shdmin(n,placeshdmin)
        do t=1,LIS_rc%npatch(n,mtype)
           if (placeshdmax(LIS_surface(n,mtype)%tile(t)%col,                  &
                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
              noah32_struc(n)%noah(t)%shdmax =                         &
                            placeshdmax(LIS_surface(n,mtype)%tile(t)%col,     &
                                        LIS_surface(n,mtype)%tile(t)%row)
           else
              noah32_struc(n)%noah(t)%shdmax = 1.0
              write(LIS_logunit,*) 'Noah3.2 maximum greenness -- ',    &
                                   'not defined for point ',t
              write(LIS_logunit,*) 'Noah3.2 maximum greenness -- ',    &
                                   'set to default value of 1.0'
           endif
           if (placeshdmin(LIS_surface(n,mtype)%tile(t)%col,                  &
                           LIS_surface(n,mtype)%tile(t)%row).ne.-9999.00) then
              noah32_struc(n)%noah(t)%shdmin =                         &
                            placeshdmin(LIS_surface(n,mtype)%tile(t)%col,     &
                                        LIS_surface(n,mtype)%tile(t)%row)
           else
              noah32_struc(n)%noah(t)%shdmin = 0.0
              write(LIS_logunit,*) 'Noah3.2 minimum greenness -- ',    &
                                   'not defined for point ',t
              write(LIS_logunit,*) 'Noah3.2 minimum greenness -- ',    &
                                   'set to default value of 0.0'
           endif
        enddo
        deallocate(placeshdmax)
        deallocate(placeshdmin)
     endif
     if (LIS_rc%usealbedomap(n).eq."none") then
        write(LIS_logunit,*) 'Calculating Noah3.2 albedo from ',       &
                             'greenness and albedo min/max in table'
     else
        write(LIS_logunit,*) 'Reading Noah3.2 albedo from ',           &
                             'albedo maps'
     endif
     if (LIS_rc%uselaimap(n).eq."none") then
        write(LIS_logunit,*) 'Calculating Noah3.2 LAI from ',          &
                             'greenness and LAI min/max in table'
     else
        write(LIS_logunit,*) 'Reading Noah3.2 LAI from ',              &
                             'LAI maps'
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        if (LIS_rc%startcode.eq."coldstart") then
           noah32_struc(n)%noah(t)%albedo = month_d(noah32_struc(n)%albedo_monthly, &
                                                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%z0 = LIS_roughness(n)%roughness(tid)
           else
              noah32_struc(n)%noah(t)%z0 = &
                   month_d(noah32_struc(n)%z0brd_monthly,      &
                   LIS_rc%yr,LIS_rc%mo,LIS_rc%da)
           endif

           noah32_struc(n)%noah(t)%z0_old = noah32_struc(n)%noah(t)%z0
        endif

        if (LIS_rc%usegreennessmap(n).eq."none") then
           gid = LIS_surface(n,mtype)%tile(t)%index
           call calc_localtime(LIS_rc%hr,LIS_domain(n)%grid(gid)%lon,  &
                               local_hour,change)
! For a true benchmark against the Noah3.2 testcase from NCAR,
! set "change = 0".  This line allows the code to _incorrectly_
! run in the same way as the Bondville testcase, which runs on
! local time instead of on UTC time. - dmm
           if (Bondvillecheck) then
              change = 0
              if (t.eq.1) &
              write(LIS_logunit,*) 'Performing benchmark of Noah3.2 ',&
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
           noah32_struc(n)%noah(t)%shdfac =                            &
                               month_d(noah32_struc(n)%shdfac_monthly, &
                               locyr,locmo,locda)
           noah32_struc(n)%noah(t)%shdmax = 0.0
           noah32_struc(n)%noah(t)%shdmin = 1.0
           do m = 1,12
              noah32_struc(n)%noah(t)%shdmax =                         &
                                  max(noah32_struc(n)%noah(t)%shdmax,  &
                                      noah32_struc(n)%shdfac_monthly(m))
              noah32_struc(n)%noah(t)%shdmin =                         &
                                  min(noah32_struc(n)%noah(t)%shdmin,  &
                                      noah32_struc(n)%shdfac_monthly(m))
           enddo
        else
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah32_struc(n)%noah(t)%shdfac =  LIS_gfrac(n)%greenness(tid)
        endif

        if (noah32_struc(n)%noah(t)%vegt.eq.LIS_rc%bareclass)          &
                 noah32_struc(n)%noah(t)%shdfac = 0.0
        if (noah32_struc(n)%noah(t)%vegt.eq.LIS_rc%snowclass)          &
                 noah32_struc(n)%noah(t)%shdfac = 0.0
        if (noah32_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass)         &
                 noah32_struc(n)%noah(t)%shdfac = 0.0

        if(noah32_struc(n)%noah(t)%shdfac.ge.noah32_struc(n)%noah(t)%shdmax) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah32_struc(n)%noah(t)%embrd =                             &
                   noah32_struc(n)%emissmax(noah32_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah32_struc(n)%noah(t)%z0brd =                             &
                   noah32_struc(n)%z0max(noah32_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah32_struc(n)%noah(t)%alb =                            &
                   noah32_struc(n)%albmin(noah32_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah32_struc(n)%noah(t)%lai =                            &
                   noah32_struc(n)%laimax(noah32_struc(n)%noah(t)%vegt) 
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        elseif(noah32_struc(n)%noah(t)%shdfac.le.noah32_struc(n)%noah(t)%shdmin) then
           if(LIS_rc%useemissmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
           else
              noah32_struc(n)%noah(t)%embrd =                             &
                   noah32_struc(n)%emissmin(noah32_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%useroughnessmap(n).ne."none") then 
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
           else
              noah32_struc(n)%noah(t)%z0brd =                             &
                   noah32_struc(n)%z0min(noah32_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
              noah32_struc(n)%noah(t)%alb =                            &
                    noah32_struc(n)%albmax(noah32_struc(n)%noah(t)%vegt)
           endif
           if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
              noah32_struc(n)%noah(t)%lai =                            &
                    noah32_struc(n)%laimin(noah32_struc(n)%noah(t)%vegt)
           else
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
           endif

        else
           if(noah32_struc(n)%noah(t)%shdmax.gt.noah32_struc(n)%noah(t)%shdmin) then
              interp_fraction =                                        &
                   (noah32_struc(n)%noah(t)%shdfac-noah32_struc(n)%noah(t)%shdmin) / &
                   (noah32_struc(n)%noah(t)%shdmax-noah32_struc(n)%noah(t)%shdmin)
              interp_fraction = min(interp_fraction,1.0)
              interp_fraction = max(interp_fraction,0.0)

              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah32_struc(n)%noah(t)%embrd = ((1.0 - interp_fraction) * &
                      noah32_struc(n)%emissmin(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%emissmax(noah32_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah32_struc(n)%noah(t)%z0brd = ((1.0 - interp_fraction) * &
                      noah32_struc(n)%z0min(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%z0max(noah32_struc(n)%noah(t)%vegt))
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah32_struc(n)%noah(t)%lai = ((1.0-interp_fraction) * &
                      noah32_struc(n)%laimin(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%laimax(noah32_struc(n)%noah(t)%vegt))
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah32_struc(n)%noah(t)%alb = ((1.0-interp_fraction) * &
                      noah32_struc(n)%albmax(noah32_struc(n)%noah(t)%vegt) + &
                      interp_fraction*noah32_struc(n)%albmin(noah32_struc(n)%noah(t)%vegt))
              endif
           else
              if(LIS_rc%useemissmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%embrd = LIS_emiss(n)%emiss(tid)
              else
                 noah32_struc(n)%noah(t)%embrd =                          &
                      0.5 * noah32_struc(n)%emissmin(noah32_struc(n)%noah(t)%vegt) + &
                      0.5 * noah32_struc(n)%emissmax(noah32_struc(n)%noah(t)%vegt)
              endif

              if(LIS_rc%useroughnessmap(n).ne."none") then 
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
              else
                 noah32_struc(n)%noah(t)%z0brd =                          &
                      0.5 * noah32_struc(n)%z0min(noah32_struc(n)%noah(t)%vegt) + &
                      0.5 * noah32_struc(n)%z0max(noah32_struc(n)%noah(t)%vegt)
              endif
              if(LIS_rc%uselaimap(n).eq."none") then ! no map specified
                 noah32_struc(n)%noah(t)%lai =                         &
                      0.5 * noah32_struc(n)%laimin(noah32_struc(n)%noah(t)%vegt) + &
                      0.5 * noah32_struc(n)%laimax(noah32_struc(n)%noah(t)%vegt)
              else
                 tid = LIS_surface(n,mtype)%tile(t)%tile_id
                 noah32_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
              endif
              if(LIS_rc%usealbedomap(n).eq."none") then ! no map specified
                 noah32_struc(n)%noah(t)%alb =                         &
                      0.5 * noah32_struc(n)%albmin(noah32_struc(n)%noah(t)%vegt) + &
                      0.5 * noah32_struc(n)%albmax(noah32_struc(n)%noah(t)%vegt)
              endif
           endif
        endif

        if ((LIS_rc%usealbedomap(n).ne."none").and.(LIS_rc%startcode.eq."coldstart")) then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah32_struc(n)%noah(t)%alb = LIS_alb(n)%albsf(tid)
           noah32_struc(n)%noah(t)%albedo = noah32_struc(n)%noah(t)%alb
        endif
        if (LIS_rc%uselaimap(n).ne."none") then
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah32_struc(n)%noah(t)%lai = LIS_lai(n)%tlai(tid)
        endif

        noah32_struc(n)%noah(t)%topt = topt_data
        noah32_struc(n)%noah(t)%cmcmax = cmcmax_data
        noah32_struc(n)%noah(t)%rsmax = rsmax_data
        noah32_struc(n)%noah(t)%cfactr = cfactr_data
        
     enddo

     if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
        write(LIS_logunit,*) 'Fixing Noah3.2 max snow albedo: ',       &
                              noah32_struc(n)%fixedmxsnalb
        if (noah32_struc(n)%fixedmxsnalb.eq.0.0) then
           write(LIS_logunit,*) 'WARNING: Noah3.2 max snow albedo ',   &
                                'set to zero!'
        endif
     else
        if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
           write(LIS_logunit,*) 'Reading Noah3.2 max snow albedo ',    &
                                'from max snow albedo maps'
        else
           write(LIS_logunit,*) 'Reading Noah3.2 max snow albedo ',    &
                                'from vegetation parameter file'
        endif
     endif

     do t=1,LIS_rc%npatch(n,mtype)
        noah32_struc(n)%noah(t)%nroot = &
             noah32_struc(n)%nrotbl(noah32_struc(n)%noah(t)%vegt)
        noah32_struc(n)%noah(t)%rsmin = &
             noah32_struc(n)%rstbl(noah32_struc(n)%noah(t)%vegt)
        noah32_struc(n)%noah(t)%rgl = &
             noah32_struc(n)%rgltbl(noah32_struc(n)%noah(t)%vegt)
        noah32_struc(n)%noah(t)%hs = &
             noah32_struc(n)%hstbl(noah32_struc(n)%noah(t)%vegt)
        noah32_struc(n)%noah(t)%snup = &
             noah32_struc(n)%snuptbl(noah32_struc(n)%noah(t)%vegt)

        if (LIS_rc%usemxsnalbmap(n).eq."fixed") then
           noah32_struc(n)%noah(t)%mxsnalb = noah32_struc(n)%fixedmxsnalb
        else
           if (LIS_rc%usemxsnalbmap(n).eq."LDT") then
              tid = LIS_surface(n,mtype)%tile(t)%tile_id
              noah32_struc(n)%noah(t)%mxsnalb = LIS_alb(n)%mxsnalb(tid)
           else
              noah32_struc(n)%noah(t)%mxsnalb = noah32_struc(n)%maxalb(noah32_struc(n)%noah(t)%vegt)/100.0
           endif
        endif
        if(LIS_rc%useroughnessmap(n).ne."none") then 
           tid = LIS_surface(n,mtype)%tile(t)%tile_id
           noah32_struc(n)%noah(t)%z0brd = LIS_roughness(n)%roughness(tid)
        else
           if (LIS_rc%startcode.eq."coldstart") then
              noah32_struc(n)%noah(t)%z0brd = noah32_struc(n)%noah(t)%z0
           endif
        endif
     enddo
  enddo

end subroutine noah32_setvegparms
