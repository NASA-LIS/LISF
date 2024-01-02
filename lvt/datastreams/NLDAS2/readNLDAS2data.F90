!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
!
! !ROUTINE: readNLDAS2data
! \label{readNLDAS2data}
!
! !INTERFACE:
  subroutine readNLDAS2data(source)
!
! !USES:
    use grib_api
    use ESMF
    use LVT_coreMod,    only : LVT_rc
    use LVT_logMod,     only : LVT_logunit, LVT_endrun, LVT_verify,  &
         LVT_getNextUnitNumber,                &
         LVT_releaseUnitNumber
    use LVT_timeMgrMod, only : LVT_tick, LVT_calendar
    use LVT_histDataMod
    use NLDAS2_dataMod, only : nldas2data

    implicit none
! !INPUT PARAMETERS:
!
    integer, intent(in)    :: source
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine provides the data reader for NLDAS-2 output data.
!  The plugin processes the land surface model output variables that
!  are specified in the MODEL OUTPUT TBL. The routine also spatially
!  interpolates the NLDAS-2 data to the model (LIS) output grid and
!  resolution.
!
! !FILES USED:
!
! !REVISION HISTORY:
!  09 Dec 2010: Sujay Kumar, Initial Specification
!  27 Jan 2014: David Mocko, Updates for NLDAS-2 SAC
!  16 May 2014: David Mocko, Updates for NLDAS-2 VIC from NASA GES DISC
!  11 Dec 2014: David Mocko, Added additional NLDAS-2 variables
!  17 Nov 2016: David Mocko, Added NLDAS-2 forcing variables
!
!EOP
!BOP
! !ARGUMENTS:
!
!EOP
    integer, parameter :: nc = 464
    integer, parameter :: nr = 224
    character*120      :: filename,sacsmname,streamflowname
    logical            :: file_exists
    integer            :: ftn,sacsmftn,streamflowftn,iret
    integer            :: c,r,t

    real :: swnet(nc*nr),lwnet(nc*nr)
    real :: qle(nc*nr),qh(nc*nr),qg(nc*nr),totalp(nc*nr)
    real :: snowf(nc*nr),rainf(nc*nr),evap(nc*nr),qs(nc*nr)
    real :: qsb(nc*nr),qsm(nc*nr),avgsurft(nc*nr),albedo(nc*nr)
    real :: sml1(nc*nr),sml2(nc*nr),sml3(nc*nr),sml4(nc*nr)
    real :: stl1(nc*nr),stl2(nc*nr),stl3(nc*nr),stl4(nc*nr)
    real :: potevap(nc*nr),ecanop(nc*nr),tveg(nc*nr),esoil(nc*nr)
    real :: rootmoist(nc*nr),canopint(nc*nr),subsnow(nc*nr)
    real :: snowcover(nc*nr),snowdepth(nc*nr),swe(nc*nr)
    real :: acond(nc*nr),ccond(nc*nr),streamflow(nc*nr)
    real :: precip_f(nc*nr),tair_f(nc*nr),swdown_f(nc*nr)
    real :: pet_f(nc*nr),cape_f(nc*nr),lwdown_f(nc*nr)
    real :: psurf_f(nc*nr),qair_f(nc*nr),crain_f(nc*nr)
    real :: nwind_f(nc*nr),ewind_f(nc*nr),fheight_f(nc*nr)
    real :: ecanopoverqle(nc*nr),tvegoverqle(nc*nr),esoiloverqle(nc*nr)
    real :: wind_f(nc*nr),tws(nc*nr)

    integer :: swnet_index,lwnet_index
    integer :: qle_index,qh_index,qg_index
    integer :: snowf_index,rainf_index,evap_index,qs_index
    integer :: qsb_index,qsm_index,avgsurft_index,albedo_index
    integer :: sml1_index,sml2_index,sml3_index,sml4_index
    integer :: stl1_index,stl2_index,stl3_index,stl4_index
    integer :: potevap_index,ecanop_index,tveg_index,esoil_index
    integer :: rootmoist_index,canopint_index,subsnow_index
    integer :: snowcover_index,snowdepth_index,swe_index
    integer :: acond_index,ccond_index
    integer :: precip_f_index,tair_f_index,swdown_f_index
    integer :: pet_f_index,cape_f_index,lwdown_f_index
    integer :: psurf_f_index,qair_f_index,crain_f_index
    integer :: nwind_f_index,ewind_f_index,fheight_f_index

    integer :: sml1_topLev,sml2_topLev,sml3_topLev,sml4_topLev
    integer :: sml1_botLev,sml2_botLev,sml3_botLev,sml4_botLev
    integer :: stl1_topLev,stl2_topLev,stl3_topLev,stl4_topLev
    integer :: stl1_botLev,stl2_botLev,stl3_botLev,stl4_botLev

    integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts

    real*8                  :: timenow
    real                    :: gmt1
    type(ESMF_Time)         :: time1

    integer                 :: nvars,index,igrib,status
    integer, allocatable    :: pid(:),tid(:),topLev(:),botLev(:)
    real, allocatable       :: var(:,:),sacsmvar(:,:),streamflowvar(:,:)
    real                    :: varfield(LVT_rc%lnc,LVT_rc%lnr)
    real                    :: rnet(LVT_rc%lnc, LVT_rc%lnr)
    real                    :: qs_ip(LVT_rc%lnc, LVT_rc%lnr)
    real                    :: qsb_ip(LVT_rc%lnc, LVT_rc%lnr)
    real                    :: trunoff_ip(LVT_rc%lnc, LVT_rc%lnr)

    ! GRIB IDs for NLDAS-2 forcing
    precip_f_index = 61
    tair_f_index = 11
    swdown_f_index = 204
    lwdown_f_index = 205
    pet_f_index = 228
    cape_f_index = 157
    psurf_f_index = 1
    qair_f_index = 51
    nwind_f_index = 34
    ewind_f_index = 33
    fheight_f_index = 7
    crain_f_index = 63

    ! GRIB IDs for all NLDAS-2 LSMs
    swnet_index = 111
    lwnet_index = 112
    qle_index = 121
    qh_index = 122
    qg_index = 155
    snowf_index = 161
    rainf_index = 162
    evap_index = 57
    qs_index = 235
    qsb_index = 234
    qsm_index = 99
    avgsurft_index = 148
    albedo_index = 84
    sml1_index = 86
    sml2_index = 86
    sml3_index = 86
    sml4_index = 86
    stl1_index = 85
    stl2_index = 85
    stl3_index = 85
    stl4_index = 85
    potevap_index = 145
    ecanop_index = 200
    tveg_index = 210
    esoil_index = 199
    rootmoist_index = 250
    canopint_index = 223
    subsnow_index = 198
    snowcover_index = 238
    snowdepth_index = 66
    swe_index = 65
    acond_index = 179
    ccond_index = 181

    ! NLDAS-2 soil levels
    if (nldas2data(source)%lsm.eq."NOAH") then
       sml1_topLev = 0
       sml1_botLev = 10
       sml2_topLev = 10
       sml2_botLev = 40
       sml3_topLev = 40
       sml3_botLev = 100
       sml4_topLev = 100
       sml4_botLev = 200
       stl1_topLev = 0
       stl1_botLev = 10
       stl2_topLev = 10
       stl2_botLev = 40
       stl3_topLev = 40
       stl3_botLev = 100
       stl4_topLev = 100
       stl4_botLev = 200
    elseif (nldas2data(source)%lsm.eq."MOS") then
       sml1_topLev = 0
       sml1_botLev = 10
       sml2_topLev = 10
       sml2_botLev = 40
       sml3_topLev = 40
       sml3_botLev = 200
       stl1_topLev = 0
       stl1_botLev = 0
    elseif (nldas2data(source)%lsm.eq."VIC") then
       sml1_topLev = 1
       sml1_botLev = 1
       sml2_topLev = 2
       sml2_botLev = 2
       sml3_topLev = 3
       sml3_botLev = 3
    endif

    ! Initialize time variables
    yr1 = LVT_rc%dyr(source)
    mo1 = LVT_rc%dmo(source)
    da1 = LVT_rc%dda(source)
    hr1 = LVT_rc%dhr(source)
    mn1 = 0
    ss1 = 0

    call ESMF_TimeSet(time1,yy=yr1,mm=mo1,dd=da1,h=hr1,m=mn1,          &
         s=ss1,calendar=LVT_calendar,rc=status)
    call LVT_verify(status)

! Tick backwards one-day for reading the NLDAS-2 monthly output.
! This is done so LVT will open the previous month's file when
! at the end of the month, so that the LVT output file will be
! the average over the previous NLDAS-2 monthly file.  - Mocko
    if (trim(nldas2data(source)%interval).eq."hourly") then
       ts = 0
    elseif (trim(nldas2data(source)%interval).eq."monthly") then
       ts = -86400
    endif
    call LVT_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts)

    ! Initialize NLDAS-2 variables
    swnet = LVT_rc%udef
    lwnet = LVT_rc%udef
    qle = LVT_rc%udef
    qh = LVT_rc%udef
    qg = LVT_rc%udef
    snowf = LVT_rc%udef
    rainf = LVT_rc%udef
    evap = LVT_rc%udef
    qs = LVT_rc%udef
    qsb = LVT_rc%udef
    qsm = LVT_rc%udef
    avgsurft = LVT_rc%udef
    albedo = LVT_rc%udef
    sml1 = LVT_rc%udef
    sml2 = LVT_rc%udef
    sml3 = LVT_rc%udef
    sml4 = LVT_rc%udef
    stl1 = LVT_rc%udef
    stl2 = LVT_rc%udef
    stl3 = LVT_rc%udef
    stl4 = LVT_rc%udef
    potevap = LVT_rc%udef
    ecanop = LVT_rc%udef
    tveg = LVT_rc%udef
    esoil = LVT_rc%udef
    rootmoist = LVT_rc%udef
    canopint = LVT_rc%udef
    subsnow = LVT_rc%udef
    snowcover = LVT_rc%udef
    snowdepth = LVT_rc%udef
    swe = LVT_rc%udef
    acond = LVT_rc%udef
    ccond = LVT_rc%udef
    streamflow = LVT_rc%udef
    precip_f = LVT_rc%udef
    tair_f = LVT_rc%udef
    swdown_f = LVT_rc%udef
    lwdown_f = LVT_rc%udef
    pet_f = LVT_rc%udef
    cape_f = LVT_rc%udef
    psurf_f = LVT_rc%udef
    qair_f = LVT_rc%udef
    crain_f = LVT_rc%udef
    nwind_f = LVT_rc%udef
    ewind_f = LVT_rc%udef
    wind_f = LVT_rc%udef ! EMK
    fheight_f = LVT_rc%udef
    ecanopoverqle = LVT_rc%udef
    tvegoverqle =LVT_rc%udef
    esoiloverqle = LVT_rc%udef
    tws = LVT_rc%udef

    if (nldas2data(source)%lsm.eq."SAC") then
       allocate(sacsmvar(5,nc*nr))
    endif
    if (LVT_MOC_STREAMFLOW(source).ge.1) then
       allocate(streamflowvar(1,nc*nr))
    endif

    if (nldas2data(source)%anlys_data_class.eq."Routing") then
       if (LVT_MOC_STREAMFLOW(source).ge.1) then
          call create_nldas2data_filename(nldas2data(source)%odir,     &
               nldas2data(source)%interval,                            &
               trim(nldas2data(source)%lsm) //                         &
               "ST",yr1,mo1,da1,hr1,doy1,streamflowname)
          write(LVT_logunit,*) '[INFO] Reading ',trim(streamflowname)
          streamflowftn = LVT_getNextUnitNumber()
          open(streamflowftn,file=trim(streamflowname),                &
               status='old',access='sequential',                       &
               form='unformatted',iostat=iret)
          if (iret.ne.0) then
             write(LVT_logunit,*)                                      &
                  '[ERR] Stopping routine: File not opened: ',iret
             call LVT_endrun
          endif
          read(streamflowftn) streamflowvar(1,:)
          close(streamflowftn)
          call LVT_releaseUnitNumber(streamflowftn)
       endif

       do r = 1,nr
          do c = 1,nc
             if (streamflowvar(1,c+(r-1)*nc).gt.-9998.0) then
                streamflow(c+(r-1)*nc) =                               &
                     streamflowvar(1,c+(r-1)*nc)
             endif
          enddo
       enddo
    endif

    if (nldas2data(source)%anlys_data_class.eq."LSM") then
       call create_nldas2data_filename(nldas2data(source)%odir,        &
                                       nldas2data(source)%interval,    &
            trim(nldas2data(source)%lsm),yr1,mo1,da1,hr1,doy1,filename)
       inquire(file=trim(filename),exist=file_exists)
       if (file_exists) then
          write(LVT_logunit,*) '[INFO] Reading ',trim(filename)
          call grib_open_file(ftn,trim(filename),'r',iret)
          if (iret.ne.0) then
             write(LVT_logunit,*)                                      &
                  '[ERR] Stopping routine: File not opened : ',iret
             call LVT_endrun
          endif
       else
          write(LVT_logunit,*) '[INFO] Cannot find ',trim(filename)
       endif

       if (((LVT_MOC_SOILMOIST(source).ge.1).or.                       &
            (LVT_MOC_ROOTMOIST(source).ge.1).or.                       &
            (LVT_MOC_TWS(source).ge.1))                                &
            .and.(nldas2data(source)%lsm.eq."SAC")) then
          call create_nldas2data_filename(nldas2data(source)%odir,     &
                                          nldas2data(source)%interval, &
                                "SACSM",yr1,mo1,da1,hr1,doy1,sacsmname)
          write(LVT_logunit,*) '[INFO] Reading ',trim(sacsmname)
          sacsmftn = LVT_getNextUnitNumber()
          open(sacsmftn,file=trim(sacsmname),status='old',             &
               access='direct',form='unformatted',                     &
               recl=nc*nr*4,iostat=iret)
          if (iret.ne.0) then
             write(LVT_logunit,*)                                      &
                  '[ERR] Stopping routine: File not opened: ',iret
             call LVT_endrun
          endif
          read(sacsmftn,rec=1) sacsmvar(1,:)
          read(sacsmftn,rec=2) sacsmvar(2,:)
          read(sacsmftn,rec=3) sacsmvar(3,:)
          read(sacsmftn,rec=4) sacsmvar(4,:)
          read(sacsmftn,rec=5) sacsmvar(5,:)
          close(sacsmftn)
          call LVT_releaseUnitNumber(sacsmftn)

          if (nldas2data(source)%smvol) then
             do r = 1,nr
                do c = 1,nc
                   if (sacsmvar(1,c+(r-1)*nc).lt.9.99E+20) then
                      sml1(c+(r-1)*nc) =                               &
                           sacsmvar(1,c+(r-1)*nc) /  100.0
                   endif
                   if (sacsmvar(2,c+(r-1)*nc).lt.9.99E+20) then
                      sml2(c+(r-1)*nc) =                               &
                           sacsmvar(2,c+(r-1)*nc) /  300.0
                   endif
                   if (sacsmvar(3,c+(r-1)*nc).lt.9.99E+20) then
                      sml3(c+(r-1)*nc) =                               &
                           sacsmvar(3,c+(r-1)*nc) /  600.0
                   endif
                   if (sacsmvar(4,c+(r-1)*nc).lt.9.99E+20) then
                      sml4(c+(r-1)*nc) =                               &
                           sacsmvar(4,c+(r-1)*nc) / 1000.0
                   endif
                   if (sacsmvar(5,c+(r-1)*nc).lt.9.99E+20) then
                      rootmoist(c+(r-1)*nc) =                          &
                           sacsmvar(5,c+(r-1)*nc) / 1000.0
                   endif
                enddo
             enddo
          else
             do r = 1,nr
                do c = 1,nc
                   if (sacsmvar(1,c+(r-1)*nc).lt.9.99E+20) then
                      sml1(c+(r-1)*nc) =                               &
                           sacsmvar(1,c+(r-1)*nc)
                   endif
                   if (sacsmvar(2,c+(r-1)*nc).lt.9.99E+20) then
                      sml2(c+(r-1)*nc) =                               &
                           sacsmvar(2,c+(r-1)*nc)
                   endif
                   if (sacsmvar(3,c+(r-1)*nc).lt.9.99E+20) then
                      sml3(c+(r-1)*nc) =                               &
                           sacsmvar(3,c+(r-1)*nc)
                   endif
                   if (sacsmvar(4,c+(r-1)*nc).lt.9.99E+20) then
                      sml4(c+(r-1)*nc) =                               &
                           sacsmvar(4,c+(r-1)*nc)
                   endif
                   if (sacsmvar(5,c+(r-1)*nc).lt.9.99E+20) then
                      rootmoist(c+(r-1)*nc) =                          &
                           sacsmvar(5,c+(r-1)*nc)
                   endif
                enddo
             enddo
          endif
       endif

       call grib_count_in_file(ftn,nvars)

       allocate(var(nvars,nc*nr))
       allocate(pid(nvars))
       allocate(tid(nvars))
       allocate(topLev(nvars))
       allocate(botLev(nvars))

       do index = 1,nvars
          call grib_new_from_file(ftn,igrib,iret)
          call LVT_verify(iret,                                        &
               'grib_new_from_file failed in readNLDAS2data')

          call grib_get(igrib,"indicatorOfParameter",pid(index),iret)
          call LVT_verify(iret,                                        &
               'grib_get failed for indicatorOfParameter in readNLDAS2data')

          call grib_get(igrib,"timeRangeIndicator",tid(index),iret)
          call LVT_verify(iret,                                        &
               'grib_get failed for timeRangeIndicator in readNLDAS2data')

          call grib_get(igrib, "topLevel",topLev(index),iret)
          call LVT_verify(iret,                                        &
               'grib_get failed for topLevel in readNLDAS2data')

          call grib_get(igrib, "bottomLevel",botLev(index),iret)
          call LVT_verify(iret,                                        &
               'grib_get failed for bottomLevel in readNLDAS2data')

          if ((pid(index).eq.precip_f_index).and.                      &
              (LVT_MOC_TOTALPRECIP(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,precip_f)
          endif

          if ((pid(index).eq.tair_f_index).and.                        &
              (LVT_MOC_TAIRFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,tair_f)
          endif

          if ((pid(index).eq.swdown_f_index).and.                      &
              (LVT_MOC_SWDOWNFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,swdown_f)
          endif

          if ((pid(index).eq.lwdown_f_index).and.                      &
              (LVT_MOC_LWDOWNFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,lwdown_f)
          endif

          if ((pid(index).eq.pet_f_index).and.                         &
              (LVT_MOC_PETFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,pet_f)
          endif

          if ((pid(index).eq.cape_f_index).and.                        &
              (LVT_MOC_CAPEFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,cape_f)
          endif

          if ((pid(index).eq.psurf_f_index).and.                       &
              (LVT_MOC_PSURFFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,psurf_f)
          endif

          if ((pid(index).eq.qair_f_index).and.                        &
              (LVT_MOC_QAIRFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,qair_f)
          endif

          if ((pid(index).eq.crain_f_index).and.                       &
              (LVT_MOC_CRAINFFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,crain_f)
          endif

          if ((pid(index).eq.nwind_f_index).and.                       &
              ((LVT_MOC_NWINDFORC(source).ge.1).or.                    &
               (LVT_MOC_WINDFORC(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,nwind_f)
          endif

          if ((pid(index).eq.ewind_f_index).and.                       &
              ((LVT_MOC_EWINDFORC(source).ge.1).or.                    &
               (LVT_MOC_WINDFORC(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,ewind_f)
          endif

          if ((pid(index).eq.fheight_f_index).and.                     &
              (LVT_MOC_FHEIGHTFORC(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,fheight_f)
          endif

          if ((pid(index).eq.swnet_index).and.                         &
              (LVT_MOC_SWNET(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,swnet)
          endif

          if ((pid(index).eq.lwnet_index).and.                         &
              (LVT_MOC_LWNET(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,lwnet)
          endif

          if ((pid(index).eq.qle_index).and.                           &
              ((LVT_MOC_QLE(source).ge.1).or.                          &
               (LVT_MOC_RNET(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,qle)
          endif

          if ((pid(index).eq.qh_index).and.                            &
              ((LVT_MOC_QH(source).ge.1).or.                           &
               (LVT_MOC_RNET(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,qh)
          endif

          if ((pid(index).eq.qg_index).and.                            &
              ((LVT_MOC_QG(source).ge.1).or.                           &
               (LVT_MOC_RNET(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,qg)
          endif

          if ((pid(index).eq.snowf_index).and.                         &
              ((LVT_MOC_SNOWF(source).ge.1).or.                        &
               (LVT_MOC_TOTALPRECIP(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,snowf)
          endif

          if ((pid(index).eq.rainf_index).and.                         &
              ((LVT_MOC_RAINF(source).ge.1).or.                        &
               (LVT_MOC_TOTALPRECIP(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,rainf)
          endif

          if ((pid(index).eq.evap_index).and.                          &
              (LVT_MOC_EVAP(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,evap)
          endif

          if ((pid(index).eq.qs_index).and.                            &
              ((LVT_MOC_QS(source).ge.1).or.                           &
               (LVT_MOC_RUNOFF(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,qs)
          endif

          if ((pid(index).eq.qsb_index).and.                           &
              ((LVT_MOC_QSB(source).ge.1).or.                          &
               (LVT_MOC_RUNOFF(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,qsb)
          endif

          if ((pid(index).eq.qsm_index).and.                           &
              (LVT_MOC_QSM(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,qsm)
          endif

          if ((pid(index).eq.avgsurft_index).and.                      &
              ((LVT_MOC_AVGSURFT(source).ge.1).or.                     &
               (LVT_MOC_RADT(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,avgsurft)
          endif

          if ((pid(index).eq.albedo_index).and.                        &
              (LVT_MOC_ALBEDO(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,albedo)
          endif

          if ((pid(index).eq.swe_index).and.                           &
              ((LVT_MOC_SWE(source).ge.1).or.                          &
               (LVT_MOC_TWS(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,swe)
             if (nldas2data(source)%lsm.eq."SAC") then !SAC has units of m 
                do r = 1,nr
                   do c = 1,nc
                      if (swe(c+(r-1)*nc).ne.9999.0) then
                         swe(c+(r-1)*nc) = swe(c+(r-1)*nc) * 1000.0
                      endif
                   enddo
                enddo
             endif
          endif

          if ((pid(index).eq.snowdepth_index).and.                     &
              (LVT_MOC_SNOWDEPTH(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,snowdepth)
          endif

          if ((pid(index).eq.snowcover_index).and.                     &
              (LVT_MOC_SNOWCOVER(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,snowcover)
          endif

          if ((pid(index).eq.potevap_index).and.                       &
              (LVT_MOC_POTEVAP(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,potevap)
          endif

          if ((pid(index).eq.ecanop_index).and.                        &
              ((LVT_MOC_ECANOP(source).ge.1).or.                       &
               (LVT_MOC_ECANOPOVERQLE(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,ecanop)

             if (nldas2data(source)%lsm.eq."MOS") then
                do r = 1,nr
                   do c = 1,nc
                      if (ecanop(c+(r-1)*nc).ne.9999.0) then
                         ecanop(c+(r-1)*nc) = -1.0*ecanop(c+(r-1)*nc)
                      endif
                   enddo
                enddo
             endif
          endif

          if ((pid(index).eq.tveg_index).and.                          &
              ((LVT_MOC_TVEG(source).ge.1).or.                         &
               (LVT_MOC_TVEGOVERQLE(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,tveg)

             if (nldas2data(source)%lsm.eq."MOS") then
                do r = 1,nr
                   do c = 1,nc
                      if (tveg(c+(r-1)*nc).ne.9999.0) then
                         tveg(c+(r-1)*nc) = -1.0*tveg(c+(r-1)*nc)
                      endif
                   enddo
                enddo
             endif
          endif

          if ((pid(index).eq.esoil_index).and.                         &
              ((LVT_MOC_ESOIL(source).ge.1).or.                        &
               (LVT_MOC_ESOILOVERQLE(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,esoil)
             if (nldas2data(source)%lsm.eq."MOS") then
                do r = 1,nr
                   do c = 1,nc
                      if (esoil(c+(r-1)*nc).ne.9999.0) then
                         esoil(c+(r-1)*nc) = -1.0*esoil(c+(r-1)*nc)
                      endif
                   enddo
                enddo
             endif
          endif

          if ((pid(index).eq.canopint_index).and.                      &
              ((LVT_MOC_CANOPINT(source).ge.1).or.                     &
               (LVT_MOC_TWS(source).ge.1))) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,canopint)
          endif

          if ((pid(index).eq.subsnow_index).and.(LVT_MOC_SUBSNOW(source).ge.1)) then
             call retrieve_nldas2data(igrib,nc,nr,nvars,index,subsnow)
          endif

          if (nldas2data(source)%lsm.ne."SAC") then
             if ((pid(index).eq.sml1_index).and.                       &
                  (topLev(index).eq.sml1_topLev).and.                  &
                  (botLev(index).eq.sml1_botLev).and.                  &
                  ((LVT_MOC_SOILMOIST(source).ge.1).or.                &
                   (LVT_MOC_ROOTMOIST(source).ge.1).or.                &
                   (LVT_MOC_TWS(source).ge.1))) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,sml1)
             endif
             if ((pid(index).eq.sml2_index).and.                       &
                  (topLev(index).eq.sml2_topLev).and.                  &
                  (botLev(index).eq.sml2_botLev).and.                  &
                  ((LVT_MOC_SOILMOIST(source).ge.1).or.                &
                   (LVT_MOC_ROOTMOIST(source).ge.1).or.                &
                   (LVT_MOC_TWS(source).ge.1))) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,sml2)
             endif
             if ((pid(index).eq.sml3_index).and.                       &
                  (topLev(index).eq.sml3_topLev).and.                  &
                  (botLev(index).eq.sml3_botLev).and.                  &
                  ((LVT_MOC_SOILMOIST(source).ge.1).or.                &
                   (LVT_MOC_ROOTMOIST(source).ge.1).or.                &
                   (LVT_MOC_TWS(source).ge.1))) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,sml3)
             endif
             if ((pid(index).eq.stl1_index).and.                       &
                  (topLev(index).eq.stl1_topLev).and.                  &
                  (botLev(index).eq.stl1_botLev).and.                  &
                  (LVT_MOC_SOILTEMP(source).ge.1)) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,stl1)
             endif
! Grab the NLDAS-2 top-1 meter soil moisture directly. - Mocko
             if ((pid(index).eq.sml1_index).and.                       &
                  (topLev(index).eq.  0).and.                          &
                  (botLev(index).eq.100).and.                          &
                  ((LVT_MOC_SOILMOIST(source).ge.1).or.                &
                   (LVT_MOC_ROOTMOIST(source).ge.1))) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,rootmoist)
             endif
          endif

          if (nldas2data(source)%lsm.eq."NOAH") then
             if ((pid(index).eq.sml4_index).and.                       &
                  (topLev(index).eq.sml4_topLev).and.                  &
                  (botLev(index).eq.sml4_botLev).and.                  &
                  ((LVT_MOC_SOILMOIST(source).ge.1).or.                &
                   (LVT_MOC_ROOTMOIST(source).ge.1).or.                &
                   (LVT_MOC_TWS(source).ge.1))) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,sml4)
             endif
             if ((pid(index).eq.stl2_index).and.                       &
                  (topLev(index).eq.stl2_topLev).and.                  &
                  (botLev(index).eq.stl2_botLev).and.                  &
                  (LVT_MOC_SOILTEMP(source).ge.1)) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,stl2)
             endif
             if ((pid(index).eq.stl3_index).and.                       &
                  (topLev(index).eq.stl3_topLev).and.                  &
                  (botLev(index).eq.stl3_botLev).and.                  &
                  (LVT_MOC_SOILTEMP(source).ge.1)) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,stl3)
             endif
             if ((pid(index).eq.stl4_index).and.                       &
                  (topLev(index).eq.stl4_topLev).and.                  &
                  (botLev(index).eq.stl4_botLev).and.                  &
                  (LVT_MOC_SOILTEMP(source).ge.1)) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,stl4)
             endif
          endif

          if (nldas2data(source)%lsm.eq."VIC") then
             if ((pid(index).eq.stl2_index).and.                       &
                  (topLev(index).eq.stl2_topLev).and.                  &
                  (botLev(index).eq.stl2_botLev).and.                  &
                  (LVT_MOC_SOILTEMP(source).ge.1)) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,stl2)
             endif
             if ((pid(index).eq.stl3_index).and.                       &
                  (topLev(index).eq.stl3_topLev).and.                  &
                  (botLev(index).eq.stl3_botLev).and.                  &
                  (LVT_MOC_SOILTEMP(source).ge.1)) then
                call retrieve_nldas2data(igrib,nc,nr,nvars,index,stl3)
             endif
          endif

          call grib_release(igrib,iret)
          call LVT_verify(iret,'grib_release failed in readNLDAS2data')
       enddo

       call grib_close_file(ftn,iret)

       deallocate(var)
       deallocate(pid)
       deallocate(tid)
       deallocate(topLev)
       deallocate(botLev)
    endif

    if (nldas2data(source)%lsm.eq."SAC") then
       deallocate(sacsmvar)
    endif
    if (LVT_MOC_STREAMFLOW(source).ge.1) then
       deallocate(streamflowvar)
    endif

    if (nldas2data(source)%lsm.ne."SAC") then
       call interp_nldas2var(source,nc,nr,swnet,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SWNET,source,varfield,  &
                vlevel=1,units="W/m2")
       call interp_nldas2var(source,nc,nr,lwnet,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_LWNET,source,varfield,  &
                vlevel=1,units="W/m2")
       call interp_nldas2var(source,nc,nr,qle,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,varfield,    &
                vlevel=1,units="W/m2")
       rnet = varfield
       call interp_nldas2var(source,nc,nr,qh,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_QH,source,varfield,     &
                vlevel=1,units="W/m2")
       ! Use LVT grid dimensions for these arrays
       do r = 1,lvt_rc%lnr
          do c = 1,lvt_rc%lnc
             if (varfield(c,r).ne.LVT_rc%udef) then
                rnet(c,r) = rnet(c,r) + varfield(c,r)
             endif
          enddo
       enddo
       call interp_nldas2var(source,nc,nr,qg,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_QG,source,varfield,     &
                vlevel=1,units="W/m2")
       ! Use LVT grid dimensions for these arrays
       do r = 1,lvt_rc%lnr
          do c = 1,lvt_rc%lnc
             if (varfield(c,r).ne.LVT_rc%udef) then
                rnet(c,r) = rnet(c,r) + varfield(c,r)
             endif
          enddo
       enddo
       call LVT_logSingleDataStreamVar(LVT_MOC_RNET,source,rnet,       &
                vlevel=1,units="W/m2")
       call interp_nldas2var(source,nc,nr,avgsurft,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,varfield,&
                vlevel=1,units="K")
       call LVT_logSingleDataStreamVar(LVT_MOC_RADT,source,varfield,   &
                vlevel=1,units="K")
       if (LVT_MOC_ALBEDO(source).ge.1) then
          do r = 1,nr
             do c = 1,nc
                if (albedo(c+(r-1)*nc).gt.-9998.0) then
                   albedo(c+(r-1)*nc) = albedo(c+(r-1)*nc) / 100.0
                endif
             enddo
          enddo
       endif
       call interp_nldas2var(source,nc,nr,albedo,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_ALBEDO,source,varfield, &
                vlevel=1,units="-")
    endif

    if ((nldas2data(source)%lsm.ne."FORCINGA").and.                    &
        (nldas2data(source)%lsm.ne."FORCINGB")) then
       if (LVT_MOC_TOTALPRECIP(source).ge.1) then
          do r = 1,nr
             do c = 1,nc
                if (rainf(c+(r-1)*nc).gt.-9998.0) then
                   totalp(c+(r-1)*nc) = rainf(c+(r-1)*nc) + snowf(c+(r-1)*nc)
                endif
             enddo
          enddo
       endif
       call interp_nldas2var(source,nc,nr,totalp,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source,varfield,&
                vlevel=1,units="kg/m2")
       if (LVT_MOC_TOTALPRECIP(source).ge.1) then
          do r = 1,nr
             do c = 1,nc
                if (totalp(c+(r-1)*nc).gt.-9998.0) then
                   totalp(c+(r-1)*nc) = totalp(c+(r-1)*nc) / 3600.0
                endif
             enddo
          enddo
       endif
       call interp_nldas2var(source,nc,nr,totalp,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source,varfield,&
                vlevel=1,units="kg/m2s")
    endif

    call interp_nldas2var(source,nc,nr,snowf,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_SNOWF,source,varfield,     &
             vlevel=1,units="kg/m2")
    if (LVT_MOC_SNOWF(source).ge.1) then
       do r = 1,nr
          do c = 1,nc
             if (snowf(c+(r-1)*nc).gt.-9998.0) then
                snowf(c+(r-1)*nc) = snowf(c+(r-1)*nc) / 3600.0
             endif
          enddo
       enddo
    endif
    call interp_nldas2var(source,nc,nr,snowf,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_SNOWF,source,varfield,     &
             vlevel=1,units="kg/m2s")
    call interp_nldas2var(source,nc,nr,rainf,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_RAINF,source,varfield,     &
             vlevel=1,units="kg/m2")
    if (LVT_MOC_RAINF(source).ge.1) then
       do r = 1,nr
          do c = 1,nc
             if (rainf(c+(r-1)*nc).gt.-9998.0) then
                rainf(c+(r-1)*nc) = rainf(c+(r-1)*nc) / 3600.0
             endif
          enddo
       enddo
    endif
    call interp_nldas2var(source,nc,nr,rainf,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_RAINF,source,varfield,     &
             vlevel=1,units="kg/m2s")

    if (LVT_MOC_EVAP(source).ge.1) then
       do r = 1,nr
          do c = 1,nc
             if (evap(c+(r-1)*nc).gt.-9998.0) then
                evap(c+(r-1)*nc) = evap(c+(r-1)*nc) / 3600.0
             endif
          enddo
       enddo
    endif
    call interp_nldas2var(source,nc,nr,evap,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_EVAP,source,varfield,      &
             vlevel=1,units="kg/m2s")

    if ((LVT_MOC_QS(source).ge.1).or.(LVT_MOC_RUNOFF(source).ge.1)) then
       do r = 1,nr
          do c = 1,nc
             if (qs(c+(r-1)*nc).gt.-9998.0) then
                qs(c+(r-1)*nc) = qs(c+(r-1)*nc) / 3600.0
             endif
          enddo
       enddo
    endif
    call interp_nldas2var(source,nc,nr,qs,qs_ip)
    call LVT_logSingleDataStreamVar(LVT_MOC_QS,source,qs_ip,           &
             vlevel=1,units="kg/m2s")
    if ((LVT_MOC_QSB(source).ge.1).or.(LVT_MOC_RUNOFF(source).ge.1)) then
       do r = 1,nr
          do c = 1,nc
             if (qsb(c+(r-1)*nc).gt.-9998.0) then
                qsb(c+(r-1)*nc) = qsb(c+(r-1)*nc) / 3600.0
             endif
          enddo
       enddo
    endif
    call interp_nldas2var(source,nc,nr,qsb,qsb_ip)
    call LVT_logSingleDataStreamVar(LVT_MOC_QSB,source,qsb_ip,         &
             vlevel=1,units="kg/m2s")
    ! Use LVT grid dimensions for these arrays
    do r = 1,lvt_rc%lnr
       do c = 1,lvt_rc%lnc
          if (qs_ip(c,r).ne.LVT_rc%udef) then
             trunoff_ip(c,r) = qs_ip(c,r) + qsb_ip(c,r)
          else
             trunoff_ip(c,r) = LVT_rc%udef
          endif
       enddo
    enddo
    call LVT_logSingleDataStreamVar(LVT_MOC_RUNOFF,source,trunoff_ip,  &
             vlevel=1,units="kg/m2s")
    ! Use LVT grid dimensions for these arrays
    do r = 1,lvt_rc%lnr
       do c = 1,lvt_rc%lnc
          if (trunoff_ip(c,r).ne.LVT_rc%udef) then
             trunoff_ip(c,r) = trunoff_ip(c,r) * 86400.0 * 30.0
          else
             trunoff_ip(c,r) = LVT_rc%udef
          endif
       enddo
    enddo
    call LVT_logSingleDataStreamVar(LVT_MOC_RUNOFF,source,trunoff_ip,  &
             vlevel=1,units="mm/month")

    if (LVT_MOC_QSM(source).ge.1) then
       do r = 1,nr
          do c = 1,nc
             if (qsm(c+(r-1)*nc).gt.-9998.0) then
                qsm(c+(r-1)*nc) = qsm(c+(r-1)*nc) / 3600.0
             endif
          enddo
       enddo
    endif

    if (LVT_MOC_ECANOPOVERQLE(source).ge.1) then
       do r = 1,nr
          do c = 1,nc
             if (ecanop(c+(r-1)*nc).gt.-9998.0) then
                if ((ecanop(c+(r-1)*nc) + tveg(c+(r-1)*nc) +           &
                     esoil(c+(r-1)*nc)).ne.0) then
                   ecanopoverqle(c+(r-1)*nc) = ecanop(c+(r-1)*nc) /    &
                        ((ecanop(c+(r-1)*nc) + tveg(c+(r-1)*nc) +      &
                          esoil(c+(r-1)*nc)))
                   if (abs(ecanopoverqle(c+(r-1)*nc)).gt.1)            &
                        ecanopoverqle(c+(r-1)*nc) = LVT_rc%udef
                endif
             endif
          enddo
       enddo
       call interp_nldas2var(source,nc,nr,ecanopoverqle,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_ECANOPOVERQLE,source,varfield,&
                vlevel=1,units="-")
    endif

    if (LVT_MOC_TVEGOVERQLE(source).ge.1) then
       do r = 1,nr
          do c = 1,nc
             if (tveg(c+(r-1)*nc).gt.-9998.0) then
                if ((ecanop(c+(r-1)*nc) + tveg(c+(r-1)*nc) +           &
                     esoil(c+(r-1)*nc)).ne.0) then
                   tvegoverqle(c+(r-1)*nc) = tveg(c+(r-1)*nc) /        &
                        ((ecanop(c+(r-1)*nc) + tveg(c+(r-1)*nc) +      &
                          esoil(c+(r-1)*nc)))
                   if (abs(tvegoverqle(c+(r-1)*nc)).gt.1)              &
                        tvegoverqle(c+(r-1)*nc) = LVT_rc%udef
                endif
             endif
          enddo
       enddo
       call interp_nldas2var(source,nc,nr,tvegoverqle,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_TVEGOVERQLE,source,varfield,&
                vlevel=1,units="-")
    endif

    if (LVT_MOC_ESOILOVERQLE(source).ge.1) then
       do r = 1,nr
          do c = 1,nc
             if (esoil(c+(r-1)*nc).gt.-9998.0) then
                if ((ecanop(c+(r-1)*nc) + tveg(c+(r-1)*nc) +           &
                     esoil(c+(r-1)*nc)).ne.0) then
                   esoiloverqle(c+(r-1)*nc) = esoil(c+(r-1)*nc) /      &
                        ((ecanop(c+(r-1)*nc) + tveg(c+(r-1)*nc) +      &
                          esoil(c+(r-1)*nc)))
                   if (abs(esoiloverqle(c+(r-1)*nc)).gt.1)             &
                        esoiloverqle(c+(r-1)*nc) = LVT_rc%udef
                endif
             endif
          enddo
       enddo
       call interp_nldas2var(source,nc,nr,esoiloverqle,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_ESOILOVERQLE,source,varfield,&
                vlevel=1,units="-")
    endif

    if (LVT_MOC_TWS(source).ge.1) then
       if (nldas2data(source)%lsm.eq."SAC") then
          do r = 1,nr
             do c = 1,nc
                if (sml1(c+(r-1)*nc).gt.0) then
                   tws(c+(r-1)*nc) = sml1(c+(r-1)*nc) +                &
                                     sml2(c+(r-1)*nc) +                &
                                     sml3(c+(r-1)*nc) +                &
                                     sml4(c+(r-1)*nc) +                &
                                     swe(c+(r-1)*nc)
                endif
             enddo
          enddo
       elseif ((nldas2data(source)%lsm.eq."MOS").or.                   &
               (nldas2data(source)%lsm.eq."VIC")) then
          do r = 1,nr
             do c = 1,nc
                if (sml1(c+(r-1)*nc).gt.0) then
                   tws(c+(r-1)*nc) = sml1(c+(r-1)*nc) +                &
                                     sml2(c+(r-1)*nc) +                &
                                     sml3(c+(r-1)*nc) +                &
                                     swe(c+(r-1)*nc)
                endif
             enddo
          enddo
       else
          do r = 1,nr
             do c = 1,nc
                if (sml1(c+(r-1)*nc).gt.0) then
                   tws(c+(r-1)*nc) = sml1(c+(r-1)*nc) +                &
                                     sml2(c+(r-1)*nc) +                &
                                     sml3(c+(r-1)*nc) +                &
                                     sml4(c+(r-1)*nc) +                &
                                     swe(c+(r-1)*nc) +                 &
                                     canopint(c+(r-1)*nc)
                endif
             enddo
          enddo
       endif
       call interp_nldas2var(source,nc,nr,tws,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_TWS,source,varfield,    &
                vlevel=1,units="mm")
    endif

    call interp_nldas2var(source,nc,nr,qsm,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_QSM,source,varfield,       &
             vlevel=1,units="kg/m2s")
    call interp_nldas2var(source,nc,nr,swe,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,varfield,       &
             vlevel=1,units="kg/m2")
    call interp_nldas2var(source,nc,nr,snowdepth,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH,source,varfield, &
             vlevel=1,units="m")
    call interp_nldas2var(source,nc,nr,snowcover,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_SNOWCOVER,source,varfield, &
             vlevel=1,units="-")
    call interp_nldas2var(source,nc,nr,potevap,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_POTEVAP,source,varfield,   &
             vlevel=1,units="W/m2")
    call interp_nldas2var(source,nc,nr,ecanop,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_ECANOP,source,varfield,    &
             vlevel=1,units="W/m2")
    call interp_nldas2var(source,nc,nr,tveg,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_TVEG,source,varfield,      &
             vlevel=1,units="W/m2")
    call interp_nldas2var(source,nc,nr,esoil,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_ESOIL,source,varfield,     &
             vlevel=1,units="W/m2")
    call interp_nldas2var(source,nc,nr,canopint,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_CANOPINT,source,varfield,  &
             vlevel=1,units="kg/m2")
    call interp_nldas2var(source,nc,nr,subsnow,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_SUBSNOW,source,varfield,   &
             vlevel=1,units="W/m2")
    if ((nldas2data(source)%lsm.eq."FORCINGA").or.                     &
        (nldas2data(source)%lsm.eq."FORCINGB")) then
       call interp_nldas2var(source,nc,nr,precip_f,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source,varfield,&
                vlevel=1,units="kg/m2")
       if (LVT_MOC_TOTALPRECIP(source).ge.1) then
          do r = 1,nr
             do c = 1,nc
                if (precip_f(c+(r-1)*nc).gt.-9998.0) then
                   precip_f(c+(r-1)*nc) = precip_f(c+(r-1)*nc) / 3600.0
                endif
             enddo
          enddo
       endif
       call interp_nldas2var(source,nc,nr,precip_f,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source,varfield,&
                vlevel=1,units="kg/m2s")
    endif
    call interp_nldas2var(source,nc,nr,tair_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,source,varfield,  &
             vlevel=1,units="K")
    call interp_nldas2var(source,nc,nr,swdown_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_SWDOWNFORC,source,varfield,&
             vlevel=1,units="W/m2")
    call interp_nldas2var(source,nc,nr,lwdown_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_LWDOWNFORC,source,varfield,&
             vlevel=1,units="W/m2")
    call interp_nldas2var(source,nc,nr,pet_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_PETFORC,source,varfield,   &
             vlevel=1,units="kg/m2")
    call interp_nldas2var(source,nc,nr,cape_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_CAPEFORC,source,varfield,  &
             vlevel=1,units="J/kg")
    call interp_nldas2var(source,nc,nr,psurf_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_PSURFFORC,source,varfield, &
             vlevel=1,units="Pa")
    call interp_nldas2var(source,nc,nr,qair_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_QAIRFORC,source,varfield,  &
             vlevel=1,units="kg/kg")
    call interp_nldas2var(source,nc,nr,crain_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_CRAINFFORC,source,varfield,&
             vlevel=1,units="kg/m2")
    call interp_nldas2var(source,nc,nr,nwind_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_NWINDFORC,source,varfield, &
             vlevel=1,units="m/s")
    call interp_nldas2var(source,nc,nr,ewind_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_EWINDFORC,source,varfield, &
             vlevel=1,units="m/s")

    ! EMK...Create total wind speed forcing in m/s
    do r = 1,nr
       do c = 1,nc
          if ((ewind_f(c+(r-1)*nc).gt.-9998).and.                      &
              (nwind_f(c+(r-1)*nc).gt.-9998)) then
             wind_f(c+(r-1)*nc) = sqrt(ewind_f(c+(r-1)*nc)**2 +        &
                                       nwind_f(c+(r-1)*nc)**2)
          endif
       enddo ! c
    enddo ! r
    call interp_nldas2var(source,nc,nr,wind_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_WINDFORC,source,varfield,vlevel=1,   &
         units="m/s")

    ! EMK...Create total wind speed forcing in km/day
    do r = 1,nr
       do c = 1,nc
          if ((ewind_f(c+(r-1)*nc).gt.-9998).and.                      &
              (nwind_f(c+(r-1)*nc).gt.-9998)) then
             wind_f(c+(r-1)*nc) = sqrt(ewind_f(c+(r-1)*nc)**2 +        &
                                       nwind_f(c+(r-1)*nc)**2)         &
                                       * 86400.0 * 0.001
          endif
       enddo ! c
    enddo ! r
    call interp_nldas2var(source,nc,nr,wind_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_WINDFORC,source,varfield,  &
             vlevel=1,units="km/day")

    call interp_nldas2var(source,nc,nr,fheight_f,varfield)
    call LVT_logSingleDataStreamVar(LVT_MOC_FHEIGHTFORC,source,varfield,&
             vlevel=1,units="m")

    if (nldas2data(source)%lsm.eq."NOAH") then
       if ((nldas2data(source)%smvol).and.                             &
           ((LVT_MOC_SOILMOIST(source).ge.1).or.                       &
            (LVT_MOC_ROOTMOIST(source).ge.1))) then
          do r = 1,nr
             do c = 1,nc
                if (sml1(c+(r-1)*nc).gt.0) then
                   sml1(c+(r-1)*nc) = sml1(c+(r-1)*nc) /  100.0
                   sml2(c+(r-1)*nc) = sml2(c+(r-1)*nc) /  300.0
                   sml3(c+(r-1)*nc) = sml3(c+(r-1)*nc) /  600.0
                   sml4(c+(r-1)*nc) = sml4(c+(r-1)*nc) / 1000.0
                   rootmoist(c+(r-1)*nc) = rootmoist(c+(r-1)*nc) / 1000.0
                endif
             enddo
          enddo
       endif
    elseif (nldas2data(source)%lsm.eq."MOS") then
       if ((nldas2data(source)%smvol).and.                             &
           ((LVT_MOC_SOILMOIST(source).ge.1).or.                       &
            (LVT_MOC_ROOTMOIST(source).ge.1))) then
          do r = 1,nr
             do c = 1,nc
                if (sml1(c+(r-1)*nc).gt.0) then
                   sml1(c+(r-1)*nc) = sml1(c+(r-1)*nc) /  100.0
                   sml2(c+(r-1)*nc) = sml2(c+(r-1)*nc) /  300.0
                   sml3(c+(r-1)*nc) = sml3(c+(r-1)*nc) / 1600.0
                   rootmoist(c+(r-1)*nc) = rootmoist(c+(r-1)*nc) / 1000.0
                endif
             enddo
          enddo
       endif
    elseif (nldas2data(source)%lsm.eq."VIC") then
       if ((nldas2data(source)%smvol).and.                             &
           ((LVT_MOC_SOILMOIST(source).ge.1).or.                       &
            (LVT_MOC_ROOTMOIST(source).ge.1))) then
          do r = 1,nr
             do c = 1,nc
                if (sml1(c+(r-1)*nc).gt.0) then
                   sml1(c+(r-1)*nc) = sml1(c+(r-1)*nc) / (1000.0 *     &
                        nldas2data(source)%vic_depth1(c,r))
                   sml2(c+(r-1)*nc) = sml2(c+(r-1)*nc) / (1000.0 *     &
                        nldas2data(source)%vic_depth2(c,r))
                   sml3(c+(r-1)*nc) = sml3(c+(r-1)*nc) / (1000.0 *     &
                        nldas2data(source)%vic_depth3(c,r))
                   rootmoist(c+(r-1)*nc) = rootmoist(c+(r-1)*nc) / 1000.0
                endif
             enddo
          enddo
       endif
    endif

    if (nldas2data(source)%smvol) then
       call interp_nldas2var(source,nc,nr,sml1,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                vlevel=1,units="m3/m3")
       call interp_nldas2var(source,nc,nr,sml2,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                vlevel=2,units="m3/m3")
       call interp_nldas2var(source,nc,nr,sml3,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                vlevel=3,units="m3/m3")
       if ((nldas2data(source)%lsm.eq."NOAH").or.                      &
           (nldas2data(source)%lsm.eq."SAC")) then
          call interp_nldas2var(source,nc,nr,sml4,varfield)
          call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                   vlevel=4,units="m3/m3")
       endif
       call interp_nldas2var(source,nc,nr,rootmoist,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_ROOTMOIST,source,varfield,&
                vlevel=1,units="m3/m3")
    else
       call interp_nldas2var(source,nc,nr,sml1,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                vlevel=1,units="kg/m2")
       call interp_nldas2var(source,nc,nr,sml2,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                vlevel=2,units="kg/m2")
       call interp_nldas2var(source,nc,nr,sml3,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                vlevel=3,units="kg/m2")
       if ((nldas2data(source)%lsm.eq."NOAH").or.                      &
           (nldas2data(source)%lsm.eq."SAC")) then
          call interp_nldas2var(source,nc,nr,sml4,varfield)
          call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,&
                   vlevel=4,units="kg/m2")
       endif
    endif

    if (nldas2data(source)%lsm.ne."SAC") then
       call interp_nldas2var(source,nc,nr,stl1,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILTEMP,source,varfield,&
                vlevel=1,units="K")
    endif
    if ((nldas2data(source)%lsm.eq."NOAH").or.                         &
        (nldas2data(source)%lsm.eq."VIC")) then
       call interp_nldas2var(source,nc,nr,stl2,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILTEMP,source,varfield,&
                vlevel=2,units="K")
       call interp_nldas2var(source,nc,nr,stl3,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILTEMP,source,varfield,&
                vlevel=3,units="K")
    endif
    if (nldas2data(source)%lsm.eq."NOAH") then
       call interp_nldas2var(source,nc,nr,stl4,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_SOILTEMP,source,varfield,&
                vlevel=4,units="K")
    endif

    if (LVT_MOC_STREAMFLOW(source).ge.1) then
       call interp_nldas2var(source,nc,nr,streamflow,varfield)
       call LVT_logSingleDataStreamVar(LVT_MOC_STREAMFLOW,source,varfield,&
            vlevel=1,units="m3/s")
    endif

  end subroutine readNLDAS2data

!BOP
!
! !ROUTINE: interp_nldas2var
!  \label{interp_nldas2var}
!
! !INTERFACE:
  subroutine interp_nldas2var(source,nc,nr,var_input,var_output)
!
! !USES:
    use LVT_coreMod,    only : LVT_rc
    use NLDAS2_dataMod, only : nldas2data

    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the NLDAS-2 variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (NLDAS-2) grid
!    \item[nr]      number of rows in the input (NLDAS-2) grid
!    \item[var_input] input variable to be interpolated
!    \item[lb]        input bitmap (true//false)
!    \item[var_output] resulting interpolated field
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
! !ARGUMENTS:
    integer            :: source
    integer            :: nc
    integer            :: nr
    real               :: var_input(nc*nr)
    logical*1          :: lb(nc*nr)
    real               :: var_output(LVT_rc%lnc, LVT_rc%lnr)
!EOP
    integer            :: iret
    integer            :: c,r
    logical*1          :: lo(LVT_rc%lnc*LVT_rc%lnr)
    real               :: go(LVT_rc%lnc*LVT_rc%lnr)

    var_output = LVT_rc%udef
    lb = .false.
    do r = 1,nr
       do c = 1,nc
          if (var_input(c+(r-1)*nc).ne.LVT_rc%udef) then
             lb(c+(r-1)*nc) = .true.
          endif
       enddo
    enddo

    call bilinear_interp(LVT_rc%gridDesc,lb,var_input,                 &
         lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,                            &
         nldas2data(source)%rlat,nldas2data(source)%rlon,              &
         nldas2data(source)%w11,nldas2data(source)%w12,                &
         nldas2data(source)%w21,nldas2data(source)%w22,                &
         nldas2data(source)%n11,nldas2data(source)%n12,                &
         nldas2data(source)%n21,nldas2data(source)%n22,                &
         LVT_rc%udef,iret)

    do r = 1,LVT_rc%lnr
       do c = 1,LVT_rc%lnc
          var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)
       enddo
    enddo

  end subroutine interp_nldas2var
!BOP
!
! !ROUTINE: create_nldas2data_filename
! \label{create_nldas2data_filename}
!
! !INTERFACE:
  subroutine create_nldas2data_filename(odir,interval,                 &
                                        lsm,yr,mo,da,hr,doy,filename)

    implicit none
!
! !USES:
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine creates the NLDAS-2 filename based on the given date
!  (year, month, day, and hour).
!
!  The arguments are:
!  \begin{description}
!   \item[odir]      NLDAS-2 base directory
!   \item[interval]  NLDAS-2 hourly or monthly data
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[hr]        hour of data
!   \item[filename]  Name of the EBBR file
!  \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
!BOP
!
! !ARGUMENTS:
    character(len=*)  :: odir
    character(len=*)  :: interval
    character(len=*)  :: lsm
    integer           :: yr
    integer           :: mo
    integer           :: da
    integer           :: hr
    integer           :: doy
    character(len=*)  :: filename
!EOP
    character*4       :: fyr
    character*2       :: fmo
    character*2       :: fda
    character*2       :: fhr
    character*3       :: fdoy

    write(unit=fyr, fmt='(i4.4)') yr
    write(unit=fmo, fmt='(i2.2)') mo
    write(unit=fda, fmt='(i2.2)') da
    write(unit=fhr, fmt='(i2.2)') hr
    write(unit=fdoy,fmt='(i3.3)') doy

    if (trim(interval).eq."hourly") then
       if ((lsm.eq."NOAH").or.(lsm.eq."MOS").or.(lsm.eq."VIC")) then
          filename = trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//     &
               '/NLDAS_'//trim(lsm)//'0125_H.A'//trim(fyr)//           &
               trim(fmo)//trim(fda)//'.'//trim(fhr)//'00.002.grb'
       elseif (lsm.eq."SAC") then
          filename = trim(odir)//'/'//trim(fyr)//'/'//                 &
               trim(fyr)//trim(fmo)//trim(fda)//'/'//trim(fyr)//       &
               trim(fmo)//trim(fda)//trim(fhr)//'.SAC.grb'
       elseif (lsm.eq."SACSM") then
          filename = trim(odir)//'/SM/'//trim(fyr)//'/'//trim(fyr)//   &
               trim(fmo)//trim(fda)//trim(fhr)//'.SAC.gdat'
       elseif (lsm.eq."NOAHST") then
          filename = trim(odir)//'/streamflow/'//trim(fyr)//'/'//      &
               trim(fyr)//trim(fmo)//trim(fda)//'/'//trim(fyr)//       &
               trim(fmo)//trim(fda)//trim(fhr)//                       &
               '.NOAH.streamflow.bin'
       elseif (lsm.eq."MOSST") then
          filename = trim(odir)//'/streamflow/'//trim(fyr)//'/'//      &
               trim(fyr)//trim(fmo)//trim(fda)//'/'//trim(fyr)//       &
               trim(fmo)//trim(fda)//trim(fhr)//                       &
               '.Mosaic.streamflow.bin'
       elseif (lsm.eq."VICST") then
          filename = trim(odir)//'/streamflow/'//trim(fyr)//'/'//      &
               trim(fyr)//trim(fmo)//trim(fda)//'/'//trim(fyr)//       &
               trim(fmo)//trim(fda)//trim(fhr)//                       &
               '.VIC.streamflow.bin'
       elseif (lsm.eq."SACST") then
          filename = trim(odir)//'/streamflow/'//trim(fyr)//'/'//      &
               trim(fyr)//trim(fmo)//trim(fda)//'/'//trim(fyr)//       &
               trim(fmo)//trim(fda)//trim(fhr)//                       &
               '.SAC.streamflow.bin'
       elseif (lsm.eq."FORCINGA") then
          filename = trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//     &
               '/NLDAS_FORA0125_H.A'//trim(fyr)//                      &
               trim(fmo)//trim(fda)//'.'//trim(fhr)//'00.002.grb'
       elseif (lsm.eq."FORCINGB") then
          filename = trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//     &
               '/NLDAS_FORB0125_H.A'//trim(fyr)//                      &
               trim(fmo)//trim(fda)//'.'//trim(fhr)//'00.002.grb'
       endif

    elseif (trim(interval).eq."monthly") then
       if ((lsm.eq."NOAH").or.(lsm.eq."MOS").or.(lsm.eq."VIC")) then
          filename = trim(odir)//'/'//trim(fyr)//                      &
               '/NLDAS_'//trim(lsm)//'0125_M.A'//trim(fyr)//           &
               trim(fmo)//'01.'//trim(fhr)//'00.002.grb'
       elseif (lsm.eq."FORCINGA") then
          filename = trim(odir)//'/'//trim(fyr)//                      &
               '/NLDAS_FORA0125_H.A'//trim(fyr)//                      &
               trim(fmo)//'01.'//trim(fhr)//'00.002.grb'
       elseif (lsm.eq."FORCINGB") then
          filename = trim(odir)//'/'//trim(fyr)//                      &
               '/NLDAS_FORB0125_H.A'//trim(fyr)//                      &
               trim(fmo)//'01.'//trim(fhr)//'00.002.grb'
       endif
    endif

  end subroutine create_nldas2data_filename
!BOP
!
! !ROUTINE: create_nldas2data_filename
! \label{create_nldas2data_filename}
!
! !INTERFACE:
  subroutine retrieve_nldas2data(igrib,nc,nr,nvars,index,nldas_var)
! !USES:
    use grib_api
    use LVT_logMod, only : LVT_verify

    implicit none
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine sums and counts up for each NLDAS-2 variable.
!
! !REVISION HISTORY:
!
!EOP
!BOP
!
! !ARGUMENTS:
    integer              :: c,r,iret
    integer, intent(in)  :: igrib,nc,nr
    integer, intent(in)  :: nvars,index
    real                 :: var(nvars,nc*nr)
    real,    intent(out) :: nldas_var(nc*nr)

    call grib_get(igrib,"values",var(index,:),iret)
    call LVT_verify(iret,                                              &
         'grib_get failed for values in readNLDAS2data')

    do r = 1,nr
       do c = 1,nc
          if (var(index,c+(r-1)*nc).ne.9999.0) then
             nldas_var(c+(r-1)*nc) = var(index,c+(r-1)*nc)
          endif
       enddo
    enddo

  end subroutine retrieve_nldas2data
