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
! !ROUTINE: readGLDAS1Obs
! \label{readGLDAS1Obs}
!
! !INTERFACE: 
subroutine readGLDAS1Obs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use GLDAS1obsMod
  use grib_api

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This plugin processes the Global Land Data Assimilation System (GLDAS)
!   version 1 data available from NASA GES-DISC.
!   
!  NOTES: 
!  Currently the NOAH model-based monthly outputs in NetCDF format is 
!  is supported. The data can be downloaded from: 
!  http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP
  integer                 :: nc 
  integer                 :: nr 
  character*120           :: fname
  logical                 :: file_exists
  logical                 :: qle_flag, qh_flag, qg_flag
  logical                 :: precip_flag, tskin_flag
  integer                 :: ftn,iret
  integer                 :: c,r,t

  integer                 :: swnet_index,lwnet_index
  integer                 :: qle_index,qh_index,qg_index
  integer                 :: evap_index,qs_index
  integer                 :: rainf_index, snowf_index
  integer                 :: qsb_index,qsm_index,avgsurft_index,albedo_index
  integer                 :: sml1_index,sml2_index,sml3_index,sml4_index
  integer                 :: sml5_index,sml6_index,sml7_index,sml8_index
  integer                 :: sml9_index,sml10_index
  integer                 :: stl1_index,stl2_index,stl3_index,stl4_index
  integer                 :: stl5_index,stl6_index,stl7_index,stl8_index
  integer                 :: stl9_index,stl10_index
  integer                 :: potevap_index,ecanop_index,tveg_index,esoil_index
  integer                 :: rootmoist_index,canopint_index,subsnow_index
  integer                 :: snowcover_index,snowdepth_index,swe_index
    
  real*8                  :: timenow
  real                    :: gmt1
  type(ESMF_Time)         :: time1
  
  integer                 :: nvars,index,igrib,status
  integer, allocatable    :: pid(:),tid(:),topLev(:),botLev(:)
  real, allocatable       :: var(:,:)
  real                    :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: rnet(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qle_ip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qh_ip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qg_ip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: rainf_ip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: snowf_ip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: precip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qs_ip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qsb_ip(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: trunoff_ip(LVT_rc%lnc, LVT_rc%lnr)

  real, allocatable       :: swnet(:),lwnet(:)
  real, allocatable       :: qle(:),qh(:),qg(:)
  real, allocatable       :: evap(:),qs(:)
  real, allocatable       :: rainf(:), snowf(:)
  real, allocatable       :: qsb(:),qsm(:),avgsurft(:),albedo(:)
  real, allocatable       :: sml1(:),sml2(:),sml3(:),sml4(:)
  real, allocatable       :: sml5(:),sml6(:),sml7(:),sml8(:)
  real, allocatable       :: sml9(:),sml10(:)
  real, allocatable       :: stl1(:),stl2(:),stl3(:),stl4(:)
  real, allocatable       :: stl5(:),stl6(:),stl7(:),stl8(:)
  real, allocatable       :: stl9(:),stl10(:)
  real, allocatable       :: swe(:)

  integer                 :: sml1_topLev,sml2_topLev,sml3_topLev,sml4_topLev
  integer                 :: sml5_topLev,sml6_topLev,sml7_topLev,sml8_topLev
  integer                 :: sml9_topLev,sml10_topLev
  integer                 :: sml1_botLev,sml2_botLev,sml3_botLev,sml4_botLev
  integer                 :: sml5_botLev,sml6_botLev,sml7_botLev,sml8_botLev
  integer                 :: sml9_botLev,sml10_botLev
  integer                 :: stl1_topLev,stl2_topLev,stl3_topLev,stl4_topLev
  integer                 :: stl5_topLev,stl6_topLev,stl7_topLev,stl8_topLev
  integer                 :: stl9_topLev,stl10_topLev
  integer                 :: stl1_botLev,stl2_botLev,stl3_botLev,stl4_botLev
  integer                 :: stl5_botLev,stl6_botLev,stl7_botLev,stl8_botLev
  integer                 :: stl9_botLev,stl10_botLev
  logical                 :: alarmCheck

  nc = gldas1obs(source)%nc
  nr = gldas1obs(source)%nr

  swnet_index = 111
  lwnet_index = 112
  qle_index = 121
  qh_index = 122
  qg_index = 155
  evap_index = 57
  qs_index = 235
  qsb_index = 234
  qsm_index = 99
  avgsurft_index = 138
  swe_index = 65
  rainf_index = 132
  snowf_index = 131

  allocate(swnet(nc*nr))
  allocate(lwnet(nc*nr))
  allocate(qle(nc*nr))
  allocate(qh(nc*nr))
  allocate(qg(nc*nr))
  allocate(evap(nc*nr))
  allocate(qs(nc*nr))
  allocate(qsb(nc*nr))
  allocate(qsm(nc*nr))
  allocate(avgsurft(nc*nr))
  allocate(albedo(nc*nr))
  allocate(swe(nc*nr))
  allocate(rainf(nc*nr))
  allocate(snowf(nc*nr))

  swnet = LVT_rc%udef
  lwnet = LVT_rc%udef
  qle = LVT_rc%udef
  qh = LVT_rc%udef
  qg = LVT_rc%udef
  evap = LVT_rc%udef
  qs = LVT_rc%udef
  qsb = LVT_rc%udef
  qsm = LVT_rc%udef
  avgsurft = LVT_rc%udef
  swe = LVT_rc%udef
  rainf = LVT_rc%udef
  snowf = LVT_rc%udef

  qle_ip = LVT_rc%udef
  qh_ip = LVT_rc%udef
  qg_ip = LVT_rc%udef
  rnet = LVT_rc%udef
  rainf_ip  = LVT_rc%udef
  snowf_ip  = LVT_rc%udef
  precip  = LVT_rc%udef
  qs_ip = LVT_rc%udef
  qsb_ip = LVT_rc%udef
  trunoff_ip = LVT_rc%udef

  if(GLDAS1obs(source)%model_name.eq."NOAH") then 

     stl1_index = 85
     stl2_index = 85
     stl3_index = 85
     stl4_index = 85
     sml1_index = 86
     sml2_index = 86
     sml3_index = 86
     sml4_index = 86

     sml1_topLev = 0
     sml1_botLev = 4
     sml2_topLev = 0
     sml2_botLev = 3
     sml3_topLev = 0
     sml3_botLev = 2
     sml4_topLev = 0
     sml4_botLev = 1
     stl1_topLev = 0
     stl1_botLev = 4
     stl2_topLev = 0
     stl2_botLev = 3
     stl3_topLev = 0
     stl3_botLev = 2
     stl4_topLev = 0
     stl4_botLev = 1

     allocate(sml1(nc*nr))
     allocate(sml2(nc*nr))
     allocate(sml3(nc*nr))
     allocate(sml4(nc*nr))

     allocate(stl1(nc*nr))
     allocate(stl2(nc*nr))
     allocate(stl3(nc*nr))
     allocate(stl4(nc*nr))

     sml1 = LVT_rc%udef
     sml2 = LVT_rc%udef
     sml3 = LVT_rc%udef
     sml4 = LVT_rc%udef
     stl1 = LVT_rc%udef
     stl2 = LVT_rc%udef
     stl3 = LVT_rc%udef
     stl4 = LVT_rc%udef
  elseif(GLDAS1obs(source)%model_name.eq."CLM") then 
     stl1_index = 85
     stl2_index = 85
     stl3_index = 85
     stl4_index = 85
     stl5_index = 85
     stl6_index = 85
     stl7_index = 85
     stl8_index = 85
     stl9_index = 85
     stl10_index = 85

     sml1_index = 86
     sml2_index = 86
     sml3_index = 86
     sml4_index = 86
     sml5_index = 86
     sml6_index = 86
     sml7_index = 86
     sml8_index = 86
     sml9_index = 86
     sml10_index = 86

     sml1_topLev = 0
     sml1_botLev = 110
     sml2_topLev = 0
     sml2_botLev = 109
     sml3_topLev = 0
     sml3_botLev = 108
     sml4_topLev = 0
     sml4_botLev = 107
     sml5_topLev = 0
     sml5_botLev = 106
     sml6_topLev = 0
     sml6_botLev = 105
     sml7_topLev = 0
     sml7_botLev = 104
     sml8_topLev = 0
     sml8_botLev = 103
     sml9_topLev = 0
     sml9_botLev = 102
     sml10_topLev = 0
     sml10_botLev = 101

     stl1_topLev = 0
     stl1_botLev = 110
     stl2_topLev = 0
     stl2_botLev = 109
     stl3_topLev = 0
     stl3_botLev = 108
     stl4_topLev = 0
     stl4_botLev = 107
     stl5_topLev = 0
     stl5_botLev = 106
     stl6_topLev = 0
     stl6_botLev = 105
     stl7_topLev = 0
     stl7_botLev = 104
     stl8_topLev = 0
     stl8_botLev = 103
     stl9_topLev = 0
     stl9_botLev = 102
     stl10_topLev = 0
     stl10_botLev = 101

     allocate(sml1(nc*nr))
     allocate(sml2(nc*nr))
     allocate(sml3(nc*nr))
     allocate(sml4(nc*nr))
     allocate(sml5(nc*nr))
     allocate(sml6(nc*nr))
     allocate(sml7(nc*nr))
     allocate(sml8(nc*nr))
     allocate(sml9(nc*nr))
     allocate(sml10(nc*nr))

     allocate(stl1(nc*nr))
     allocate(stl2(nc*nr))
     allocate(stl3(nc*nr))
     allocate(stl4(nc*nr))
     allocate(stl5(nc*nr))
     allocate(stl6(nc*nr))
     allocate(stl7(nc*nr))
     allocate(stl8(nc*nr))
     allocate(stl9(nc*nr))
     allocate(stl10(nc*nr))


     sml1 = LVT_rc%udef
     sml2 = LVT_rc%udef
     sml3 = LVT_rc%udef
     sml4 = LVT_rc%udef
     sml5 = LVT_rc%udef
     sml6 = LVT_rc%udef
     sml7 = LVT_rc%udef
     sml8 = LVT_rc%udef
     sml9 = LVT_rc%udef
     sml10 = LVT_rc%udef

     stl1 = LVT_rc%udef
     stl2 = LVT_rc%udef
     stl3 = LVT_rc%udef
     stl4 = LVT_rc%udef
     stl5 = LVT_rc%udef
     stl6 = LVT_rc%udef
     stl7 = LVT_rc%udef
     stl8 = LVT_rc%udef
     stl9 = LVT_rc%udef
     stl10 = LVT_rc%udef


  elseif(GLDAS1obs(source)%model_name.eq."MOS") then 
     stl1_index = 85

     sml1_index = 86
     sml2_index = 86
     sml3_index = 86

     stl1_topLev = 195
     stl1_botLev = 81

     sml1_topLev = 195
     sml1_botLev = 81
     sml2_topLev = 156
     sml2_botLev = 66
     sml3_topLev = 117
     sml3_botLev = 51

     allocate(sml1(nc*nr))
     allocate(sml2(nc*nr))
     allocate(sml3(nc*nr))

     allocate(stl1(nc*nr))

     sml1 = LVT_rc%udef
     sml2 = LVT_rc%udef
     sml3 = LVT_rc%udef

     stl1 = LVT_rc%udef
  elseif(GLDAS1obs(source)%model_name.eq."VIC") then 

     sml1_index = 86
     sml2_index = 86
     sml3_index = 86

     sml1_topLev = 0 
     sml1_botLev = 10 
     sml2_topLev = 10
     sml2_botLev = 160
     sml3_topLev = 160
     sml3_botLev = 190

     allocate(sml1(nc*nr))
     allocate(sml2(nc*nr))
     allocate(sml3(nc*nr))

     sml1 = LVT_rc%udef
     sml2 = LVT_rc%udef
     sml3 = LVT_rc%udef

  endif

  if(LVT_MOC_QLE(source).ge.1.or.&
       LVT_MOC_RNET(source).ge.1) then 
     qle_flag = .true. 
  else
     qle_flag = .false. 
  endif
  if(LVT_MOC_QH(source).ge.1.or.&
       LVT_MOC_RNET(source).ge.1) then 
     qh_flag = .true. 
  else
     qh_flag = .false. 
  endif

  if(LVT_MOC_QG(source).ge.1.or.&
       LVT_MOC_RNET(source).ge.1) then 
     qg_flag = .true. 
  else
     qg_flag = .false. 
  endif

  if(LVT_MOC_AVGSURFT(source).ge.1.or.&
       LVT_MOC_RADT(source).ge.1) then 
     tskin_flag = .true. 
  else
     tskin_flag = .false. 
  endif

  if(LVT_MOC_TOTALPRECIP(source).ge.1.or.&
       LVT_MOC_RAINF(source).ge.1.or.&
       LVT_MOC_RAINFFORC(source).ge.1) then 
     precip_flag = .true. 
  else
     precip_flag = .false. 
  endif
  
  timenow = float(LVT_rc%dhr(source))*3600+ &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmCheck = (mod(timenow,10800.0).eq.0)

  if(GLDAS1obs(source)%startFlag.or.alarmCheck) then

     
     if(GLDAS1obs(source)%startFlag) then 
        GLDAS1obs(source)%startFlag = .false. 
     endif

     call create_GLDAS1_filename(GLDAS1obs(source)%odir,&
          GLDAS1obs(source)%model_name, &
          GLDAS1obs(source)%datares, &
          LVT_rc%dyr(source),&
          LVT_rc%dmo(source),&
          LVT_rc%ddoy(source),& 
          LVT_rc%dhr(source),& 
          fname)
     
     inquire(file=trim(fname),exist=file_exists) 

     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading GLDAS1 file ',trim(fname)
        
        call grib_open_file(ftn,trim(fname),'r',iret)
        if(iret.ne.0) then 
           write(LVT_logunit,*) & 
                '[ERR] File not opened ',trim(fname)
           call LVT_endrun()
        endif
        
        call grib_count_in_file(ftn,nvars)

        allocate(var(nvars,nc*nr))
        allocate(pid(nvars))
        allocate(tid(nvars))
        allocate(topLev(nvars))
        allocate(botLev(nvars))

        do index=1,nvars
           call grib_new_from_file(ftn,igrib,iret)
           call LVT_verify(iret,&
                'grib_new_from_file failed in readGLDAS1obs')
           
           call grib_get(igrib,"indicatorOfParameter",pid(index),iret)
           call LVT_verify(iret,&
                'grib_get failed for indicatorOfParameter in readGLDAS1obs')
           
           call grib_get(igrib, "timeRangeIndicator",tid(index), iret)
           call LVT_verify(iret, &
                'grib_get failed for timeRangeIndicator in readGLDAS1obs')

           call grib_get(igrib,"topLevel", topLev(index), iret)
           call LVT_verify(iret, &
                'grib_get failed for topLevel in readGLDAS1obs')

           call grib_get(igrib,"bottomLevel", botLev(index), iret)
           call LVT_verify(iret, &
                'grib_get failed for bottomLevel in readGLDAS1obs')

           if((pid(index).eq.swnet_index).and.(LVT_MOC_SWNET(source).ge.1)) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,swnet)
           endif

           if((pid(index).eq.lwnet_index).and.(LVT_MOC_LWNET(source).ge.1)) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,lwnet)
           endif

           if((pid(index).eq.qle_index).and.qle_flag) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,qle)
           endif
           
           if((pid(index).eq.qh_index).and.qh_flag) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,qh)
           endif

           if((pid(index).eq.qg_index).and.qg_flag) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,qg)
           endif
           
           if((pid(index).eq.evap_index).and.(LVT_MOC_EVAP(source).ge.1)) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,evap)
           endif

           if((pid(index).eq.qs_index).and.((LVT_MOC_QS(source).ge.1).or.&
                (LVT_MOC_RUNOFF(source).ge.1))) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,qs)
           endif

           if((pid(index).eq.qsb_index).and.((LVT_MOC_QSB(source).ge.1).or.&
                (LVT_MOC_RUNOFF(source).ge.1))) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,qsb)
           endif

           if((pid(index).eq.qsm_index).and.(LVT_MOC_QSM(source).ge.1)) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,qsm)
           endif

           if((pid(index).eq.rainf_index).and.precip_flag) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,rainf)
           endif
           if((pid(index).eq.snowf_index).and.precip_flag) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,snowf)
           endif
           
           if(GLDAS1obs(source)%model_name.ne."VIC") then 
              if((pid(index).eq.avgsurft_index).and.tskin_flag) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,avgsurft)
              endif
           endif

           if((pid(index).eq.swe_index).and.(LVT_MOC_SWE(source).ge.1)) then 
              call retrieve_gldas1data(igrib,nc,nr,nvars,index,swe)
           endif

           if(GLDAS1obs(source)%model_name.eq."NOAH") then 
              if((pid(index).eq.sml1_index).and.       &
                   (topLev(index).eq.sml1_topLev).and. &
                   (botLev(index).eq.sml1_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml1)
              endif
              
              if((pid(index).eq.sml2_index).and.       &
                   (topLev(index).eq.sml2_topLev).and. &
                   (botLev(index).eq.sml2_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml2)
              endif
              
              if((pid(index).eq.sml3_index).and.       &
                   (topLev(index).eq.sml3_topLev).and. &
                   (botLev(index).eq.sml3_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml3)
              endif
              
              if((pid(index).eq.sml4_index).and.       &
                   (topLev(index).eq.sml4_topLev).and. &
                   (botLev(index).eq.sml4_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml4)
              endif
              

              if((pid(index).eq.stl1_index).and.       &
                   (topLev(index).eq.stl1_topLev).and. &
                   (botLev(index).eq.stl1_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl1)
              endif
              
              if((pid(index).eq.stl2_index).and.       &
                   (topLev(index).eq.stl2_topLev).and. &
                   (botLev(index).eq.stl2_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl2)
              endif
              
              if((pid(index).eq.stl3_index).and.       &
                   (topLev(index).eq.stl3_topLev).and. &
                   (botLev(index).eq.stl3_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl3)
              endif
              
              if((pid(index).eq.stl4_index).and.       &
                   (topLev(index).eq.stl4_topLev).and. &
                   (botLev(index).eq.stl4_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl4)
              endif
           elseif(GLDAS1obs(source)%model_name.eq."CLM") then 
              if((pid(index).eq.sml1_index).and.       &
                   (topLev(index).eq.sml1_topLev).and. &
                   (botLev(index).eq.sml1_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml1)
              endif
              
              if((pid(index).eq.sml2_index).and.       &
                   (topLev(index).eq.sml2_topLev).and. &
                   (botLev(index).eq.sml2_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml2)
              endif
              
              if((pid(index).eq.sml3_index).and.       &
                   (topLev(index).eq.sml3_topLev).and. &
                   (botLev(index).eq.sml3_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml3)
              endif
              
              if((pid(index).eq.sml4_index).and.       &
                   (topLev(index).eq.sml4_topLev).and. &
                   (botLev(index).eq.sml4_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml4)
              endif
              
              if((pid(index).eq.sml5_index).and.       &
                   (topLev(index).eq.sml5_topLev).and. &
                   (botLev(index).eq.sml5_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml5)
              endif

              if((pid(index).eq.sml6_index).and.       &
                   (topLev(index).eq.sml6_topLev).and. &
                   (botLev(index).eq.sml6_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml6)
              endif

              if((pid(index).eq.sml7_index).and.       &
                   (topLev(index).eq.sml7_topLev).and. &
                   (botLev(index).eq.sml7_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml7)
              endif
              
              if((pid(index).eq.sml8_index).and.       &
                   (topLev(index).eq.sml8_topLev).and. &
                   (botLev(index).eq.sml8_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml8)
              endif
              
              if((pid(index).eq.sml9_index).and.       &
                   (topLev(index).eq.sml9_topLev).and. &
                   (botLev(index).eq.sml9_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml9)
              endif
              if((pid(index).eq.sml10_index).and.       &
                   (topLev(index).eq.sml10_topLev).and. &
                   (botLev(index).eq.sml10_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml10)
              endif

              if((pid(index).eq.stl1_index).and.       &
                   (topLev(index).eq.stl1_topLev).and. &
                   (botLev(index).eq.stl1_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl1)
              endif
              
              if((pid(index).eq.stl2_index).and.       &
                   (topLev(index).eq.stl2_topLev).and. &
                   (botLev(index).eq.stl2_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl2)
              endif
              
              if((pid(index).eq.stl3_index).and.       &
                   (topLev(index).eq.stl3_topLev).and. &
                   (botLev(index).eq.stl3_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl3)
              endif
              
              if((pid(index).eq.stl4_index).and.       &
                   (topLev(index).eq.stl4_topLev).and. &
                   (botLev(index).eq.stl4_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl4)
              endif
              
              if((pid(index).eq.stl5_index).and.       &
                   (topLev(index).eq.stl5_topLev).and. &
                   (botLev(index).eq.stl5_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl5)
              endif

              if((pid(index).eq.stl6_index).and.       &
                   (topLev(index).eq.stl6_topLev).and. &
                   (botLev(index).eq.stl6_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl6)
              endif

              if((pid(index).eq.stl7_index).and.       &
                   (topLev(index).eq.stl7_topLev).and. &
                   (botLev(index).eq.stl7_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl7)
              endif
              
              if((pid(index).eq.stl8_index).and.       &
                   (topLev(index).eq.stl8_topLev).and. &
                   (botLev(index).eq.stl8_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl8)
              endif
              
              if((pid(index).eq.stl9_index).and.       &
                   (topLev(index).eq.stl9_topLev).and. &
                   (botLev(index).eq.stl9_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl9)
              endif
              if((pid(index).eq.stl10_index).and.       &
                   (topLev(index).eq.stl10_topLev).and. &
                   (botLev(index).eq.stl10_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl10)
              endif
           elseif(GLDAS1obs(source)%model_name.eq."MOS") then 
              if((pid(index).eq.sml1_index).and.       &
                   (topLev(index).eq.sml1_topLev).and. &
                   (botLev(index).eq.sml1_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml1)
              endif
              
              if((pid(index).eq.sml2_index).and.       &
                   (topLev(index).eq.sml2_topLev).and. &
                   (botLev(index).eq.sml2_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml2)
              endif
              
              if((pid(index).eq.sml3_index).and.       &
                   (topLev(index).eq.sml3_topLev).and. &
                   (botLev(index).eq.sml3_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml3)
              endif
              
              if((pid(index).eq.stl1_index).and.       &
                   (topLev(index).eq.stl1_topLev).and. &
                   (botLev(index).eq.stl1_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,stl1)
              endif
           elseif(GLDAS1obs(source)%model_name.eq."VIC") then 
              if((pid(index).eq.sml1_index).and.       &
                   (topLev(index).eq.sml1_topLev).and. &
                   (botLev(index).eq.sml1_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml1)
              endif
              
              if((pid(index).eq.sml2_index).and.       &
                   (topLev(index).eq.sml2_topLev).and. &
                   (botLev(index).eq.sml2_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml2)
              endif
              
              if((pid(index).eq.sml3_index).and.       &
                   (topLev(index).eq.sml3_topLev).and. &
                   (botLev(index).eq.sml3_botLev).and. &
                   (LVT_MOC_SOILMOIST(source).ge.1)) then 
                 call retrieve_gldas1data(igrib,nc,nr,nvars,index,sml3)
              endif
              
           endif

           call grib_release(igrib,iret)
           call LVT_verify(iret, &
                'grib_release failed in readGLDAS1obs')
           
        enddo
        call grib_close_file(ftn,iret)
        
        deallocate(var)
        deallocate(pid)
        deallocate(tid)
        deallocate(topLev)
        deallocate(botLev)
     endif

  end if

  call interp_gldas1var(source,nc,nr,swnet,varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_SWNET, source, varfield, &
       vlevel=1,units="W/m2")


  call interp_gldas1var(source,nc,nr,lwnet,varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_LWNET, source, varfield, &
       vlevel=1,units="W/m2")

  call interp_gldas1var(source,nc,nr,qle,qle_ip)
  call LVT_logSingleDataStreamVar(LVT_MOC_QLE, source, qle_ip, &
       vlevel=1,units="W/m2")

  call interp_gldas1var(source,nc,nr,qh,qh_ip)
  call LVT_logSingleDataStreamVar(LVT_MOC_QH, source, qh_ip, &
       vlevel=1,units="W/m2")
  
  call interp_gldas1var(source,nc,nr,qg,qg_ip)
  call LVT_logSingleDataStreamVar(LVT_MOC_QG, source, qg_ip, &
       vlevel=1,units="W/m2")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qle_ip(c,r).ne.LVT_rc%udef) then 
           rnet(c,r) = qle_ip(c,r)+qh_ip(c,r)+qg_ip(c,r)
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_RNET, source, rnet, &
       vlevel=1,units="W/m2")

  call interp_gldas1var(source,nc,nr,evap,varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_EVAP, source, varfield, &
       vlevel=1,units="kg/m2")

  call interp_gldas1var(source,nc,nr,qs,qs_ip)
  call LVT_logSingleDataStreamVar(LVT_MOC_QS, source, qs_ip, &
       vlevel=1,units="kg/m2")

  call interp_gldas1var(source,nc,nr,qsb,qsb_ip)
  call LVT_logSingleDataStreamVar(LVT_MOC_QSB, source, qsb_ip, &
       vlevel=1,units="kg/m2")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qs_ip(c,r).ne.LVT_rc%udef) then 
           trunoff_ip(c,r) = qs_ip(c,r)+qsb_ip(c,r)
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_RUNOFF, source, trunoff_ip, &
       vlevel=1,units="kg/m2")

  call interp_gldas1var(source,nc,nr,qsm,varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_QSM, source, varfield, &
       vlevel=1,units="kg/m2")

  call interp_gldas1var(source,nc,nr,avgsurft,varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT, source, varfield, &
       vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_RADT, source, varfield, &
       vlevel=1,units="K")

  call interp_gldas1var(source,nc,nr,swe,varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_SWE, source, varfield, &
       vlevel=1,units="kg/m2")

  call interp_gldas1var(source,nc,nr,rainf,rainf_ip)
  call interp_gldas1var(source,nc,nr,snowf,snowf_ip)
  
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(rainf_ip(c,r).ne.LVT_rc%udef) then 
           precip(c,r) = rainf_ip(c,r)+snowf_ip(c,r)
        endif
     enddo
  enddo


  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP, source, precip, &
       vlevel=1,units="kg/m2s")

  if((GLDAS1obs(source)%model_name.eq."NOAH").and.&
       LVT_MOC_SOILMOIST(source).ge.1) then 
     do r=1,nr
        do c=1,nc
            if (sml1(c+(r-1)*nc).gt.0) then
               sml1(c+(r-1)*nc) = sml1(c+(r-1)*nc) /  100.0
               sml2(c+(r-1)*nc) = sml2(c+(r-1)*nc) /  300.0
               sml3(c+(r-1)*nc) = sml3(c+(r-1)*nc) /  600.0
               sml4(c+(r-1)*nc) = sml4(c+(r-1)*nc) / 1000.0
            endif
         enddo
      enddo
      call interp_gldas1var(source,nc,nr, sml1, varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=1, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml2,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=2, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml3,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=3, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml4,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=4, &
           units="m3/m3")

      call interp_gldas1var(source,nc,nr, stl1, varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=1, &
           units="K")
      
      call interp_gldas1var(source,nc,nr, stl2,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=2, &
           units="K")
      
      call interp_gldas1var(source,nc,nr, stl3,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=3, &
           units="K")
      
      call interp_gldas1var(source,nc,nr, stl4,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=4, &
           units="K")
      
   elseif((GLDAS1obs(source)%model_name.eq."CLM").and.&
       LVT_MOC_SOILMOIST(source).ge.1) then 
      do r=1,nr
         do c=1,nc
            if (sml1(c+(r-1)*nc).gt.0) then
               sml1(c+(r-1)*nc) = sml1(c+(r-1)*nc) / 18.0
               sml2(c+(r-1)*nc) = sml2(c+(r-1)*nc) / 27.0
               sml3(c+(r-1)*nc) = sml3(c+(r-1)*nc) / 46.0
               sml4(c+(r-1)*nc) = sml4(c+(r-1)*nc) / 75.0
               sml5(c+(r-1)*nc) = sml5(c+(r-1)*nc) / 123.0
               sml6(c+(r-1)*nc) = sml6(c+(r-1)*nc) / 204.0
               sml7(c+(r-1)*nc) = sml7(c+(r-1)*nc) / 336.0
               sml8(c+(r-1)*nc) = sml8(c+(r-1)*nc) / 554.0
               sml9(c+(r-1)*nc) = sml9(c+(r-1)*nc) / 913.0
               sml10(c+(r-1)*nc) = sml10(c+(r-1)*nc) / 1137.0
            endif
         enddo
      enddo
      
      call interp_gldas1var(source,nc,nr, sml1, varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=1, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml2,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=2, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml3,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=3, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml4,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=4, &
           units="m3/m3")

      call interp_gldas1var(source,nc,nr, sml5, varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=5, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml6,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=6, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml7,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=7, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml8,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=8, &
           units="m3/m3")

      call interp_gldas1var(source,nc,nr, sml9,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=9, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml10,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=10, &
           units="m3/m3")

   elseif((GLDAS1obs(source)%model_name.eq."MOS").and.&
       LVT_MOC_SOILMOIST(source).ge.1) then 
     do r=1,nr
        do c=1,nc
            if (sml1(c+(r-1)*nc).gt.0) then
               sml1(c+(r-1)*nc) = sml1(c+(r-1)*nc) /  20.0
               sml2(c+(r-1)*nc) = sml2(c+(r-1)*nc) / 1480.0 
               sml3(c+(r-1)*nc) = sml3(c+(r-1)*nc) / 2000.0 
            endif
         enddo
      enddo
      call interp_gldas1var(source,nc,nr, sml1, varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=1, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml2,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=2, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml3,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=3, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, stl1, varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=1, &
           units="K")

   elseif((GLDAS1obs(source)%model_name.eq."VIC").and.&
       LVT_MOC_SOILMOIST(source).ge.1) then 
     do r=1,nr
        do c=1,nc
            if (sml1(c+(r-1)*nc).gt.0) then
               sml1(c+(r-1)*nc) = sml1(c+(r-1)*nc) /  100.0
               sml2(c+(r-1)*nc) = sml2(c+(r-1)*nc) / 1500.0
               sml3(c+(r-1)*nc) = sml3(c+(r-1)*nc) /  300.0
            endif
         enddo
      enddo
      call interp_gldas1var(source,nc,nr, sml1, varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=1, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml2,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=2, &
           units="m3/m3")
      
      call interp_gldas1var(source,nc,nr, sml3,  varfield)
      call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,varfield,vlevel=3, &
           units="m3/m3")
   endif

   deallocate(swnet)
   deallocate(lwnet)
   deallocate(qle)
   deallocate(qh)
   deallocate(qg)
   deallocate(evap)
   deallocate(qs)
   deallocate(qsb)
   deallocate(qsm)
   deallocate(avgsurft)
   deallocate(albedo)
   deallocate(swe)

   if(GLDAS1obs(source)%model_name.eq."NOAH") then 
      deallocate(sml1)
      deallocate(sml2)
      deallocate(sml3)
      deallocate(sml4)
      
      deallocate(stl1)
      deallocate(stl2)
      deallocate(stl3)
      deallocate(stl4)

   elseif(GLDAS1obs(source)%model_name.eq."CLM") then 
        
      deallocate(sml1)
      deallocate(sml2)
      deallocate(sml3)
      deallocate(sml4)
      deallocate(sml5)
      deallocate(sml6)
      deallocate(sml7)
      deallocate(sml8)
      deallocate(sml9)
      deallocate(sml10)
      
      deallocate(stl1)
      deallocate(stl2)
      deallocate(stl3)
      deallocate(stl4)
      deallocate(stl5)
      deallocate(stl6)
      deallocate(stl7)
      deallocate(stl8)
      deallocate(stl9)
      deallocate(stl10)
      
   elseif(GLDAS1obs(source)%model_name.eq."MOS") then 
      deallocate(sml1)
      deallocate(sml2)
      deallocate(sml3)

      deallocate(stl1)

   elseif(GLDAS1obs(source)%model_name.eq."VIC") then 
      deallocate(sml1)
      deallocate(sml2)
      deallocate(sml3)

   endif
   
 end subroutine readGLDAS1Obs
 
!BOP
!
! !ROUTINE: retrieve_gldas1data
! \label{retrieve_gldas1data}
!
! !INTERFACE:
  subroutine retrieve_gldas1data(igrib,nc,nr,nvars,index,nldas_var)
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
    call LVT_verify(iret,                                            &
         'grib_get failed for values in readGLDAS1data')

    do r = 1,nr
       do c = 1,nc
          if (var(index,c+(r-1)*nc).ne.9999.0) then
             nldas_var(c+(r-1)*nc) = var(index,c+(r-1)*nc)
          endif
       enddo
    enddo

  end subroutine retrieve_gldas1data



!BOP
!
! !ROUTINE: interp_gldas1var
!  \label{interp_gldas1var}
!
! !INTERFACE:
  subroutine interp_gldas1var(source, nc,nr,var_input,var_output)
!
! !USES:
    use LVT_coreMod
    use GLDAS1obsMod
      
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
!EOP
!BOP
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

    if(LVT_isAtAfinerResolution(gldas1obs(source)%datares)) then
       call neighbor_interp(LVT_rc%gridDesc,lb,var_input,  &
            lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,             &
            gldas1obs(source)%rlat,gldas1obs(source)%rlon, &
            gldas1obs(source)%n11,                         & 
            LVT_rc%udef,iret)
    else
       call upscaleByAveraging(&
            nc*nr, &
            LVT_rc%lnc*LVT_rc%lnr, &
            LVT_rc%udef, &
            gldas1obs(source)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LVT_rc%lnr
       do c = 1,LVT_rc%lnc
          var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)
       enddo
    enddo

  end subroutine interp_gldas1var


!BOP
! 
! !ROUTINE: create_GLDAS1_filename
! \label{create_GLDAS1_filename}
!
! !INTERFACE: 
subroutine create_GLDAS1_filename(odir, model_name, datares, &
     yr,mo,doy, hr, filename)
! 
! !USES:   
  use LVT_logMod

  implicit none
  !
! !ARGUMENTS: 
  character(len=*)             :: odir
  real                         :: datares
  character(len=*)             :: model_name
  integer                      :: yr
  integer                      :: mo
  integer                      :: doy
  integer                      :: hr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the GLDAS1 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            GLDAS1 base directory
!   \item[model\_name]     name of the model used in the GLDAS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the GLDAS1 file
!  \end{description}
! 
!EOP
  
  integer                 :: ftn
  integer                 :: ierr
  character*4             :: fyr
  character*3             :: fdoy
  character*2             :: fmo
  character*2             :: fhr

  character*100           :: list_name

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fhr, fmt='(i2.2)') hr

  if(datares .eq. 0.25) then   
     list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
          '/GLDAS_'//trim(model_name)//'025SUBP_3H.A'//&
          trim(fyr)//trim(fdoy)//'.'//trim(fhr)//&
          '*grb > GLDAS1_file'
  elseif(datares.eq. 1.0) then 
     if(model_name.ne."VIC") then 
        list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
             '/GLDAS_'//trim(model_name)//'10SUBP_3H.A'//&
             trim(fyr)//trim(fdoy)//'.'//trim(fhr)//&
             '*grb > GLDAS1_file'
     else
        list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
             '/GLDAS_'//trim(model_name)//'10_3H.A'//&
             trim(fyr)//trim(fdoy)//'.'//trim(fhr)//&
             '*grb > GLDAS1_file'
     endif
  endif
  call system(trim(list_name))         
     
  ftn = LVT_getNextUnitNumber()
  open(ftn,file='GLDAS1_file',status='old',iostat=ierr)
  do while(ierr.eq.0) 
     read(ftn,'(a)',iostat=ierr) filename
     if(ierr.ne.0) then 
        exit
     endif
  enddo
  call LVT_releaseUnitNumber(ftn)

end subroutine create_GLDAS1_filename


