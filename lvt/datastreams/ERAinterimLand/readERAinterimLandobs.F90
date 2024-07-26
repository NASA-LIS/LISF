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
! !ROUTINE: readERAinterimLandObs
! \label{readERAinterimLandObs}
!
! !INTERFACE: 
subroutine readERAinterimLandObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use ERAinterimLandobsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This plugin processes the ERA Interim Land data 
!   data available from http://apps.ecmwf.int/datasets/
!   
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP

  real                    :: timenow
  logical                 :: alarmCheck, alarmCheck1
  integer                 :: c,r, k,nc,nr
  integer                 :: yr, mo, da, hr, mn, ss, doy
  real                    :: gmt
  integer                 :: st,et,t
  real*8                  :: lis_prevtime
  type(ESMF_Time)         :: merra2time1, merra2time2, initTime
  type(ESMF_TimeInterval) :: dayInterval
  character(len=100)      :: fcst1_filename, fcst2_filename, anlys_filename

  real                    :: sm1(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: sm2(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: sm3(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: sm4(LVT_rc%lnc, LVT_rc%lnr)

  real                    :: st1(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: st2(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: st3(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: st4(LVT_rc%lnc, LVT_rc%lnr)

  real                    :: sd(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: swe(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: skt(LVT_rc%lnc, LVT_rc%lnr)

  real                    :: qle(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qh(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qsm(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: evap(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: tp(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: rnet(LVT_rc%lnc, LVT_rc%lnr)

  real                    :: srunoff(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: baseflow(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: trunoff(LVT_rc%lnc,LVT_rc%lnr)

  real                    :: depth(4), sf_wt(4), rz_wt(4)
  real                    :: sfsm(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: rzsm(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: sfst(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: rzst(LVT_rc%lnc, LVT_rc%lnr)
  integer                 :: status

  nc = ERAIlandobs(source)%nc
  nr = ERAIlandobs(source)%nr

  sm1 = LVT_rc%udef
  sm2 = LVT_rc%udef
  sm3 = LVT_rc%udef
  sm4 = LVT_rc%udef

  st1 = LVT_rc%udef
  st2 = LVT_rc%udef
  st3 = LVT_rc%udef
  st4 = LVT_rc%udef
  
  sd = LVT_rc%udef
  swe = LVT_rc%udef
  skt = LVT_rc%udef

  srunoff = LVT_rc%udef
  baseflow = LVT_rc%udef
  trunoff = LVT_rc%udef

  qle = LVT_rc%udef
  qh  = LVT_rc%udef
  tp  = LVT_rc%udef
  qsm = LVT_rc%udef
  evap = LVT_rc%udef

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 10800.0).eq.0)
         
  if(alarmCheck.or.(LVT_rc%dda(source).ne.&
       ERAIlandobs(source)%da).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     ERAIlandobs(source)%da = LVT_rc%dda(source)

     call create_ERAinterimLand_anlys_filename(&
          ERAIlandobs(source)%odir,LVT_rc%dyr(source),anlys_filename)
     
     call process_ERAIlandANLYSdata(source, anlys_filename, &
          sm1, sm2, sm3, sm4, &
          st1, st2, st3, st4, &
          sd, swe, skt)

     call create_ERAinterimLand_fcst1_filename(&
          ERAIlandobs(source)%odir,LVT_rc%dyr(source),fcst1_filename)

     call process_ERAIland_fcst1_data(source, fcst1_filename, &
          qle, qh, rnet, tp, qsm, evap)
        
     timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
          + LVT_rc%dss(source)
     alarmcheck1 = (mod(timenow, 86400.0).eq.0)
     
     if(alarmCheck1) then 
        call create_ERAinterimLand_fcst2_filename(&
             ERAIlandobs(source)%odir,LVT_rc%dyr(source),fcst2_filename)
        
        call process_ERAIland_runoff_data(source, fcst2_filename, &
             srunoff, baseflow)
        
     endif
     
  endif

  depth(1) = 0.07
  depth(2) = 0.28
  depth(3) = 1.00
  depth(4) = 2.89
  
  call compute_vinterp_weights(4, LVT_rc%lis_sf_d, &
       LVT_rc%lis_rz_d, & 
       depth(1:4), sf_wt, rz_wt)

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(sm1(c,r).gt.0) then 
           rzsm(c,r) = rz_wt(1)*sm1(c,r) + & 
                rz_wt(2)*sm2(c,r) + & 
                rz_wt(3)*sm3(c,r) + & 
                rz_wt(4)*sm4(c,r) 
           
           sfsm(c,r) = sf_wt(1)*sm1(c,r) + & 
                sf_wt(2)*sm2(c,r) + & 
                sf_wt(3)*sm3(c,r) + & 
                sf_wt(4)*sm4(c,r) 

           rzst(c,r) = rz_wt(1)*st1(c,r) + & 
                rz_wt(2)*st2(c,r) + & 
                rz_wt(3)*st3(c,r) + & 
                rz_wt(4)*st4(c,r) 
           
           sfst(c,r) = sf_wt(1)*st1(c,r) + & 
                sf_wt(2)*st2(c,r) + & 
                sf_wt(3)*st3(c,r) + & 
                sf_wt(4)*st4(c,r) 

        else
           sfsm(c,r) = LVT_rc%udef
           rzsm(c,r) = LVT_rc%udef
           sfst(c,r) = LVT_rc%udef
           rzst(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  
  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, &
       source, sfsm,&
       vlevel=1,units="m3/m3")

  call LVT_logSingleDataStreamVar(LVT_MOC_ROOTMOIST, &
       source, rzsm,&
       vlevel=1,units="m3/m3")

  call LVT_logSingleDataStreamVar(LVT_MOC_SOILTEMP, &
       source, sfst,&
       vlevel=1,units="K")

  call LVT_logSingleDataStreamVar(LVT_MOC_ROOTTEMP, &
       source, rzst,&
       vlevel=1,units="K")

  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH, &
       source, sd,&
       vlevel=1,units="m")

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE, &
       source, swe,&
       vlevel=1,units="kg/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT, &
       source, skt,&
       vlevel=1,units="K")

  call LVT_logSingleDataStreamVar(LVT_MOC_RADT, &
       source, skt,&
       vlevel=1,units="K")

  call LVT_logSingleDataStreamVar(LVT_MOC_QS, &
       source, srunoff,&
       vlevel=1,units="kg/m2s")

  call LVT_logSingleDataStreamVar(LVT_MOC_QSB, &
       source, baseflow,&
       vlevel=1,units="kg/m2s")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(srunoff(c,r).ne.LVT_rc%udef) then 
           trunoff(c,r) = srunoff(c,r) + baseflow(c,r)
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_RUNOFF, &
       source, trunoff,&
       vlevel=1,units="kg/m2s")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qle(c,r).ne.LVT_rc%udef) then 
           qle(c,r) = -1.0*qle(c,r)
        endif
        if(qh(c,r).ne.LVT_rc%udef) then 
           qh(c,r) = -1.0*qh(c,r)
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE, &
       source, qle,&
       vlevel=1,units="W/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_QH, &
       source, qh,&
       vlevel=1,units="W/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_RNET, &
       source, rnet,&
       vlevel=1,units="W/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP, &
       source, tp,&
       vlevel=1,units="kg/m2s")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(tp(c,r).ne.LVT_rc%udef) then 
           tp(c,r) = tp(c,r)*3600.0
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP, &
       source, tp,&
       vlevel=1,units="kg/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_QSM, &
       source, qsm,&
       vlevel=1,units="kg/m2s")

  call LVT_logSingleDataStreamVar(LVT_MOC_EVAP, &
       source, evap,&
       vlevel=1,units="kg/m2s")

end subroutine readERAinterimLandObs

subroutine process_ERAIlandANLYSdata(source, anlysfile, &
     sm1, sm2, sm3, sm4, st1, st2, st3, st4, sd, swe, skt)

  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use ERAinterimLandobsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  
  integer                :: source
  character(len=*)       :: anlysfile

  integer                :: yr
  integer                :: mo
  integer                :: da


  real                    :: sm1(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: sm2(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: sm3(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: sm4(LVT_rc%lnc, LVT_rc%lnr)

  real                    :: st1(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: st2(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: st3(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: st4(LVT_rc%lnc, LVT_rc%lnr)

  real                    :: sd(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: swe(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: skt(LVT_rc%lnc, LVT_rc%lnr)

  real                    :: rsn(LVT_rc%lnc, LVT_rc%lnr)

  integer            :: ftn
  character*100      :: fname
  logical            :: file_exists
  integer            :: timeid
  integer            :: swvl1id,swvl2id,swvl3id,swvl4id
  integer            :: stl1id, stl2id, stl3id, stl4id
  integer            :: sdid, rsnid, sktid

  real  , pointer    :: swvl1(:,:,:)
  real  , pointer    :: swvl2(:,:,:)
  real  , pointer    :: swvl3(:,:,:)
  real  , pointer    :: swvl4(:,:,:)

  real  , pointer    :: stl1(:,:,:)
  real  , pointer    :: stl2(:,:,:)
  real  , pointer    :: stl3(:,:,:)
  real  , pointer    :: stl4(:,:,:)

  real  , pointer    :: sd_in(:,:,:)
  real  , pointer    :: rsn_in(:,:,:)
  real  , pointer    :: skt_in(:,:,:)

  integer            :: c,r,k,iret,status
  real                    :: scale_sm_f,offset_sm_v
  real                    :: scale_st_f,offset_st_v
  real                    :: scale_sd_f,offset_sd_v
  real                    :: scale_rsn_f,offset_rsn_v
  real                    :: scale_skt_f,offset_skt_v

  type(ESMF_TimeInterval) :: deltaT
  type(ESMF_Time)         :: currTime, startYearTime
  integer                 :: time(1)
  integer                 :: hr_elapsed, tindex
 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=trim(anlysfile),exist=file_exists) 

  if(file_exists) then   

     allocate(swvl1(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(swvl2(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(swvl3(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(swvl4(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))

     allocate(stl1(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(stl2(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(stl3(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(stl4(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))

     allocate(sd_in(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(rsn_in(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))

     allocate(skt_in(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     
     
     iret = nf90_open(path=trim(anlysfile),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 

        call ESMF_TimeSet(currTime,yy=LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h = LVT_rc%dhr(source), &
             m = LVT_rc%dmn(source), &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (1): readERAinterimLandobs')

        call ESMF_TimeSet(startYearTime,yy=LVT_rc%dyr(source), &
             mm = 1, &
             dd = 1, &
             h = 0, &
             m = 0, &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (2): readERAinterimLandobs')
        
        deltaT = startYearTime - ERAIlandobs(source)%startTime
        call ESMF_TimeIntervalGet(deltaT, h=hr_elapsed)
        
        call LVT_verify(nf90_inq_varid(ftn,"time",timeid),&
             'nf90_inq_varid failed for time')
        
        call LVT_verify(nf90_get_var(ftn,timeid,time,&
             start=(/1/),count=(/1/)),&
             'Error in nf90_get_var time')

        if(hr_elapsed.ne.time(1)) then

           write(LVT_logunit,*) '[ERR] The starting time in the ERA interim anlys file'
           write(LVT_logunit,*) '[ERR] do not match the computed starting time'
           call LVT_endrun()
        endif

        deltaT = currTime - startYearTime
        call ESMF_TimeIntervalGet(deltaT, h=tindex)
        
        if(mod(tindex,6).eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading ERAinterimLand file ',trim(anlysfile)
           tindex = tindex/6 + 1
           
           call LVT_verify(nf90_inq_varid(ftn,"swvl1",swvl1id),&
                'nf90_inq_varid failed for swvl1')
           call LVT_verify(nf90_inq_varid(ftn,"swvl2",swvl2id),&
                'nf90_inq_varid failed for swvl2')
           call LVT_verify(nf90_inq_varid(ftn,"swvl3",swvl3id),&
                'nf90_inq_varid failed for swvl3')
           call LVT_verify(nf90_inq_varid(ftn,"swvl4",swvl4id),&
                'nf90_inq_varid failed for swvl4')           

           call LVT_verify(nf90_inq_varid(ftn,"stl1",stl1id),&
                'nf90_inq_varid failed for stl1')
           call LVT_verify(nf90_inq_varid(ftn,"stl2",stl2id),&
                'nf90_inq_varid failed for stl2')
           call LVT_verify(nf90_inq_varid(ftn,"stl3",stl3id),&
                'nf90_inq_varid failed for stl3')
           call LVT_verify(nf90_inq_varid(ftn,"stl4",stl4id),&
                'nf90_inq_varid failed for stl4')

           call LVT_verify(nf90_inq_varid(ftn,"sd",sdid),&
                'nf90_inq_varid failed for sd')           
           call LVT_verify(nf90_inq_varid(ftn,"rsn",rsnid),&
                'nf90_inq_varid failed for rsn')           
           call LVT_verify(nf90_inq_varid(ftn,"skt",sktid),&
                'nf90_inq_varid failed for skt')           

           call LVT_verify(nf90_get_var(ftn,swvl1id,swvl1,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var swvl1')
           call LVT_verify(nf90_get_var(ftn,swvl2id,swvl2,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var swvl2')
           call LVT_verify(nf90_get_var(ftn,swvl2id,swvl3,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var swvl3')
           call LVT_verify(nf90_get_var(ftn,swvl4id,swvl4,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var swvl4')           
           call LVT_verify(nf90_get_att(ftn,swvl1id, "scale_factor",scale_sm_f),&
                'Error in nf90_get_att scale_factor')          
           call LVT_verify(nf90_get_att(ftn,swvl1id, "add_offset",offset_sm_v),&
                'Error in nf90_get_att add_offset')


           call LVT_verify(nf90_get_var(ftn,stl1id,stl1,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var stl1')
           call LVT_verify(nf90_get_var(ftn,stl2id,stl2,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var stl2')
           call LVT_verify(nf90_get_var(ftn,stl3id,stl3,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var stl2')
           call LVT_verify(nf90_get_var(ftn,stl4id,stl4,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var stl4')
           call LVT_verify(nf90_get_att(ftn,stl1id, "scale_factor",scale_st_f),&
                'Error in nf90_get_att scale_factor')          
           call LVT_verify(nf90_get_att(ftn,stl1id, "add_offset",offset_st_v),&
                'Error in nf90_get_att add_offset')           

           call LVT_verify(nf90_get_var(ftn,sdid,sd_in,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var sd')           
           call LVT_verify(nf90_get_att(ftn,sdid, "scale_factor",scale_sd_f),&
                'Error in nf90_get_att scale_factor')          
           call LVT_verify(nf90_get_att(ftn,sdid, "add_offset",offset_sd_v),&
                'Error in nf90_get_att add_offset')

           call LVT_verify(nf90_get_var(ftn,rsnid,rsn_in,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var rsn')           
           call LVT_verify(nf90_get_att(ftn,rsnid, "scale_factor",scale_rsn_f),&
                'Error in nf90_get_att scale_factor')          
           call LVT_verify(nf90_get_att(ftn,rsnid, "add_offset",offset_rsn_v),&
                'Error in nf90_get_att add_offset')

           call LVT_verify(nf90_get_var(ftn,sktid,skt_in,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var skt')           
           call LVT_verify(nf90_get_att(ftn,sktid, "scale_factor",scale_skt_f),&
                'Error in nf90_get_att scale_factor')          
           call LVT_verify(nf90_get_att(ftn,sktid, "add_offset",offset_skt_v),&
                'Error in nf90_get_att add_offset')

           
           call interp_ERAIlandvar2d(source,scale_sm_f, offset_sm_v, swvl1,&
                sm1)
           call interp_ERAIlandvar2d(source,scale_sm_f, offset_sm_v, swvl2,&
                sm2)
           call interp_ERAIlandvar2d(source,scale_sm_f, offset_sm_v, swvl3,&
                sm3)
           call interp_ERAIlandvar2d(source,scale_sm_f, offset_sm_v, swvl4,&
                sm4)

           call interp_ERAIlandvar2d(source,scale_st_f, offset_st_v, stl1,&
                st1)
           call interp_ERAIlandvar2d(source,scale_st_f, offset_st_v, stl2,&
                st2)
           call interp_ERAIlandvar2d(source,scale_st_f, offset_st_v, stl3,&
                st3)
           call interp_ERAIlandvar2d(source,scale_st_f, offset_st_v, stl4,&
                st4)

           call interp_ERAIlandvar2d(source,scale_sd_f, offset_sd_v, sd_in,&
                sd)
           call interp_ERAIlandvar2d(source,scale_rsn_f, offset_rsn_v, rsn_in,&
                rsn)
           call interp_ERAIlandvar2d(source,scale_skt_f, offset_skt_v, skt_in,&
                skt)
           
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(sd(c,r).ne.LVT_rc%udef) then
                    swe(c,r) = sd(c,r)*rsn(c,r)
                 endif
              enddo
           enddo
        endif
     endif
     call LVT_verify(nf90_close(ftn))

     deallocate(swvl1)
     deallocate(swvl2)
     deallocate(swvl3)
     deallocate(swvl4)

     deallocate(stl1)
     deallocate(stl2)
     deallocate(stl3)
     deallocate(stl4)

     deallocate(sd_in)
     deallocate(rsn_in)
     deallocate(skt_in)
  end if
#endif
end subroutine process_ERAIlandANLYSdata


subroutine process_ERAIland_runoff_data(source, fcstfile, &
     srunoff, baseflow)

  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use LVT_constantsMod
  use ERAinterimLandobsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  
  integer                :: source
  character(len=*)       :: fcstfile

  integer                :: yr
  integer                :: mo
  integer                :: da


  real                    :: srunoff(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: baseflow(LVT_rc%lnc, LVT_rc%lnr)

  integer            :: ftn
  character*100      :: fname
  logical            :: file_exists
  integer            :: timeid
  integer            :: sroid, ssroid

  real  , pointer    :: sro(:,:,:)
  real  , pointer    :: ssro(:,:,:)

  integer            :: c,r,k,iret,status
  real                    :: scale_sro_f,offset_sro_v
  real                    :: scale_ssro_f,offset_ssro_v

  type(ESMF_TimeInterval) :: deltaT
  type(ESMF_Time)         :: currTime, startYearTime
  integer                 :: time(1)
  integer                 :: hr_elapsed, tindex
 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=trim(fcstfile),exist=file_exists) 

  if(file_exists) then   

     allocate(sro(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))
     allocate(ssro(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1))

     iret = nf90_open(path=trim(fcstfile),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 

        call ESMF_TimeSet(currTime,yy=LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h = LVT_rc%dhr(source), &
             m = LVT_rc%dmn(source), &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (3): readERAinterimLandobs')

        call ESMF_TimeSet(startYearTime,yy=LVT_rc%dyr(source), &
             mm = 1, &
             dd = 2, &
             h = 0, &
             m = 0, &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (4): readERAinterimLandobs')
        
        deltaT = startYearTime - ERAIlandobs(source)%startTime
        call ESMF_TimeIntervalGet(deltaT, h=hr_elapsed)
        
        call LVT_verify(nf90_inq_varid(ftn,"time",timeid),&
             'nf90_inq_varid failed for time')
        
        call LVT_verify(nf90_get_var(ftn,timeid,time,&
             start=(/1/),count=(/1/)),&
             'Error in nf90_get_var time')

        if(hr_elapsed.ne.time(1)) then

           write(LVT_logunit,*) '[ERR] The starting time in the ERA interim fcst file'
           write(LVT_logunit,*) '[ERR] do not match the computed starting time'
           call LVT_endrun()
        endif

        call ESMF_TimeSet(startYearTime,yy=LVT_rc%dyr(source), &
             mm = 1, &
             dd = 1, &
             h = 0, &
             m = 0, &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (5): readERAinterimLandobs')

        deltaT = currTime - startYearTime
        call ESMF_TimeIntervalGet(deltaT, h=tindex)
        
        if(mod(tindex,24).eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading ERAinterimLand file ',trim(fcstfile)
     

           tindex = tindex/24 + 1
           
           call LVT_verify(nf90_inq_varid(ftn,"sro",sroid),&
                'nf90_inq_varid failed for sro')
           call LVT_verify(nf90_inq_varid(ftn,"ssro",ssroid),&
                'nf90_inq_varid failed for ssro')
          
           call LVT_verify(nf90_get_var(ftn,sroid,sro,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var sro')
           call LVT_verify(nf90_get_var(ftn,ssroid,ssro,&
                start=(/1,1,tindex/), &
                count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,1/)),&
                'Error in nf90_get_var ssro')

           call LVT_verify(nf90_get_att(ftn,sroid, "scale_factor",scale_sro_f),&
                'Error in nf90_get_att scale_factor')          
           call LVT_verify(nf90_get_att(ftn,sroid, "add_offset",offset_sro_v),&
                'Error in nf90_get_att add_offset')

           call LVT_verify(nf90_get_att(ftn,ssroid, "scale_factor",scale_ssro_f),&
                'Error in nf90_get_att scale_factor')          
           call LVT_verify(nf90_get_att(ftn,ssroid, "add_offset",offset_ssro_v),&
                'Error in nf90_get_att add_offset')           

           call interp_ERAIlandvar2d(source,scale_sro_f, offset_sro_v, sro,&
                srunoff)
           call interp_ERAIlandvar2d(source,scale_ssro_f, offset_ssro_v, ssro,&
                baseflow)

        endif
        call LVT_verify(nf90_close(ftn))           

     endif
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(srunoff(c,r).ne.LVT_rc%udef) then
              srunoff(c,r) = srunoff(c,r)*LVT_CONST_RHOFW/(24.0*3600.0)   !to kg/m2s
           endif
           if(baseflow(c,r).ne.LVT_rc%udef) then 
              baseflow(c,r) = baseflow(c,r)*LVT_CONST_RHOFW/(24.0*3600.0)   !to kg/m2s
           endif
        enddo
     enddo
    
     deallocate(sro)
     deallocate(ssro)
  end if
#endif
end subroutine process_ERAIland_runoff_data

subroutine process_ERAIland_fcst1_data(source, fcstfile, &
     qle_out, qh_out, rnet_out, tp_out, qsm_out, evap_out)

  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use LVT_constantsMod
  use ERAinterimLandobsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  
  integer                :: source
  character(len=*)       :: fcstfile

  integer                :: yr
  integer                :: mo
  integer                :: da


  real                    :: qle_out(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qh_out(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: rnet_out(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: tp_out(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: qsm_out(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: ssr_out(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: str_out(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: evap_out(LVT_rc%lnc, LVT_rc%lnr)

  integer            :: ftn
  character*100      :: fname
  logical            :: file_exists
  integer            :: timeid
  integer            :: qleid, qhid, tpid, qsmid,ssrid, strid, evapid

  real  , pointer    :: qle(:,:,:)
  real  , pointer    :: qh(:,:,:)
  real  , pointer    :: tp(:,:,:)
  real  , pointer    :: qsm(:,:,:)
  real  , pointer    :: evap(:,:,:)
  real  , pointer    :: ssr(:,:,:)
  real  , pointer    :: str(:,:,:)

  integer            :: c,r,k,iret,status
  real                    :: scale_qle_f,offset_qle_v
  real                    :: scale_qh_f,offset_qh_v
  real                    :: scale_tp_f,offset_tp_v
  real                    :: scale_qsm_f,offset_qsm_v
  real                    :: scale_evap_f,offset_evap_v
  real                    :: scale_ssr_f,offset_ssr_v
  real                    :: scale_str_f,offset_str_v

  type(ESMF_TimeInterval) :: deltaT
  type(ESMF_Time)         :: currTime, startTime, startYearTime
  integer                 :: time(1)
  real                    :: dt
  integer                 :: ncount, stid
  integer                 :: hr_elapsed, tindex, dindex
 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ssr_out = LVT_rc%udef
  str_out = LVT_rc%udef

  inquire(file=trim(fcstfile),exist=file_exists) 

  if(file_exists) then   

     allocate(qle(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,2))
     allocate(qh(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,2))
     allocate(tp(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,2))
     allocate(qsm(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,2))
     allocate(evap(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,2))

     allocate(ssr(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,2))
     allocate(str(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,2))

     
     iret = nf90_open(path=trim(fcstfile),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 

        call ESMF_TimeSet(currTime,yy=LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h = LVT_rc%dhr(source), &
             m = LVT_rc%dmn(source), &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (6): readERAinterimLandobs')

        call ESMF_TimeSet(startYearTime,yy=LVT_rc%dyr(source), &
             mm = 1, &
             dd = 1, &
             h = 3, &
             m = 0, &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (7): readERAinterimLandobs')
        
        deltaT = startYearTime - ERAIlandobs(source)%startTime
        call ESMF_TimeIntervalGet(deltaT, h=hr_elapsed)
        
        call LVT_verify(nf90_inq_varid(ftn,"time",timeid),&
             'nf90_inq_varid failed for time')
        
        call LVT_verify(nf90_get_var(ftn,timeid,time,&
             start=(/1/),count=(/1/)),&
             'Error in nf90_get_var time')

        if(hr_elapsed.ne.time(1)) then
           print*, 'h_fcst1', hr_elapsed, time(1)
           write(LVT_logunit,*) '[ERR] The starting time in the ERA interim fcst file'
           write(LVT_logunit,*) '[ERR] do not match the computed starting time'
           call LVT_endrun()
        endif

        call ESMF_TimeSet(startTime,yy=LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h = 0, &
             m = 0, &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (8): readERAinterimLandobs')

        call ESMF_TimeSet(startYearTime,yy=LVT_rc%dyr(source), &
             mm = 1, &
             dd = 1, &
             h = 0, &
             m = 0, &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'Error in ESMF_TimeSet (9): readERAinterimLandobs')

        deltaT = currTime - startTime
        call ESMF_TimeIntervalGet(deltaT, h=tindex)

        deltaT = currTime - startYearTime
        call ESMF_TimeIntervalGet(deltaT, d=dindex)

        if(mod(tindex,3).eq.0) then 
           write(LVT_logunit,*) '[INFO] Reading ERAinterimLand file ',trim(fcstfile)

           if(LVT_rc%dhr(source).eq.0) then 
              tindex = dindex*6
              dt = 6*3600.0
           elseif(LVT_rc%dhr(source).le.12) then 
              tindex = tindex/3 + dindex*6
              dt = 3*3600.0
           elseif(LVT_rc%dhr(source).eq.15) then 
              tindex = -1
              dt = -1
           elseif(LVT_rc%dhr(source).eq.18) then 
              tindex = tindex/3-1 + dindex*6
              dt = 6*3600.0
           elseif(LVT_rc%dhr(source).eq.21) then 
              tindex = -1
              dt = -1
           endif
              
           if(tindex.gt.0) then 
              if(tindex.gt.1) then 
                 stid = tindex - 1
                 ncount = 2
              else
                 stid = tindex
                 ncount = 1
              endif
              if(LVT_rc%dhr(source).eq.3) then 
                 stid = tindex
                 ncount = 1
              endif

              call LVT_verify(nf90_inq_varid(ftn,"slhf",qleid),&
                   'nf90_inq_varid failed for qle')
              call LVT_verify(nf90_inq_varid(ftn,"sshf",qhid),&
                   'nf90_inq_varid failed for sshf')
              call LVT_verify(nf90_inq_varid(ftn,"tp",tpid),&
                   'nf90_inq_varid failed for tp')
              call LVT_verify(nf90_inq_varid(ftn,"smlt",qsmid),&
                   'nf90_inq_varid failed for smlt')
              call LVT_verify(nf90_inq_varid(ftn,"e",evapid),&
                   'nf90_inq_varid failed for evap')
              call LVT_verify(nf90_inq_varid(ftn,"ssr",ssrid),&
                   'nf90_inq_varid failed for ssr')
              call LVT_verify(nf90_inq_varid(ftn,"str",strid),&
                   'nf90_inq_varid failed for str')
              
              if(ncount.eq.2) then 
                 call LVT_verify(nf90_get_var(ftn,qleid,qle(:,:,:),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var qle')
                 call LVT_verify(nf90_get_var(ftn,qhid,qh(:,:,:),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var qh')
                 call LVT_verify(nf90_get_var(ftn,tpid,tp(:,:,:),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                   'Error in nf90_get_var tp')
                 call LVT_verify(nf90_get_var(ftn, qsmid,qsm(:,:,:),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var smlt')
                 call LVT_verify(nf90_get_var(ftn, evapid,evap(:,:,:),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var evap')

                 call LVT_verify(nf90_get_var(ftn, ssrid,ssr(:,:,:),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var ssr')

                 call LVT_verify(nf90_get_var(ftn, strid,str(:,:,:),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var str')

                 do r=1,ERAIlandobs(source)%nr
                    do c=1,ERAIlandobs(source)%nc
                       if(qle(c,r,2).ne.-32767) then 
                          qle(c,r,1) = qle(c,r,2) - qle(c,r,1)
                       endif
                       if(qh(c,r,2).ne.-32767) then 
                          qh(c,r,1) = qh(c,r,2) - qh(c,r,1)
                       endif
                       if(tp(c,r,2).ne.-32767) then 
                          tp(c,r,1) = tp(c,r,2) - tp(c,r,1)
                       endif
                       if(qsm(c,r,2).ne.-32767) then 
                          qsm(c,r,1) = qsm(c,r,2) - qsm(c,r,1)
                       endif
                       if(evap(c,r,2).ne.-32767) then 
                          evap(c,r,1) = evap(c,r,2) - evap(c,r,1)
                       endif
                       if(ssr(c,r,2).ne.-32767) then 
                          ssr(c,r,1) = ssr(c,r,2) - ssr(c,r,1)
                       endif
                       if(str(c,r,2).ne.-32767) then 
                          str(c,r,1) = str(c,r,2) - str(c,r,1)
                       endif
                    enddo
                 enddo                 
                 
              else
                 call LVT_verify(nf90_get_var(ftn,qleid,qle(:,:,1),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var qle')
                 call LVT_verify(nf90_get_var(ftn,qhid,qh(:,:,1),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var qh')
                 call LVT_verify(nf90_get_var(ftn,tpid,tp(:,:,1),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                   'Error in nf90_get_var tp')
                 call LVT_verify(nf90_get_var(ftn, qsmid,qsm(:,:,1),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var smlt')
                 call LVT_verify(nf90_get_var(ftn, evapid,evap(:,:,1),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var evap')
                 call LVT_verify(nf90_get_var(ftn, ssrid,ssr(:,:,1),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var ssr')

                 call LVT_verify(nf90_get_var(ftn, strid,str(:,:,1),&
                      start=(/1,1,stid/), &
                      count=(/ERAIlandobs(source)%nc,ERAIlandobs(source)%nr,ncount/)),&
                      'Error in nf90_get_var str')
              endif

              call LVT_verify(nf90_get_att(ftn,qleid, "scale_factor",scale_qle_f),&
                   'Error in nf90_get_att scale_factor')          
              call LVT_verify(nf90_get_att(ftn,qleid, "add_offset",offset_qle_v),&
                   'Error in nf90_get_att add_offset')
              
              call LVT_verify(nf90_get_att(ftn,qhid, "scale_factor",scale_qh_f),&
                   'Error in nf90_get_att scale_factor')          
              call LVT_verify(nf90_get_att(ftn,qhid, "add_offset",offset_qh_v),&
                   'Error in nf90_get_att add_offset')           
              
              call LVT_verify(nf90_get_att(ftn,tpid, "scale_factor",scale_tp_f),&
                   'Error in nf90_get_att scale_factor')          
              call LVT_verify(nf90_get_att(ftn,tpid, "add_offset",offset_tp_v),&
                   'Error in nf90_get_att add_offset')           
              
              call LVT_verify(nf90_get_att(ftn,qsmid, "scale_factor",scale_qsm_f),&
                   'Error in nf90_get_att scale_factor')          
              call LVT_verify(nf90_get_att(ftn,qsmid, "add_offset",offset_qsm_v),&
                   'Error in nf90_get_att add_offset')           

              call LVT_verify(nf90_get_att(ftn,evapid, "scale_factor",scale_evap_f),&
                   'Error in nf90_get_att scale_factor')          
              call LVT_verify(nf90_get_att(ftn,evapid, "add_offset",offset_evap_v),&
                   'Error in nf90_get_att add_offset')           

              call LVT_verify(nf90_get_att(ftn,ssrid, "scale_factor",scale_ssr_f),&
                   'Error in nf90_get_att scale_factor')          
              call LVT_verify(nf90_get_att(ftn,ssrid, "add_offset",offset_ssr_v),&
                   'Error in nf90_get_att add_offset')           

              call LVT_verify(nf90_get_att(ftn,strid, "scale_factor",scale_str_f),&
                   'Error in nf90_get_att scale_factor')          
              call LVT_verify(nf90_get_att(ftn,strid, "add_offset",offset_str_v),&
                   'Error in nf90_get_att add_offset')           
              
              
              if(ncount.eq.2) then 
                 offset_qle_v = 0.0
                 offset_qh_v  = 0.0
                 offset_tp_v  = 0.0
                 offset_qsm_v = 0.0
                 offset_evap_v = 0.0
                 offset_ssr_v = 0.0
                 offset_str_v = 0.0
              endif
               
              
              call interp_ERAIlandvar2d(source,scale_qle_f,&
                   offset_qle_v, qle(:,:,1),&
                   qle_out)
              call interp_ERAIlandvar2d(source,scale_qh_f, &
                   offset_qh_v, qh(:,:,1),&
                   qh_out)
              call interp_ERAIlandvar2d(source,scale_tp_f, &
                   offset_tp_v, tp(:,:,1),&
                   tp_out)
              call interp_ERAIlandvar2d(source,scale_qsm_f, &
                   offset_qsm_v, qsm(:,:,1),&
                   qsm_out)

              call interp_ERAIlandvar2d(source,scale_evap_f, &
                   offset_evap_v, evap(:,:,1),&
                   evap_out)

              call interp_ERAIlandvar2d(source,scale_ssr_f, &
                   offset_ssr_v, ssr(:,:,1),&
                   ssr_out)
              call interp_ERAIlandvar2d(source,scale_str_f, &
                   offset_str_v, str(:,:,1),&
                   str_out)
              
           endif
        endif

        call LVT_verify(nf90_close(ftn))
              
     endif
     
     deallocate(qle)
     deallocate(qh)
     deallocate(tp)
     deallocate(qsm)
     deallocate(evap)
     deallocate(ssr)
     deallocate(str)

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(qle_out(c,r).ne.LVT_rc%udef) then
              qle_out(c,r) = qle_out(c,r)/dt
           endif
           if(qh_out(c,r).ne.LVT_rc%udef) then
              qh_out(c,r) = qh_out(c,r)/dt
           endif
           if(tp_out(c,r).ne.LVT_rc%udef) then
              tp_out(c,r) = tp_out(c,r)*LVT_CONST_RHOFW/dt
           endif
           if(qsm_out(c,r).ne.LVT_rc%udef) then
              qsm_out(c,r) = qsm_out(c,r)*LVT_CONST_RHOFW/dt
           endif
           if(evap_out(c,r).ne.LVT_rc%udef) then
              evap_out(c,r) = evap_out(c,r)*LVT_CONST_RHOFW/dt
           endif
           if(ssr_out(c,r).ne.LVT_rc%udef.and.&
                str_out(c,r).ne.LVT_rc%udef) then
              rnet_out(c,r) = (ssr_out(c,r) -str_out(c,r))/dt
           endif
        enddo
     enddo

  end if
#endif
end subroutine process_ERAIland_fcst1_data


!BOP
!
! !ROUTINE: interp_ERAIlandvar2d
! \label{interp_ERAIlandvar2d}
!
! !INTERFACE: 
subroutine interp_ERAIlandvar2d(source,scale_f,offset_v,var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use ERAinterimLandobsMod
! !ARGUMENTS: 
  integer      :: source
  real         :: scale_f
  real         :: offset_v
  real         :: var_inp(ERAIlandobs(source)%nc,ERAIlandobs(source)%nr)
  real         :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the MERRA fields to the 
!  target LVT domain
!
!EOP

  real          :: var_inp_1d(ERAIlandobs(source)%nc*ERAIlandobs(source)%nr)
  logical*1     :: input_bitmap(ERAIlandobs(source)%nc*ERAIlandobs(source)%nr)
  real          :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1     :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer       :: nc, nr, c,c1,r,k
  integer       :: iret


  nc = ERAIlandobs(source)%nc
  nr = ERAIlandobs(source)%nr

  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        if(c.ge.1.and.c.le.240) then 
           c1 = c+240
        else
           c1 = c-240
        endif
        if(var_inp(c,r).ne.-32767) then            
           var_inp_1d(c1+((nr-r+1)-1)*nc) = var_inp(c,r)*scale_f + offset_v
           input_bitmap(c1+((nr-r+1)-1)*nc) = .true. 
        else
           var_inp_1d(c1+((nr-r+1)-1)*nc) = LVT_rc%udef
        endif
     enddo
  enddo

  if(LVT_isAtAfinerResolution(ERAIlandobs(source)%datares)) then
     call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d, &
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          ERAIlandobs(source)%rlat, & 
          ERAIlandobs(source)%rlon, &
          ERAIlandobs(source)%n11, &
          LVT_rc%udef, iret)
     
  else
     call upscaleByAveraging(&
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          ERAIlandobs(source)%n11, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d)
     
  endif
  
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
           var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
        else
           var_out(c,r) = LVT_rc%udef
        endif

     enddo
  enddo

end subroutine interp_ERAIlandvar2d


!BOP
! 
! !ROUTINE: tavg_single_ERAinterimLand_var
! \label{tavg_single_ERAinterimLand_var}
!
! !INTERFACE: 
subroutine tavg_single_ERAinterimLand_var(var,nvar)
! 
! !USES: 
  use LVT_coreMod, only : LVT_rc

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine temporally averages a single variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  integer       :: source
  real          :: var(LVT_rc%lnc,LVT_rc%lnr)
  integer       :: nvar(LVT_rc%lnc,LVT_rc%lnr)
!EOP  
  
  integer       :: c,r

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nvar(c,r).gt.0) then 
           var(c,r) = var(c,r)/nvar(c,r)
        else
           var(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

end subroutine tavg_single_ERAinterimLand_var



!BOP
! 
! !ROUTINE: create_ERAinterimLand_fcst2_filename
! \label{create_ERAinterimLand_fcst2_filename}
!
! !INTERFACE: 
subroutine create_ERAinterimLand_fcst2_filename(odir,yr,filename)
! 
! !USES:   
  use LVT_logMod

  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the ERAinterimLand data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            ERAinterimLand base directory
!   \item[yr]              year of data
!   \item[filename]        Name of the ERAinterimLand file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/era-interim_stream_fcst2_'//trim(fyr)//'.nc'
  
end subroutine create_ERAinterimLand_fcst2_filename

!BOP
! 
! !ROUTINE: create_ERAinterimLand_fcst1_filename
! \label{create_ERAinterimLand_fcst1_filename}
!
! !INTERFACE: 
subroutine create_ERAinterimLand_fcst1_filename(odir,yr,filename)
! 
! !USES:   
  use LVT_logMod

  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the ERAinterimLand data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            ERAinterimLand base directory
!   \item[yr]              year of data
!   \item[filename]        Name of the ERAinterimLand file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/era-interim_stream_fcst1_'//trim(fyr)//'.nc'
  
end subroutine create_ERAinterimLand_fcst1_filename

!BOP
! 
! !ROUTINE: create_ERAinterimLand_anlys_filename
! \label{create_ERAinterimLand_anlys_filename}
!
! !INTERFACE: 
subroutine create_ERAinterimLand_anlys_filename(odir,yr,filename)
! 
! !USES:   
  use LVT_logMod

  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the ERAinterimLand data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            ERAinterimLand base directory
!   \item[yr]              year of data
!   \item[filename]        Name of the ERAinterimLand file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/era-interim_stream_anlys_'//trim(fyr)//'.nc'
  
end subroutine create_ERAinterimLand_anlys_filename


