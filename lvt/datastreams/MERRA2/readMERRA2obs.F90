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
! !ROUTINE: readMERRA2Obs
! \label{readMERRA2Obs}
!
! !INTERFACE: 
subroutine readMERRA2Obs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use MERRA2obsMod
          
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
!   This plugin processes the Global Land Data Assimilation System (MERRA)
!   version 2 data available from NASA GES-DISC.
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

  real                    :: timenow
  logical                 :: alarmCheck
  integer                 :: c,r, k,nc,nr
  integer                 :: yr, mo, da, hr, mn, ss, doy
  real                    :: gmt
  integer                 :: t
  type(ESMF_Time)         :: merra2time1, merra2time2, initTime
  type(ESMF_TimeInterval) :: dayInterval
  character(len=100)      :: var_name
  real*8                  :: lis_prevtime
  real                    :: qs(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: qsb(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: trunoff(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: swnet(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: qle(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: qh(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: frsno(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: snod(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: swe(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: qg(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: sfsm(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: rzsm(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: prcp(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: tair(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: tskin(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: rnet(LVT_rc%lnc,LVT_rc%lnr)
 
  integer                 :: status

  qs = 0.0
  qsb = 0.0
  swnet = 0.0
  qle = 0.0
  qh = 0.0
  frsno = 0.0
  snod = 0.0
  swe = 0.0
  qg = 0.0
  sfsm = 0.0
  rzsm = 0.0
  prcp = 0.0
  tair = 0.0
  tskin = 0.0
  rnet = LVT_rc%udef
  trunoff = LVT_rc%udef

  nc = MERRA2obs(source)%nc
  nr = MERRA2obs(source)%nr
  
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
         
  if(alarmCheck.or.(LVT_rc%dda(source).ne.&
       MERRA2obs(source)%da).or.LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(merra2obs(source)%starttime, yy=LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readMerra2obs')
     MERRA2obs(source)%da = LVT_rc%dda(source)
     
     yr = LVT_rc%dyr(source)
     mo = LVT_rc%dmo(source)
     da = LVT_rc%dda(source)
     hr = LVT_rc%dhr(source)
     mn = LVT_rc%dmn(source)
     ss = LVT_rc%dss(source)
     
     !set back by one day. 
     call LVT_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,-86400)
     
     call process_MERRA2data(source, yr, mo, da)
     
  endif
  
  call ESMF_TimeSet(merra2time1,yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readMERRA2obs')
  

  t = nint((merra2time1 - merra2obs(source)%starttime)/&
       merra2obs(source)%ts) + 1
  
  call aggregate_merra2var(source, t, qs, &
       merra2obs(source)%qs)
  call aggregate_merra2var(source, t, qsb, &
       merra2obs(source)%qsb)
  call aggregate_merra2var(source, t, swnet, &
       merra2obs(source)%swnet)
  call aggregate_merra2var(source, t, qle,&
       merra2obs(source)%qle)
  call aggregate_merra2var(source, t, qh, &
       merra2obs(source)%qh)
  call aggregate_merra2var(source, t, frsno,&
       merra2obs(source)%frsno)
  call aggregate_merra2var(source, t, snod, &
       merra2obs(source)%snod)
  call aggregate_merra2var(source, t, swe, &
       merra2obs(source)%swe)
  call aggregate_merra2var(source, t, qg, &
       merra2obs(source)%qg)
  call aggregate_merra2var(source, t, sfsm, &
       merra2obs(source)%sfsm)
  call aggregate_merra2var(source, t, rzsm, &
       merra2obs(source)%rzsm)
  call aggregate_merra2var(source, t, prcp,  &
       merra2obs(source)%prcp)
  call aggregate_merra2var(source, t, tair,  &
       merra2obs(source)%tair)
  call aggregate_merra2var(source, t, tskin,&
       merra2obs(source)%tskin)
  
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qle(c,r).ne.LVT_rc%udef) then 
           rnet(c,r) = qle(c,r) + qh(c,r)+qg(c,r)
        endif
     enddo
  enddo
  
  call LVT_logSingleDataStreamVar(LVT_MOC_QS,source,qs,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_QSB,source,qsb,&
       vlevel=1,units="kg/m2s")
  
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qs(c,r).ne.LVT_rc%udef) then 
           trunoff(c,r) = qs(c,r) + qsb(c,r)
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_RUNOFF,source,trunoff,&
       vlevel=1,units="kg/m2s")

  call LVT_logSingleDataStreamVar(LVT_MOC_SWNET,source,swnet,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,qle,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_QH,source,qh,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWCOVER,source,frsno,&
       vlevel=1,units="-")
  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH,source,snod,&
       vlevel=1,units="m")
  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,swe,&
       vlevel=1,units="kg/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_QG,source,qg,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_RNET,source,rnet,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,sfsm,&
       vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_ROOTMOIST,source,rzsm,&
       vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source, prcp,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,source, tair,&
       vlevel=1,units="K")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(prcp(c,r).ne.LVT_rc%udef) then 
           prcp(c,r) = prcp(c,r) * 3600.0
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source, prcp,&
       vlevel=1,units="kg/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,tskin,&
       vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_RADT,source,tskin,&
       vlevel=1,units="K")

end subroutine readMERRA2Obs

subroutine process_MERRA2data(source, yr, mo, da)

  use LVT_coreMod
  use LVT_logMod
  use MERRA2obsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  
  integer                :: source
  integer                :: yr
  integer                :: mo
  integer                :: da

  integer                :: ftn, ftn_flx
  character*100          :: fname, flxname
  logical                :: file_exists
  integer                :: qsid,qsbid, swnetid, qleid, qhid, frsnoid, snodid
  integer                :: sweid, qgid, sfsmid, rzsmid, prcpid, tairid, tskinid
  real                   :: qs(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: qsb(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: swnet(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: qle(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: qh(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: frsno(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: snod(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: swe(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: qg(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: sfsm(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: rzsm(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: prcp(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: tair(merra2obs(source)%nc, merra2obs(source)%nr,24)
  real                   :: tskin(merra2obs(source)%nc, merra2obs(source)%nr,24)
  integer                :: k,iret
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  call create_MERRA2_filename(MERRA2obs(source)%odir,&
       yr, mo, da, fname, flxname)
  
  inquire(file=trim(fname),exist=file_exists) 
  
  if(file_exists) then 
     write(LVT_logunit,*) '[INFO] Reading MERRA2 file ',trim(fname)
     
     iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 
        
        call LVT_verify(nf90_inq_varid(ftn,"RUNOFF",qsid),&
             'nf90_inq_varid failed for RUNOFF')
        call LVT_verify(nf90_inq_varid(ftn,"BASEFLOW",qsbid),&
             'nf90_inq_varid failed for BASEFLOW')
        call LVT_verify(nf90_inq_varid(ftn,"SWLAND",swnetid),&
             'nf90_inq_varid failed for SWLAND')
        call LVT_verify(nf90_inq_varid(ftn,"LHLAND",qleid),&
             'nf90_inq_varid failed for LHLAND')
        call LVT_verify(nf90_inq_varid(ftn,"SHLAND",qhid),&
             'nf90_inq_varid failed for SHLAND')
        call LVT_verify(nf90_inq_varid(ftn,"FRSNO",frsnoid),&
             'nf90_inq_varid failed for FRSNO')
        call LVT_verify(nf90_inq_varid(ftn,"SNODP",snodid),&
             'nf90_inq_varid failed for SNODP')
        call LVT_verify(nf90_inq_varid(ftn,"SNOMAS",sweid),&
             'nf90_inq_varid failed for SNOMAS')
        call LVT_verify(nf90_inq_varid(ftn,"GHLAND",qgid),&
             'nf90_inq_varid failed for GHLAND')
        call LVT_verify(nf90_inq_varid(ftn,"SFMC",sfsmid),&
             'nf90_inq_varid failed for SFMC')
        call LVT_verify(nf90_inq_varid(ftn,"RZMC",rzsmid),&
             'nf90_inq_varid failed for RZMC')

! Yeosang Yoon 
!        call LVT_verify(nf90_inq_varid(ftn,"PRECTOTLAND",prcpid),&
!             'nf90_inq_varid failed for PRECTOTLAND')
        call LVT_verify(nf90_inq_varid(ftn,"TSURF",tskinid),&
             'nf90_inq_varid failed for TSURF')
        
        call LVT_verify(nf90_get_var(ftn,qsid,qs),&
             'Error in nf90_get_var RUNOFF')
        call LVT_verify(nf90_get_var(ftn,qsbid,qsb),&
             'Error in nf90_get_var BASEFLOW')
        call LVT_verify(nf90_get_var(ftn,swnetid,swnet),&
             'Error in nf90_get_var SWLAND')
        call LVT_verify(nf90_get_var(ftn,qleid,qle),&
             'Error in nf90_get_var LHLAND')
        call LVT_verify(nf90_get_var(ftn,qhid,qh),&
             'Error in nf90_get_var SHLAND')
        call LVT_verify(nf90_get_var(ftn,frsnoid,frsno),&
             'Error in nf90_get_var FRSNO')
        call LVT_verify(nf90_get_var(ftn,snodid,snod),&
             'Error in nf90_get_var SNODP')
        call LVT_verify(nf90_get_var(ftn,sweid,swe),&
             'Error in nf90_get_var SWE')
        call LVT_verify(nf90_get_var(ftn,qgid,qg),&
             'Error in nf90_get_var GHLAND')
        call LVT_verify(nf90_get_var(ftn,sfsmid,sfsm),&
             'Error in nf90_get_var SFMC')
        call LVT_verify(nf90_get_var(ftn,rzsmid,rzsm),&
             'Error in nf90_get_var RZMC')
! Yeosang Yoon      
!        call LVT_verify(nf90_get_var(ftn,prcpid,prcp),&
!             'Error in nf90_get_var PRECTOTLAND')
        call LVT_verify(nf90_get_var(ftn,tskinid,tskin),&
             'Error in nf90_get_var TSURF')
        
        call LVT_verify(nf90_close(ftn))

        do k=1,24
           call interp_merra2var2d(source,qs(:,:,k),&
                merra2obs(source)%qs(:,:,k))
           call interp_merra2var2d(source,qsb(:,:,k),&
                merra2obs(source)%qsb(:,:,k))
           call interp_merra2var2d(source,swnet(:,:,k),&
                merra2obs(source)%swnet(:,:,k))
           call interp_merra2var2d(source,qle(:,:,k),&
                merra2obs(source)%qle(:,:,k))
           call interp_merra2var2d(source,qh(:,:,k),&
                merra2obs(source)%qh(:,:,k))
           call interp_merra2var2d(source,frsno(:,:,k),&
                merra2obs(source)%frsno(:,:,k))
           call interp_merra2var2d(source,snod(:,:,k),&
                merra2obs(source)%snod(:,:,k))
           call interp_merra2var2d(source,swe(:,:,k),&
                merra2obs(source)%swe(:,:,k))
           call interp_merra2var2d(source,qg(:,:,k),&
                merra2obs(source)%qg(:,:,k))
           call interp_merra2var2d(source,sfsm(:,:,k),&
                merra2obs(source)%sfsm(:,:,k))
           call interp_merra2var2d(source,rzsm(:,:,k),&
                merra2obs(source)%rzsm(:,:,k))
           call interp_merra2var2d(source,prcp(:,:,k),&
                merra2obs(source)%prcp(:,:,k))
           call interp_merra2var2d(source,tskin(:,:,k),&
                merra2obs(source)%tskin(:,:,k))
        enddo
     endif
  end if

! Yeosang Yoon -----------------------------------------------------
! Read in the flux file fields (*flx):
 inquire(file=trim(flxname),exist=file_exists)

  if(file_exists) then
     write(LVT_logunit,*) '[INFO] Reading MERRA2 file ',trim(flxname)

     iret = nf90_open(path=trim(flxname),mode=NF90_NOWRITE, &
          ncid = ftn_flx)
     if(iret.eq.0) then

        ! Read in the *corrected* total precipitation field:
        if(MERRA2obs(source)%usecorr.eq.1) then
           call LVT_verify(nf90_inq_varid(ftn_flx,'PRECTOTCORR',prcpid), &
                'nf90_inq_varid failed for prectotcorr (flx) in read_merra2')
           call LVT_verify(nf90_get_var(ftn_flx,prcpid, prcp), &
                'nf90_get_var failed for prectotcorr (flx) in read_merra2')

        else  ! Original total precip:
           call LVT_verify(nf90_inq_varid(ftn_flx,'PRECTOT',prcpid), &
                'nf90_inq_varid failed for prectot (flx) in read_merra2')
           call LVT_verify(nf90_get_var(ftn_flx,prcpid, prcp), &
                'nf90_get_var failed for prectot (flx) in read_merra2')
        endif
        
        call LVT_verify(nf90_inq_varid(ftn_flx,'TLML',tairid),&
             'nf90_inq_varid failed for TLML in read_merra2')        
        call LVT_verify(nf90_get_var(ftn_flx,tairid,tair),&
             'nf90_get_var failed for tlml in read_merra2')

        call LVT_verify(nf90_close(ftn_flx))

        do k=1,24
           call interp_merra2var2d(source,tair(:,:,k),&
                merra2obs(source)%tair(:,:,k))

           call interp_merra2var2d_pcp(source,prcp(:,:,k),&
                merra2obs(source)%prcp(:,:,k),&
                MERRA2obs(source)%usescalef)
        enddo
     endif
  end if
 !-----------------------------------------------------------------

#endif
end subroutine process_MERRA2data

subroutine aggregate_merra2var(source, t, var, merra2var)

  use LVT_coreMod

  integer       :: source
  integer       :: t
  real          :: var(LVT_rc%lnc, LVT_rc%lnr)
  real          :: merra2var(LVT_rc%lnc, LVT_rc%lnr, 24)
  integer       :: c,r

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(var(c,r).ne.LVT_rc%udef) then 
           var(c,r) = merra2var(c,r,t)
        endif
     enddo
  enddo
  
end subroutine aggregate_merra2var



!BOP
!
! !ROUTINE: interp_merra2var2d
! \label{interp_merra2var2d}
!
! !INTERFACE: 
subroutine interp_merra2var2d(source, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use MERRA2obsMod
  use LVT_logMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)       
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(merra2obs(source)%nc,merra2obs(source)%nr)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the MERRA fields to the 
!  target LVT domain
!
!EOP
  integer, parameter :: nbins = 100
  integer           :: ftn,pcp1id,pcp2id,pcp3id,pcp4id

  real              :: var_inp_1d(merra2obs(source)%nc*merra2obs(source)%nr)
  logical*1         :: input_bitmap(merra2obs(source)%nc*merra2obs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r,doy
  real              :: pcp_tmp
  integer           :: iret

  nc = merra2obs(source)%nc
  nr = merra2obs(source)%nr
  
  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        if(var_inp(c,r).ne.1e15) then 
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
           input_bitmap(c+(r-1)*nc) = .true. 
        else
           var_inp(c,r) = LVT_rc%udef
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
        endif
     enddo
  enddo
  
    if(LVT_isAtAfinerResolution(merra2obs(source)%datares)) then
       call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
            var_inp_1d, output_bitmap, var_out_1d, &
            nc*nr, &
            LVT_rc%lnc*LVT_rc%lnr, &
            merra2obs(source)%rlat, & 
            merra2obs(source)%rlon, &
            merra2obs(source)%n11, &
            LVT_rc%udef, iret)
       
    else
       call upscaleByAveraging(&
            nc*nr, &
            LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
            merra2obs(source)%n11, input_bitmap, &
            var_inp_1d, output_bitmap, var_out_1d)
       
    endif

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
             var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
          endif
       enddo
    enddo
    
     
end subroutine interp_merra2var2d


!BOP
!
! !ROUTINE: interp_merra2var2d_pcp
! \label{interp_merra2var2d_pcp}
!
! !INTERFACE: 
subroutine interp_merra2var2d_pcp(source, var_inp,var_out,usescalef)
! !USES: 
  use LVT_coreMod
  use MERRA2obsMod
  use LVT_logMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)       
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(merra2obs(source)%nc,merra2obs(source)%nr)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
  integer           :: usescalef
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the MERRA fields to the 
!  target LVT domain
!
!EOP
  integer, parameter :: nbins = 50
  integer           :: ftn,pcp1id,pcp2id,pcp3id,pcp4id

  real              :: var_inp_1d(merra2obs(source)%nc*merra2obs(source)%nr)
  logical*1         :: input_bitmap(merra2obs(source)%nc*merra2obs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r,doy
  integer           :: scale_pcp_choice
  real              :: pcp_tmp
  integer           :: iret

  scale_pcp_choice = usescalef

  nc = merra2obs(source)%nc
  nr = merra2obs(source)%nr
  
  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        if(var_inp(c,r).ne.1e15) then 
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
           input_bitmap(c+(r-1)*nc) = .true. 
        else
           var_inp(c,r) = LVT_rc%udef
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
        endif
     enddo
  enddo

  if(scale_pcp_choice.eq.1) then 

!     open(100,file='test.bin',form='unformatted')
!     write(100) var_inp_1d

     call finddoy(doy, &
          LVT_rc%dyr(source), &
          LVT_rc%dmo(source), &
          LVT_rc%dda(source))

#if(defined USE_NETCDF3 || defined USE_NETCDF4)      

     call LVT_verify(nf90_open(path=merra2obs(source)%scaleffile,&
            mode=nf90_nowrite,ncid=ftn),&
            'failed to open MERRA2 precip scaling factor input file')

       call LVT_verify(nf90_inq_varid(ftn,"REF_XRANGE",pcp1id),&
            'nf90_inq_varid failed for REF_XRANGE')
       call LVT_verify(nf90_get_var(ftn,pcp1id,merra2obs(source)%refxrange,&
            start=(/1,1,LVT_rc%dmo(source),1/),&
            count=(/merra2obs(source)%nc,merra2obs(source)%nr,1,nbins/)),&
            'nf90_get_var failed for REF_XRANGE')

       call LVT_verify(nf90_inq_varid(ftn,"REF_CDF",pcp2id),&
            'nf90_inq_varid failed for REF_CDF')
       call LVT_verify(nf90_get_var(ftn,pcp2id,merra2obs(source)%refcdf,&
            start=(/1,1,LVT_rc%dmo(source),1/),&
            count=(/merra2obs(source)%nc,merra2obs(source)%nr,1,nbins/)),&
            'nf90_get_var failed for PCP_REF_STD')

       call LVT_verify(nf90_inq_varid(ftn,"MERRA2_XRANGE",pcp3id),&
            'nf90_inq_varid failed for PCP_MERRA2_XRANGE')
       call LVT_verify(nf90_get_var(ftn,pcp3id,merra2obs(source)%merraxrange,&
            start=(/1,1,LVT_rc%dmo(source),1/),&
            count=(/merra2obs(source)%nc,merra2obs(source)%nr,1,nbins/)),&
            'nf90_get_var failed for PCP_MERRA2_XRANGE')

       call LVT_verify(nf90_inq_varid(ftn,"MERRA2_CDF",pcp4id),&
            'nf90_inq_varid failed for PCP_MERRA2_CDF')
       call LVT_verify(nf90_get_var(ftn,pcp4id,merra2obs(source)%merracdf,&
            start=(/1,1,LVT_rc%dmo(source),1/),&
            count=(/merra2obs(source)%nc,merra2obs(source)%nr,1,nbins/)),&
            'nf90_get_var failed for PCP_MERRA2_CDF')

       call LVT_verify(nf90_close(ftn),&
            'nf90_close failed in readMERRA2obs')
#endif
       call rescaleWithCDFmatching(&
            merra2obs(source)%nc,&
            merra2obs(source)%nr,&
            nbins,&
            merra2obs(source)%refxrange,&
            merra2obs(source)%merraxrange,&
            merra2obs(source)%refcdf,&
            merra2obs(source)%merracdf,&
            var_inp_1d)
       
!       write(100) var_inp_1d
!       close(100)
!       print*,'here'
!       stop

    endif
    
  
    if(LVT_isAtAfinerResolution(merra2obs(source)%datares)) then
       call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
            var_inp_1d, output_bitmap, var_out_1d, &
            nc*nr, &
            LVT_rc%lnc*LVT_rc%lnr, &
            merra2obs(source)%rlat, & 
            merra2obs(source)%rlon, &
            merra2obs(source)%n11, &
            LVT_rc%udef, iret)
       
    else
       call upscaleByAveraging(&
            nc*nr, &
            LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
            merra2obs(source)%n11, input_bitmap, &
            var_inp_1d, output_bitmap, var_out_1d)
       
    endif

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
             var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
          endif
       enddo
    enddo
  end subroutine interp_merra2var2d_pcp


!BOP
! 
! !ROUTINE: rescaleWithCDFmatching
! \label{rescaleWithCDFmatching}
!
! !INTERFACE:
  subroutine rescaleWithCDFmatching(&
       nc,            & 
       nr,            & 
       nbins,         & 
       ref_xrange,    & 
       merra_xrange,  &       
       ref_cdf,       &
       merra_cdf,     &
       out_value)

    use LVT_logMod
       
    implicit none
! 
! !ARGUMENTS: 
    integer             :: nc
    integer             :: nr
    integer             :: nbins
    real                :: ref_xrange(nc,nr,1,nbins)
    real                :: merra_xrange(nc,nr,1,nbins)
    real                :: ref_cdf(nc,nr,1,nbins)
    real                :: merra_cdf(nc,nr,1,nbins)
    real                :: out_value(nc*nr)
!
! !DESCRIPTION: 
! 
!   This routine rescales the input merraervation data to the ref's
!   climatology so that the cumulative distribution functions (CDFs)
!   of the merraervations and the ref match (for each grid point). 
! 
!   Ref: Reichle and Koster, 2004, Bias reduction in short records of 
!   satellite soil moisture, Geophys. Res. Lett. 31, L19501, 
!   doi:10.1029/2004GL020938. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]               index of the nest
!  \item[nbins]           number of bins used to compute the ref and merra CDFs
!  \item[max\_merra\_value] maximum allowable value of merraervation
!  \item[min\_merra\_value] minimum allowable value of merraervation
!  \item[ref\_xrange]   x-axis values corresponding to the ref CDF
!  \item[merra\_xrange]     x-axis values corresponding to the merra CDF
!  \item[ref\_cdf]      y-axis (CDF) values corresponding to the ref CDF
!  \item[merra\_cdf]        y-axis (CDF) values corresponding to the merra CDF
!  \item[merra\_value]      merraervation value to be rescaled. 
! \end{description}
!EOP

    real                 :: max_merra_value
    real                 :: min_merra_value

    real                 :: ref_delta(nc,nr)
    real                 :: merra_delta(nc,nr)

    integer              :: c,r,i
    integer              :: binval
    integer              :: col,row
    real                 :: cdf_merraval
    real                 :: merra_tmp, merra_in
    integer,dimension(1) :: index_25 , index_75
    real                 :: Lb_xrange, Ub_xrange, iqr_merra, iqr_ref

    ref_delta   = 0 
    merra_delta = 0
   
    do r=1,nr
       do c=1,nc
          if(ref_xrange(c,r,1,2).ne.-9999.0.and.&
               ref_xrange(c,r,1,1).ne.-9999.0.and.&
               merra_xrange(c,r,1,2).ne.-9999.0.and.&
               merra_xrange(c,r,1,1).ne.-9999.0) then 
             ref_delta(c,r) = ref_xrange(c,r,1,2)-ref_xrange(c,r,1,1)
             merra_delta(c,r)   = merra_xrange(c,r,1,2)-merra_xrange(c,r,1,1)
          endif
       enddo
    enddo
    do r=1,nr
       do c=1,nc
          if(ref_delta(c,r).gt.1E-14.and.merra_delta(c,r).gt.1E-14) then 

!             if(c.eq.516.and.r.eq.163) then 
!                do i=1,nbins
!                   print*, merra_xrange(c,r,1,i),merra_cdf(c,r,1,i),&
!                        ref_xrange(c,r,1,i),ref_cdf(c,r,1,i)
!                enddo
!             endif
             min_merra_value = ref_xrange(c,r,1,1)
             max_merra_value = ref_xrange(c,r,1,nbins)
             index_25 = minloc(abs(merra_cdf(c,r,1,:) - 0.25))
             index_75 = minloc(abs(merra_cdf(c,r,1,:) - 0.75))
             Lb_xrange = merra_xrange(c,r,1,index_25(1))
             Ub_xrange = merra_xrange(c,r,1,index_75(1))
             iqr_merra = Ub_xrange - Lb_xrange
             
             index_25 = minloc(abs(ref_cdf(c,r,1,:) - 0.25))
             index_75 = minloc(abs(ref_cdf(c,r,1,:) - 0.75))
             Lb_xrange = ref_xrange(c,r,1,index_25(1))
             Ub_xrange = ref_xrange(c,r,1,index_75(1))
             iqr_ref = Ub_xrange - Lb_xrange
             
! In a normal distribution 50% of data will fall between +/- 0.67448 sigma (1.134896 sigma)
! and 99.7% of data are between the +/-3 sigma (6 sigma) By dividing these two values we
! get 1.134896 / 6 = 0.189. That means IQR is about the 0.19 of the 
! x-range (~ dynamic range of the variable).
! We can say if IQR is less than 0.05 of the x-range then CDF is too 
! steep and it is better to ignore that for CDF matching.  
! NOTE: more tests are needed to determine the best threshold

!             if( iqr_merra .lt. 0.05 * (merra_xrange(c,r,1,nbins)-merra_xrange(c,r,1,1)) .or. &
!                  iqr_ref .lt. 0.05 * (ref_xrange(c,r,1,nbins)-ref_xrange(c,r,1,1)) ) then 
!                merra_delta(c,r) = 0
!             endif
             if(out_value(c+(r-1)*nc).lt.merra_xrange(c,r,1,nbins)) then 
             
                if(out_value(c+(r-1)*nc).ne.-9999.0) then 
                   merra_in = out_value(c+(r-1)*nc)
!                print*,c,r,merra_in,merra_delta(c,r)
                   if(merra_in.gt.1E-14) then 
!                      print*, c,r,out_value(c+(r-1)*nc),merra_xrange(c,r,1,1),&
!                           merra_xrange(c,r,1,nbins),&
!                           merra_delta(c,r)                   
                      binval = nint((out_value(c+(r-1)*nc)-merra_xrange(c,r,1,1))/&
                           merra_delta(c,r))+1
                      
                      if(binval.gt.nbins) binval = nbins
                      if(binval.le.0) binval = 1
                      cdf_merraval = merra_cdf(c,r,1,binval)
                      if(cdf_merraval.gt.1.0) cdf_merraval = 1.0
                      i=1
                      do while((ref_cdf(c,r,1,i).lt.cdf_merraval).and.&
                           (i.le.nbins))
                         i = i+1
                         if(i.gt.nbins) exit
                      enddo
                      if(i.gt.nbins) i = i-1
                      merra_tmp = ref_xrange(c,r,1,i)
!                      if(c.eq.516.and.r.eq.163) then 
!                         print*, merra_in,binval, cdf_merraval, i, merra_tmp
!                      endif
                      if(merra_tmp.gt.max_merra_value) then 
                         !                      print*, 'max const',merra_tmp
                         merra_tmp = -9999.0
                      endif
                      
                      if(merra_tmp.lt.min_merra_value) then 
                         merra_tmp = -9999.0
                      endif
                      out_value(c+(r-1)*nc) = merra_tmp
                   else
                      out_value(c+(r-1)*nc) = merra_in !too low value, do not scale
                   endif
                   if(out_value(c+(r-1)*nc).lt.min_merra_value.and.&
                        out_value(c+(r-1)*nc).ne.-9999.0) then 
                      write(LVT_logunit,*) '[ERR] Problem in CDF scaling '
                      call LVT_endrun()
                   endif
                else
                   out_value(c+(r-1)*nc) = -9999.0
                endif

                if(merra_in.ne.-9999.0.and.out_value(c+(r-1)*nc).eq.-9999.0) then 
                   print*,'mismatch', c,r,merra_in, out_value(c+(r-1)*nc)
                   stop
                endif
             endif
          endif
    enddo
 enddo

end subroutine rescaleWithCDFmatching

!BOP
!
! !ROUTINE: finddoy
! \label{finddoy} 
! 
! !INTERFACE:
  subroutine finddoy(doy,yr,mo,da)

    implicit none
! !ARGUMENTS:
    integer,intent(in)   :: yr,mo,da
    integer,intent(out)  :: doy
!
! !DESCRIPTION:
! 
!  Determines the time, time in GMT, and the day of the year
!  based on the value of year, month, day of month, hour of 
!  the day, minute and second. This method is the inverse of 
!  time2date.
!
!   NOTE: This routine has been known to give round off error
!   problems when attempting to retrieving minutes and seconds
!   from the given number. Use at your own risk!
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[time]
!   lis time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    integer :: yrdays,days(13),k
!    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/
    data days /31,29,31,30,31,30,31,31,30,31,30,31,30/

    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then             !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    
    doy=0
    do k=1,(mo-1)
       doy=doy+days(k)
    enddo
    doy=doy+da
    
!    if(yrdays.eq.366.and.mo.gt.2)doy=doy+1
    
    return
  end subroutine finddoy
!BOP
! 
! !ROUTINE: create_MERRA2_filename
! \label{create_MERRA2_filename}
!
! !INTERFACE: 
subroutine create_MERRA2_filename(odir,yr,mo,da,filename,flxname)
! 
! !USES:   
  use LVT_logMod

  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  character(len=*)             :: filename, flxname
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the MERRA2 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            MERRA2 base directory
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[da]              day of data
!   \item[filename]        Name of the MERRA2 file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  character*50            :: prefix

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  if (yr==1979 .and. mo>=2) then
     prefix = 'MERRA2_100'
  elseif (yr>1979 .and. yr<=1991) then
     prefix = 'MERRA2_100'
  elseif ( yr >= 1992 .and. yr <= 2000 ) then ! Since 2000 is last full year
     prefix = 'MERRA2_200'
!  elseif ( yr >= 2001 .and. yr <= 2009 ) then ! Since 2009 is last full year
  elseif ( yr >= 2001 .and. yr <= 2010 ) then ! Yeosang Yoon (match with merra2 reader in LVT) 
     prefix = 'MERRA2_300'
!  elseif ( yr >= 2010 ) then
  elseif ( yr >= 2011 ) then ! Yeosang Yoon
     prefix = 'MERRA2_400'
  else
!     write(LVT_logunit,*) '[ERR] merra2files: date out of range'
!     write(LVT_logunit,*) '[ERR] Supported years are from 1979-2-1 through ...'
!     call LVT_endrun()	
     filename = "none"
  endif
      
  filename = trim(odir)//'/'//trim(prefix)//'/Y'//trim(fyr)//&
       '/M'//trim(fmo)//'/'//trim(prefix)//'.tavg1_2d_lnd_Nx.'//&
       trim(fyr)//trim(fmo)//trim(fda)//'.nc4'
 
  ! add Yeosang Yoon (read flux file) 
  flxname = trim(odir)//'/'//trim(prefix)//'/Y'//trim(fyr)//&
       '/M'//trim(fmo)//'/'//trim(prefix)//'.tavg1_2d_flx_Nx.'//&
       trim(fyr)//trim(fmo)//trim(fda)//'.nc4'

end subroutine create_MERRA2_filename


