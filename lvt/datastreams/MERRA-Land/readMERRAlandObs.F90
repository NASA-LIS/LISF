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
! !ROUTINE: readMERRAlandObs
! \label{readMERRAlandObs}
!
! !INTERFACE: 
subroutine readMERRAlandObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use MERRAlandobsMod
          
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
  type(ESMF_Time)         :: merralandtime1, merralandtime2, initTime
  type(ESMF_TimeInterval) :: dayInterval
  character(len=100)      :: var_name
  real*8                  :: lis_prevtime
  real                    :: qs(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: qsb(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: qsm(LVT_rc%lnc,LVT_rc%lnr)
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
  real                    :: tskin(LVT_rc%lnc,LVT_rc%lnr)
 
  integer                 :: status

  qs = 0.0
  qsb = 0.0
  qsm = 0.0
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
  tskin = 0.0


  nc = MERRAlandobs(source)%nc
  nr = MERRAlandobs(source)%nr
  
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
         
  if(alarmCheck.or.(LVT_rc%dda(source).ne.&
       MERRAlandobs(source)%da).or.LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(merralandobs(source)%starttime, yy=LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readMerralandobs')
     MERRAlandobs(source)%da = LVT_rc%dda(source)
     
     yr = LVT_rc%dyr(source)
     mo = LVT_rc%dmo(source)
     da = LVT_rc%dda(source)
     hr = LVT_rc%dhr(source)
     mn = LVT_rc%dmn(source)
     ss = LVT_rc%dss(source)
     
     !set back by one day. 
     call LVT_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,-86400)
     
     call process_MERRAlanddata(source, yr, mo, da)
     
  endif
  
  call ESMF_TimeSet(merralandtime1,yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readMERRAlandobs')
  

  t = nint((merralandtime1 - merralandobs(source)%starttime)/&
       merralandobs(source)%ts) + 1
  
  call aggregate_merralandvar(source, t, qs, &
       merralandobs(source)%qs)
  call aggregate_merralandvar(source, t, qsb, &
       merralandobs(source)%qsb)
  call aggregate_merralandvar(source, t, swnet, &
       merralandobs(source)%swnet)
  call aggregate_merralandvar(source, t, qle,&
       merralandobs(source)%qle)
  call aggregate_merralandvar(source, t, qh, &
       merralandobs(source)%qh)
  call aggregate_merralandvar(source, t, frsno,&
       merralandobs(source)%frsno)
  call aggregate_merralandvar(source, t, snod, &
       merralandobs(source)%snod)
  call aggregate_merralandvar(source, t, swe, &
       merralandobs(source)%swe)
  call aggregate_merralandvar(source, t, qg, &
       merralandobs(source)%qg)
  call aggregate_merralandvar(source, t, sfsm, &
       merralandobs(source)%sfsm)
  call aggregate_merralandvar(source, t, rzsm, &
       merralandobs(source)%rzsm)
  call aggregate_merralandvar(source, t, prcp,  &
       merralandobs(source)%prcp)
  call aggregate_merralandvar(source, t, tskin,&
       merralandobs(source)%tskin)
  
  call LVT_logSingleDataStreamVar(LVT_MOC_QS,source,qs,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_QSB,source,qsb,&
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
  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,sfsm,&
       vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_ROOTMOIST,source,rzsm,&
       vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source, prcp,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source,tskin,&
       vlevel=1,units="K")

end subroutine readMERRAlandObs

subroutine process_MERRAlanddata(source, yr, mo, da)

  use LVT_coreMod
  use LVT_logMod
  use MERRAlandobsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  
  integer                :: source
  integer                :: yr
  integer                :: mo
  integer                :: da

  integer                :: ftn
  character*100          :: fname
  logical                :: file_exists
  integer                :: qsid,qsbid, swnetid, qleid, qhid, frsnoid, snodid
  integer                :: sweid, qgid, sfsmid, rzsmid, prcpid, tskinid
  real                   :: qs(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: qsb(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: swnet(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: qle(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: qh(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: frsno(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: snod(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: swe(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: qg(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: sfsm(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: rzsm(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: prcp(merralandobs(source)%nc, merralandobs(source)%nr,24)
  real                   :: tskin(merralandobs(source)%nc, merralandobs(source)%nr,24)
  integer                :: k,iret
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  call create_MERRAland_filename(MERRAlandobs(source)%odir,&
       yr, mo, da, fname)
  
  inquire(file=trim(fname),exist=file_exists) 
  
  if(file_exists) then 
     write(LVT_logunit,*) '[INFO] Reading MERRAland file ',trim(fname)
     
     iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 
        
        call LVT_verify(nf90_inq_varid(ftn,"runoff",qsid),&
             'nf90_inq_varid failed for runoff')
        call LVT_verify(nf90_inq_varid(ftn,"baseflow",qsbid),&
             'nf90_inq_varid failed for baseflow')
        call LVT_verify(nf90_inq_varid(ftn,"swland",swnetid),&
             'nf90_inq_varid failed for swland')
        call LVT_verify(nf90_inq_varid(ftn,"lhland",qleid),&
             'nf90_inq_varid failed for lhland')
!        call LVT_verify(nf90_inq_varid(ftn,"shland",qhid),&
!             'nf90_inq_varid failed for shland')
        call LVT_verify(nf90_inq_varid(ftn,"frsno",frsnoid),&
             'nf90_inq_varid failed for frsno')
        call LVT_verify(nf90_inq_varid(ftn,"snodp",snodid),&
             'nf90_inq_varid failed for snodp')
        call LVT_verify(nf90_inq_varid(ftn,"snomas",sweid),&
             'nf90_inq_varid failed for snomas')
        call LVT_verify(nf90_inq_varid(ftn,"ghland",qgid),&
             'nf90_inq_varid failed for ghland')
        call LVT_verify(nf90_inq_varid(ftn,"sfmc",sfsmid),&
             'nf90_inq_varid failed for sfmc')
!        call LVT_verify(nf90_inq_varid(ftn,"rzmc",rzsmid),&
!             'nf90_inq_varid failed for rzmc')
        call LVT_verify(nf90_inq_varid(ftn,"tsurf",tskinid),&
             'nf90_inq_varid failed for tsurf')
        
        call LVT_verify(nf90_get_var(ftn,qsid,qs),&
             'Error in nf90_get_var runoff')
        call LVT_verify(nf90_get_var(ftn,qsbid,qsb),&
             'Error in nf90_get_var baseflow')
        call LVT_verify(nf90_get_var(ftn,swnetid,swnet),&
             'Error in nf90_get_var swland')
        call LVT_verify(nf90_get_var(ftn,qleid,qle),&
             'Error in nf90_get_var lhland')
!        call LVT_verify(nf90_get_var(ftn,qhid,qh),&
!             'Error in nf90_get_var shland')
        call LVT_verify(nf90_get_var(ftn,frsnoid,frsno),&
             'Error in nf90_get_var frsno')
        call LVT_verify(nf90_get_var(ftn,snodid,snod),&
             'Error in nf90_get_var snodp')
        call LVT_verify(nf90_get_var(ftn,sweid,swe),&
             'Error in nf90_get_var swe')
        call LVT_verify(nf90_get_var(ftn,qgid,qg),&
             'Error in nf90_get_var ghland')
        call LVT_verify(nf90_get_var(ftn,sfsmid,sfsm),&
             'Error in nf90_get_var sfmc')
!        call LVT_verify(nf90_get_var(ftn,rzsmid,rzsm),&
!             'Error in nf90_get_var rzmc')
        call LVT_verify(nf90_get_var(ftn,tskinid,tskin),&
             'Error in nf90_get_var tsurf')
        
        call LVT_verify(nf90_close(ftn))

!        open(100,file='test_inp.bin',form='unformatted')
!        write(100) qle(:,:,1)
!        close(100)
!        stop

        do k=1,24
           call interp_merralandvar2d(source,qs(:,:,k),&
                merralandobs(source)%qs(:,:,k))
           call interp_merralandvar2d(source,qsb(:,:,k),&
                merralandobs(source)%qsb(:,:,k))
           call interp_merralandvar2d(source,swnet(:,:,k),&
                merralandobs(source)%swnet(:,:,k))
           call interp_merralandvar2d(source,qle(:,:,k),&
                merralandobs(source)%qle(:,:,k))
!           call interp_merralandvar2d(source,qh(:,:,k),&
!                merralandobs(source)%qh(:,:,k))
           call interp_merralandvar2d(source,frsno(:,:,k),&
                merralandobs(source)%frsno(:,:,k))
           call interp_merralandvar2d(source,snod(:,:,k),&
                merralandobs(source)%snod(:,:,k))
           call interp_merralandvar2d(source,swe(:,:,k),&
                merralandobs(source)%swe(:,:,k))
           call interp_merralandvar2d(source,qg(:,:,k),&
                merralandobs(source)%qg(:,:,k))
           call interp_merralandvar2d(source,sfsm(:,:,k),&
                merralandobs(source)%sfsm(:,:,k))
!           call interp_merralandvar2d(source,rzsm(:,:,k),&
!                merralandobs(source)%rzsm(:,:,k))
           call interp_merralandvar2d(source,prcp(:,:,k),&
                merralandobs(source)%prcp(:,:,k))
           call interp_merralandvar2d(source,tskin(:,:,k),&
                merralandobs(source)%tskin(:,:,k))
        enddo

     endif
  end if
#endif
end subroutine process_MERRAlanddata

subroutine aggregate_merralandvar(source, t, var, merralandvar)

  use LVT_coreMod

  integer       :: source
  integer       :: t
  real          :: var(LVT_rc%lnc, LVT_rc%lnr)
  real          :: merralandvar(LVT_rc%lnc, LVT_rc%lnr, 24)
  integer       :: c,r

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(var(c,r).ne.LVT_rc%udef) then 
           var(c,r) = merralandvar(c,r,t)
        endif
     enddo
  enddo
  
end subroutine aggregate_merralandvar



!BOP
!
! !ROUTINE: interp_merralandvar2d
! \label{interp_merralandvar2d}
!
! !INTERFACE: 
subroutine interp_merralandvar2d(source, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use MERRAlandobsMod
! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(merralandobs(source)%nc,merralandobs(source)%nr)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the MERRA fields to the 
!  target LVT domain
!
!EOP

  real              :: var_inp_1d(merralandobs(source)%nc*merralandobs(source)%nr)
  logical*1         :: input_bitmap(merralandobs(source)%nc*merralandobs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r
  integer           :: iret

  nc = merralandobs(source)%nc
  nr = merralandobs(source)%nr
  
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
  
  if(LVT_isAtAfinerResolution(merralandobs(source)%datares)) then
     call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d, &
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          merralandobs(source)%rlat, & 
          merralandobs(source)%rlon, &
          merralandobs(source)%n11, &
          LVT_rc%udef, iret)
     
  else
     call upscaleByAveraging(&
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          merralandobs(source)%n11, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d)
     
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
           var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
        endif
     enddo
  enddo

end subroutine interp_merralandvar2d


!BOP
! 
! !ROUTINE: create_MERRAland_filename
! \label{create_MERRAland_filename}
!
! !INTERFACE: 
subroutine create_MERRAland_filename(odir,yr,mo,da,filename)
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
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the MERRAland data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            MERRAland base directory
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[da]              day of data
!   \item[filename]        Name of the MERRAland file
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

  prefix = 'MERRA300'

  filename = trim(odir)//'/'//trim(fyr)//trim(fmo)//'/'&
       //trim(prefix)//'.prod.simul.tavg1_2d_mld_Nx.'//&
       trim(fyr)//trim(fmo)//trim(fda)//'.SUB.nc'
  
end subroutine create_MERRAland_filename


