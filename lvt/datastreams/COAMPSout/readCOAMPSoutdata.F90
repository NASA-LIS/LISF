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

!------------------------------------------------------------------------------
! !REVISION HISTORY:
!  09 Dec 2022: Mahdi Navari; Initial Specification in LVT

subroutine readCOAMPSoutdata(source)

! Imports
   use ESMF
   use COAMPSout_dataMod
   use LVT_coreMod
   use LVT_histDataMod
   use LVT_logMod
   use LVT_timeMgrMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

! Defaults
   implicit none
   
! Arguments
   integer,intent(in) :: source

! !DESCRIPTION:
!  Opens and reads 1-hourly COAMPS output forcing. 

!EOP

! Local variables
  character*100           :: fname
  integer                 :: fcsthr, yr, mo, da, hr, mn, ss
  integer                 :: c,r,gindex
  integer                 :: ftn, ios
  logical                 :: file_exists, rainc_exists
  integer                 :: t2Id, q2Id, swdownId, glwId
  integer                 :: u10Id, v10Id, psfcId, rainncId
  integer                 :: raincID
  real                    :: gvar(LVT_rc%lnc,LVT_rc%lnr,1)
  real                    :: t2(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: q2(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: swdown(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: glw(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: u10(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: v10(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: psfc(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: rainnc(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: rainc(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: t2scal,q2scal,u10scal,v10scal,swdscal,glwscal
  real                    :: psfcscal,rainncscal
  real                    :: t2offset,q2offset,u10offset,v10offset
  real                    :: psfcoffset,rainncoffset,swdoffset,glwoffset
  integer                 :: status

gvar = 0.0
t2 = 0.0
q2 = 0.0
swdown = 0.0
glw = 0.0
u10 = 0.0
v10 = 0.0
psfc = 0.0
rainnc = 0.0
rainc = 0.0

! read hourly forcing data
yr = LVT_rc%dyr(source)
mo = LVT_rc%dmo(source)
da = LVT_rc%dda(source)
hr = LVT_rc%dhr(source)
mn = 0
ss = 0

fcsthr = (hr/6)*6
call COAMPSoutfile(fname,COAMPSoutdata(source)%odir,COAMPSoutdata(source)%COAMPSnest_id,fcsthr,&
      yr,mo,da,hr,mn,ss)

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  inquire (file=trim(fname), exist=file_exists)
  if (file_exists) then
     !ferror = 1
     call LVT_verify(nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn), 'nf90_open failed in read_COAMPSout')
     write(LVT_logunit,*) '[INFO] Reading COAMPSout ',trim(fname)
     call LVT_verify(nf90_inq_varid(ftn,'air_temp_2m',t2Id),&
          'nf90_inq_varid failed for air_temp_2m in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,t2Id,"scale_factor",t2scal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,t2Id,"add_offset",t2offset),&
          'nf90_get_att failed for nf90_get_att')

     call LVT_verify(nf90_inq_varid(ftn,'spec_hum_2m',q2Id),&
          'nf90_inq_varid failed for spec_hum_2m in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,q2Id,"scale_factor",q2scal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,q2Id,"add_offset",q2offset),&
          'nf90_get_att failed for nf90_get_att')

     !call LVT_verify(nf90_inq_varid(ftn,'sw_rad_down',swdownId),&
     !     'nf90_inq_varid failed for sw_rad_down in read_COAMPSout')
     ios = nf90_inq_varid(ftn,'sw_rad_down',swdownId)
     if ( ios /= 0) then
         ios = nf90_inq_varid(ftn,'sw_flux_dn',swdownId)
     endif
     call LVT_verify(ios,'nf90_inq_varid failed for sw_rad_down or sw_flux_dn in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,swdownId,"scale_factor",swdscal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,swdownId,"add_offset",swdoffset),&
          'nf90_get_att failed for nf90_get_att')

     !call LVT_verify(nf90_inq_varid(ftn,'lw_rad_down',glwId),&
     !     'nf90_inq_varid failed for lw_rad_down in read_COAMPSout')
     ios = nf90_inq_varid(ftn,'lw_rad_down',glwId)
     if ( ios /= 0) then
         ios = nf90_inq_varid(ftn,'lw_flux_dn',glwId)
     endif
     call LVT_verify(ios,'nf90_inq_varid failed for lw_rad_down or lw_flux_dn in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,glwId,"scale_factor",glwscal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,glwId,"add_offset",glwoffset),&
          'nf90_get_att failed for nf90_get_att')

     call LVT_verify(nf90_inq_varid(ftn,'wind_10m_x',u10Id),&
          'nf90_inq_varid failed for wind_10m_x in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,u10Id,"scale_factor",u10scal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,u10Id,"add_offset",u10offset),&
          'nf90_get_att failed for nf90_get_att')

     call LVT_verify(nf90_inq_varid(ftn,'wind_10m_y',v10Id),&
          'nf90_inq_varid failed for wind_10m_y in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,v10Id,"scale_factor",v10scal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,v10Id,"add_offset",v10offset),&
          'nf90_get_att failed for nf90_get_att')

     call LVT_verify(nf90_inq_varid(ftn,'surf_atm_pres',psfcId),&
          'nf90_inq_varid failed for surf_atm_pres in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,psfcId,"scale_factor",psfcscal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,psfcId,"add_offset",psfcoffset),&
          'nf90_get_att failed for nf90_get_att')

     call LVT_verify(nf90_inq_varid(ftn,'ttl_prcp',rainncId),&
          'nf90_inq_varid failed for ttl_prcp in read_COAMPSout')
     call LVT_verify(nf90_get_att(ftn,rainncId,"scale_factor",rainncscal),&
          'nf90_get_att failed for nf90_get_att')
     call LVT_verify(nf90_get_att(ftn,rainncId,"add_offset",rainncoffset),&
          'nf90_get_att failed for nf90_get_att')

     rainc_exists = .false.
     rainc = 0.0 ! initialize to 0 
     call LVT_verify(nf90_get_var(ftn,t2id,gvar),&
          'nf90_get_var failed for t2 in read_COAMPSout')
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              t2(c,r) = gvar(c,r,1)*t2scal+t2offset
           endif
        enddo
     enddo

     call LVT_verify(nf90_get_var(ftn,q2id,gvar),&
          'nf90_get_var failed for q2 in read_COAMPSout')

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              q2(c,r) = gvar(c,r,1)*q2scal+q2offset
           endif
        enddo
     enddo

     call LVT_verify(nf90_get_var(ftn,swdownid,gvar),&
          'nf90_get_var failed for swdown in read_COAMPSout')

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              swdown(c,r) = gvar(c,r,1)*swdscal+swdoffset
           endif
        enddo
     enddo

     call LVT_verify(nf90_get_var(ftn,glwid,gvar),&
          'nf90_get_var failed for glw in read_COAMPSout')

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              glw(c,r) = gvar(c,r,1)*glwscal+glwoffset
           endif
        enddo
     enddo

     call LVT_verify(nf90_get_var(ftn,u10id,gvar),&
          'nf90_get_var failed for u10 in read_COAMPSout')
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              u10(c,r) = gvar(c,r,1)*u10scal+u10offset
           endif
        enddo
     enddo

     call LVT_verify(nf90_get_var(ftn,v10id,gvar),&
          'nf90_get_var failed for v10 in read_COAMPSout')
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              v10(c,r) = gvar(c,r,1)*v10scal+v10offset
           endif
        enddo
     enddo

     call LVT_verify(nf90_get_var(ftn,psfcid,gvar),&
          'nf90_get_var failed for psfc in read_COAMPSout')
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              psfc(c,r) = gvar(c,r,1)*psfcscal+psfcoffset
           endif
        enddo
     enddo
     call LVT_verify(nf90_get_var(ftn,rainncid,gvar),&
          'nf90_get_var failed for rainnc in read_COAMPSout')
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(gvar(c,r,1).gt.-32768) then
              rainnc(c,r) = gvar(c,r,1)*rainncscal+rainncoffset
           endif
        enddo
     enddo

     call LVT_verify(nf90_close(ftn), &
          'failed to close file in read_COAMPSout')
     write(LVT_logunit,*) '[INFO] Successfully processed ',trim(fname)

  call LVT_logSingleDataStreamVar(LVT_MOC_swdownforc,source,swdown,vlevel=1,units="W/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_lwdownforc,source,glw,vlevel=1,units="W/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,rainnc,vlevel=1,units="kg/m2")

  call LVT_logSingleDataStreamVar(LVT_MOC_tairforc,source,t2,vlevel=1,units="K")

  call LVT_logSingleDataStreamVar(LVT_MOC_qairforc,source,q2,vlevel=1,&
       units="kg/kg")

  call LVT_logSingleDataStreamVar(LVT_MOC_PSURFFORC,source,psfc, &
       vlevel=1,units="Pa")

  call LVT_logSingleDataStreamVar(LVT_MOC_NWINDFORC,source,v10, &
       vlevel=1,units="m/s")

  call LVT_logSingleDataStreamVar(LVT_MOC_EWINDFORC,source,u10, &
       vlevel=1,units="m/s")

  else
     write(LVT_logunit,*) '[ERR] Forcing file '//trim(fname)//' not found'
  endif
#else
  write(LVT_logunit,*) '[ERR] read_COAMPSout requires NetCDF'
  write(LVT_logunit,*) '[ERR] please recompile LIS'
  call LVT_endrun
#endif

   
end subroutine readCOAMPSoutdata

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: COAMPSoutfile
! \label{COAMPSoutfile}
!
! !INTERFACE:
 subroutine COAMPSoutfile(filename,coampsdir,nest,fcsthr,&
      yr,mo,da,hr,mn,ss)

   implicit none
! !ARGUMENTS: 
   character(len=*), intent(out) :: filename
   character(len=*), intent(in)  :: coampsdir
   integer, intent(in)       :: nest
   integer, intent(in)       :: fcsthr
   integer, intent(in)       :: yr,mo,da,hr,mn,ss

! !DESCRIPTION:
!
!EOP

   integer          :: hr1
   character*10     :: ftime1
   character*1      :: fnest
   character*4      :: fyr
   character*2      :: fmo
   character*2      :: fda
   character*2      :: fhr,fhr1
   character*2      :: fmn
   character*2      :: fss

   hr1 = hr - fcsthr
   write(unit=ftime1, fmt='(i4.4,i2.2,i2.2,i2.2)') yr,mo,da,hr
   write(unit=fnest,fmt='(i1.1)') nest
   write(unit=fyr,fmt='(i4.4)') yr
   write(unit=fmo,fmt='(i2.2)') mo
   write(unit=fda,fmt='(i2.2)') da
   write(unit=fhr,fmt='(i2.2)') hr1
   write(unit=fhr1,fmt='(i2.2)') fcsthr
   write(unit=fmn,fmt='(i2.2)') mn
   write(unit=fss,fmt='(i2.2)') ss

   filename = trim(coampsdir)//'/'//&
        trim(fyr)//trim(fmo)//trim(fda)//trim(fhr1)//&
        '/coamps_'//trim(fnest)//'_'//trim(fyr)//trim(fmo)//&
        trim(fda)//trim(fhr1)//'_00'//trim(fhr)//'0000.nc'

 end subroutine COAMPSoutfile

