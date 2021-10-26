!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGEOSland
! \label(readGEOSland)
!
! !INTERFACE:
subroutine readGEOSland(source)
! !USES:   
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use map_utils
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use GEOSlandMod,    only : geoslandoutput

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!

!EOP

  integer, parameter :: ntiles = 149837
  character*200    :: fname 
  logical          :: file_exists
  real             :: lat(ntiles)
  real             :: lon(ntiles)
  real             :: lh(ntiles)
  real             :: lh2d(LVT_rc%lnc,LVT_rc%lnr)
  integer          :: nlh2d(LVT_rc%lnc,LVT_rc%lnr)

  real             :: sfsm(ntiles)
  real             :: sfsm2d(LVT_rc%lnc,LVT_rc%lnr)
  integer          :: nsfsm2d(LVT_rc%lnc,LVT_rc%lnr)

  real             :: rzsm(ntiles)
  real             :: rzsm2d(LVT_rc%lnc,LVT_rc%lnr)
  integer          :: nrzsm2d(LVT_rc%lnc,LVT_rc%lnr)
  
  integer          :: latid
  integer          :: lonid
  integer          :: lhid,sfsmid,rzsmid
  integer          :: t,iret
  integer          :: ftn
  real             :: col,row
  integer          :: c,r

  character*100      :: cdate, cdate1,cdate2

  
  lh2d  = 0
  nlh2d = 0

  sfsm2d = 0
  nsfsm2d = 0

  rzsm2d = 0
  nrzsm2d = 0 
  
  write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
       LVT_rc%dyr(source), LVT_rc%dmo(source), &
       LVT_rc%dda(source)
  
  write(unit=cdate1, fmt='(i4.4)') &
       LVT_rc%dyr(source)
  
  write(unit=cdate2, fmt='(i2.2)') &
       LVT_rc%dmo(source)
  
  
  fname = trim(geoslandoutput(source)%odir)//'/Y' & 
       //trim(cdate1)//'/M'&
       //trim(cdate2)&
       //'/irrig_m09_t0m0_tile1_0.4.tavg24_1d_lnd_Nt.'&
       //trim(cdate)//'_1200z.nc4'
  
  inquire(file=trim(fname),exist=file_exists)
  
  if(file_exists) then
     write(LVT_logunit,*) '[INFO] reading GEOS land output ',trim(fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 
        
        call LVT_verify(nf90_inq_varid(ftn,"lat",latid),&
             'nf90_inq_varid failed for lat')
        call LVT_verify(nf90_inq_varid(ftn,"lon",lonid),&
             'nf90_inq_varid failed for lon')
        call LVT_verify(nf90_inq_varid(ftn,"LHLAND",lhid),&
             'nf90_inq_varid failed for LHLAND')
        call LVT_verify(nf90_inq_varid(ftn,"SFMC",sfsmid),&
             'nf90_inq_varid failed for SFMC')
        call LVT_verify(nf90_inq_varid(ftn,"RZMC",rzsmid),&
             'nf90_inq_varid failed for RZMC')
        
        call LVT_verify(nf90_get_var(ftn,latid, lat),&
             'Error in nf90_get_var lat')
        call LVT_verify(nf90_get_var(ftn,lonid, lon),&
             'Error in nf90_get_var lon')
        call LVT_verify(nf90_get_var(ftn,lhid, lh),&
             'Error in nf90_get_var LHLAND')
        call LVT_verify(nf90_get_var(ftn,sfsmid, sfsm),&
             'Error in nf90_get_var SFMC')
        call LVT_verify(nf90_get_var(ftn,rzsmid, rzsm),&
             'Error in nf90_get_var RZMC')
        
        do t=1,ntiles
           call latlon_to_ij(LVT_domain%lvtproj,&
                lat(t),lon(t),col,row)
           c=nint(col)
           r=nint(row)
           if(c.ge.1.and.c.le.LVT_rc%lnc.and.&
                r.ge.1.and.r.le.LVT_rc%lnr) then 
              lh2d(c,r)  = lh2d(c,r) + lh(t)
              nlh2d(c,r) = nlh2d(c,r)+1

              sfsm2d(c,r) = sfsm2d(c,r) + sfsm(t)
              nsfsm2d(c,r) = nsfsm2d(c,r) + 1

              rzsm2d(c,r) = rzsm2d(c,r) + rzsm(t)
              nrzsm2d(c,r) = nrzsm2d(c,r) + 1
              
           endif
        enddo
        
     endif
     iret = nf90_close(ftn)
#endif
  else
     write(LVT_logunit,*) '[WARN] Warning: GEOS land file ',trim(fname),' does not exist'
     lh = -9999.0
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nlh2d(c,r).gt.0) then
           lh2d(c,r) = lh2d(c,r)/nlh2d(c,r)
           sfsm2d(c,r) = sfsm2d(c,r)/nsfsm2d(c,r)
           rzsm2d(c,r) = rzsm2d(c,r)/nrzsm2d(c,r)           
        else
           lh2d(c,r) = LVT_rc%udef
           sfsm2d(c,r) = LVT_rc%udef
           rzsm2d(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  
  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source,&
       sfsm2d,vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_ROOTMOIST,source,&
       rzsm2d,vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,&
       lh2d,vlevel=1,units="W/m2")
end subroutine readGEOSland

