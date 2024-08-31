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
! !ROUTINE: readAPHROPRCPobs
! \label{readAPHROPRCPobs}
!
! !INTERFACE: 
subroutine readAPHROPRCPobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_verify, LVT_endrun
  use LVT_timeMgrMod, only : LVT_calendar, LVT_tick
  use APHROPRCP_obsMod,    only : aphroprcpobs
  use map_utils

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer,  intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads and processes the gridded APHRO precip
!  data. The data for a given year is read into memory at the 
!  start of a year and indexed into during each day.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 AUG 2017: Sujay Kumar, Initial Specification
!  23 SEP 2019: Yeosang Yoon, Update file name for product for 2007-2015 
!
!EOP

  integer             :: ftn 
  character*100       :: filename
  character*4         :: fyr
  logical             :: file_exists
  real, allocatable   :: rainf(:,:,:)
  real                :: prcp_in(aphroprcpobs(source)%nc*aphroprcpobs(source)%nr)
  logical*1           :: lb(aphroprcpobs(source)%nc*aphroprcpobs(source)%nr)
  logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                :: prcp_out(LVT_rc%lnc*LVT_rc%lnr)
  real                :: prcp_final(LVT_rc%lnc, LVT_rc%lnr)
  integer             :: k,c,r,nt
  integer             :: nid,varid
  integer             :: iret
  type(ESMF_Time)     :: aphrotime1
  integer             :: status
  real                :: timenow
  logical             :: alarmCheck

  prcp_out = LVT_rc%udef
  prcp_final = LVT_rc%udef

  if((aphroprcpobs(source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(aphroprcpobs(source)%startTime,&
          yy=LVT_rc%dyr(source), &
          mm = 1, &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting scan start time')

     aphroprcpobs(source)%yr = LVT_rc%dyr(source)

     if((mod(LVT_rc%dyr(source),4) .eq. 0 .and. &
          mod(LVT_rc%dyr(source), 100).ne.0) &!leap year
          .or.(mod(LVT_rc%dyr(source),400) .eq.0)) then
        nt = 366
     else
        nt = 365
     endif
     
     allocate(rainf(aphroprcpobs(source)%nc,&
          aphroprcpobs(source)%nr,&
          nt))

    if(aphroprcpobs(source)%loc.eq."MA") then 
       write(fyr,fmt='(i4.4)') LVT_rc%dyr(source)
       if(LVT_rc%dyr(source)>=2007) then           ! Yeosang Yoon
          filename = trim(aphroprcpobs(source)%odir)//&
               '/APHRO_MA_025deg_V1101_EXR1.'//trim(fyr)//'.nc' 
       else
          filename = trim(aphroprcpobs(source)%odir)//&
               '/APHRO_MA_025deg_V1101.'//trim(fyr)//'.nc'
       endif
    endif

     inquire(file=trim(filename), exist=file_exists) 
     
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading APHRO data ',&
             trim(filename)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        iret = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
        call LVT_verify(iret, 'Error opening file'//trim(filename))

        iret = nf90_inq_varid(nid, 'precip',varid)
        call LVT_verify(iret, 'Error nf90_inq_varid: precip')

        iret = nf90_get_var(nid,varid, Rainf)
        call LVT_verify(iret, 'Error nf90_get_var: precip')
        
        iret = nf90_close(nid)
        call LVT_verify(iret, 'Error nf90_close')
#endif
        
     endif
     aphroprcpobs(source)%rainf = LVT_rc%udef
     aphroprcpobs(source)%rainf(:,:,1:nt) = rainf(:,:,1:nt)
     deallocate(rainf)

     write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)

  endif

  timenow = float(LVT_rc%dhr(source))*3600 + &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  if(alarmCheck) then 
     call ESMF_TimeSet(aphrotime1, yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
          h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'aphrotime1 set failed')
     
     k = nint((aphrotime1 - aphroprcpobs(source)%starttime)/&
          aphroprcpobs(source)%timestep)+1
     
     lb = .false. 
     prcp_in = LVT_rc%udef
     do r=1,aphroprcpobs(source)%nr
        do c=1,aphroprcpobs(source)%nc
           if(aphroprcpobs(source)%rainf(c,r,k).ge.0) then
              prcp_in(c+(r-1)*aphroprcpobs(source)%nc) = &
                   aphroprcpobs(source)%rainf(c,r,k)
              lb(c+(r-1)*aphroprcpobs(source)%nc) = .true. 
           endif
        enddo
     enddo
     
     call bilinear_interp(LVT_rc%gridDesc,lb,prcp_in, &
          lo,prcp_out, &
          aphroprcpobs(source)%nc*aphroprcpobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr,   &
          aphroprcpobs(source)%rlat, &
          aphroprcpobs(source)%rlon, &
          aphroprcpobs(source)%w11,  &
          aphroprcpobs(source)%w12,  &
          aphroprcpobs(source)%w21,  &
          aphroprcpobs(source)%w22,  &
          aphroprcpobs(source)%n11,  &
          aphroprcpobs(source)%n12,  &
          aphroprcpobs(source)%n21,  &
          aphroprcpobs(source)%n22,  &
          LVT_rc%udef, iret)     
     
     do r=1,LVT_rc%lnr
        do c=1, LVT_rc%lnc
           prcp_final(c,r) = prcp_out(c+(r-1)*LVT_rc%lnc)
        enddo
     enddo
     
     ! Convert mm/day to kg/m2s (note that 1 mm = 1 kg/m2)
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(prcp_final(c,r).ge.0) then 
              prcp_final(c,r) = prcp_final(c,r)/86400.0 !kg/m2s
           else
              prcp_final(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source,prcp_final,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source,prcp_final,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,vlevel=1,units='kg/m2s')
     
  ! Now convert from kg/m2s to kg/m2
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(prcp_final(c,r).ge.0) then 
           prcp_final(c,r) = prcp_final(c,r)*86400.0 !kg/m2
        else
           prcp_final(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  
  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source,prcp_final,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source,prcp_final,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,vlevel=1,units='kg/m2') 
     
end subroutine readAPHROPRCPobs

