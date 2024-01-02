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
! !ROUTINE: readLISdaAsObs
! \label(readLISdaAsObs)
!
! !INTERFACE:

subroutine readLVTpercentileAsObs(source)
! !USES:   
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod
  use LVTpercentile_obsMod

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

  character*100    :: fname 
  logical          :: file_exists
  real             :: obsData_inp(lvtpercobs(source)%nc,lvtpercobs(source)%nr)
  real             :: obsData_inp_1d(lvtpercobs(source)%nc*lvtpercobs(source)%nr)
  real             :: obsData_out(LVT_rc%lnc,LVT_rc%lnr)
  real             :: dr_category(LVT_rc%lnc,LVT_rc%lnr)
  real             :: obsData_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1        :: li(lvtpercobs(source)%nc*lvtpercobs(source)%nr)
  logical*1        :: lo(LVT_rc%lnc*LVT_rc%lnr)
  integer          :: t
  integer          :: ftn
  integer          :: varid
  integer          :: c,r
  integer          :: ios,iret

  character*100      :: cdate

   
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source), &
        LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)

   fname = trim(lvtpercobs(source)%odir)// & 
        '/Percentile_TS.'&
        //trim(cdate)//'.d01.nc'
     
  inquire(file=trim(fname),exist=file_exists)
  
  dr_category = -9999.0


  if(file_exists) then
  write(LVT_logunit,*) '[INFO] reading LVT output ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then
        call LVT_verify(nf90_inq_varid(ftn,&
             trim(lvtpercobs(source)%var_name),varid),&
             'nf90_inq_varid failed for '//trim(lvtpercobs(source)%var_name))
        call LVT_verify(nf90_get_var(ftn,varid,obsData_inp),&
             'Error in nf90_get_var for '//trim(lvtpercobs(source)%var_name))
     endif
     iret = nf90_close(ftn)
#endif
     
     obsData_inp_1d = LVT_rc%udef
     li = .false. 
     do r=1,lvtpercobs(source)%nr
        do c=1,lvtpercobs(source)%nc
           if(obsData_inp(c,r).ne.LVT_rc%udef) then 
              obsData_inp_1d(c+(r-1)*lvtpercobs(source)%nc) = &
                   obsData_inp(c,r)
              li(c+(r-1)*lvtpercobs(source)%nc) = .true. 
           endif
        enddo
     enddo
     
     call neighbor_interp(LVT_rc%gridDesc,li,obsData_inp_1d,&
          lo, obsData_out_1d, lvtpercobs(source)%nc*lvtpercobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr,&
          lvtpercobs(source)%rlat, lvtpercobs(source)%rlon, &
          lvtpercobs(source)%n11,LVT_rc%udef, ios)
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           obsData_out(c,r) = obsData_out_1d(c+(r-1)*LVT_rc%lnc)
        enddo
     enddo
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(obsData_out(c,r).ne.-9999.0) then 
              if(obsData_out(c,r).le.0.02) then 
                 dr_category(c,r) = 5   !D4
              elseif(obsData_out(c,r).le.0.05) then 
                 dr_category(c,r) = 4   !D3
              elseif(obsData_out(c,r).le.0.10) then 
                 dr_category(c,r) = 3   !D2
              elseif(obsData_out(c,r).le.0.20) then 
                 dr_category(c,r) = 2   !D1
              elseif(obsData_out(c,r).le.0.30) then 
                 dr_category(c,r) = 1   !D0
              endif
           endif
        enddo
     enddo
     
!          open(100,file='test.bin',form='unformatted')
!          write(100) dr_category
!         close(100)
!          stop


  else
     write(LVT_logunit,*) '[WARN] Warning: Percentile file ',trim(fname),' does not exist'
     obsData_out = -9999.0
  endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_PERCENTILE, source, &
       obsData_out,vlevel=1,units="-")

   call LVT_logSingleDataStreamVar(LVT_MOC_DR_CATEGORY,source,dr_category,&
       vlevel=1,units="-")

end subroutine readLVTpercentileAsObs

