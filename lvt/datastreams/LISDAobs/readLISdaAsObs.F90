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
subroutine readLISdaAsObs(source)
! !USES:   
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod
  use LISda_obsMod,    only : lisdaobs

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
  real             :: obsData_inp(lisdaobs(source)%nc,lisdaobs(source)%nr)
  real             :: obsData_inp_1d(lisdaobs(source)%nc*lisdaobs(source)%nr)
  real             :: obsData_out(LVT_rc%lnc,LVT_rc%lnr)
  real             :: obsData_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1        :: li(lisdaobs(source)%nc*lisdaobs(source)%nr)
  logical*1        :: lo(LVT_rc%lnc*LVT_rc%lnr)
  integer          :: t
  integer          :: ftn
  integer          :: c,r
  integer          :: ios

  character*100      :: cdate, cdate1, cda

   
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source), &
        LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)

   write(unit=cdate1, fmt='(i4.4, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source)

   write(unit=cda, fmt='(i2.2)') &
        lisdaobs(source)%obsinstance


   fname = trim(lisdaobs(source)%odir)//'/' & 
        //trim(cdate1)//'/LISDAOBS_'&
        //trim(cdate)//'.a'//trim(cda)//'.d01.1gs4r'
     
  inquire(file=trim(fname),exist=file_exists)
  
  write(LVT_logunit,*) '[INFO] reading DA obs output ',trim(fname)
  if(file_exists) then
     ftn = 11
     open(ftn,file=trim(fname), form='unformatted')
     read(ftn) obsData_inp
     if(lisdaobs(source)%scal.eq.1) then 
        read(ftn) obsData_inp
     endif
     close(ftn)

     obsData_inp_1d = LVT_rc%udef
     li = .false. 
     do r=1,lisdaobs(source)%nr
        do c=1,lisdaobs(source)%nc
           if(obsData_inp(c,r).ne.LVT_rc%udef) then 
              obsData_inp_1d(c+(r-1)*lisdaobs(source)%nc) = &
                   obsData_inp(c,r)
              li(c+(r-1)*lisdaobs(source)%nc) = .true. 
           endif
        enddo
     enddo
     
!     print*, lisdaobs(source)%nc,lisdaobs(source)%nr
!     open(100,file='test_inp.bin',form='unformatted')
!     write(100) obsData_inp
!     close(100)

           
     call neighbor_interp(LVT_rc%gridDesc,li,obsData_inp_1d,&
          lo, obsData_out_1d, lisdaobs(source)%nc*lisdaobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr,&
          lisdaobs(source)%rlat, lisdaobs(source)%rlon, &
          lisdaobs(source)%n11,LVT_rc%udef, ios)

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           obsData_out(c,r) = obsData_out_1d(c+(r-1)*LVT_rc%lnc)
        enddo
     enddo

!     open(100,file='test.bin',form='unformatted')
!     write(100) obsData_out
!     close(100)
!     stop
  else
     write(LVT_logunit,*) '[WARN] Warning: DAobs file ',trim(fname),' does not exist'
     obsData_out = -9999.0
  endif

  if(lisdaobs(source)%obstype.eq.1) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
          obsData_out,vlevel=1,units="m3/m3")
  elseif(lisdaobs(source)%obstype.eq.2) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH, source, &
          obsData_out,vlevel=1,units="m")
  elseif(lisdaobs(source)%obstype.eq.3) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_SWE, source, &
          obsData_out,vlevel=1,units="kg/m2")
  elseif(lisdaobs(source)%obstype.eq.4) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_LAI, source, &
          obsData_out,vlevel=1,units="-")
  endif

end subroutine readLISdaAsObs

