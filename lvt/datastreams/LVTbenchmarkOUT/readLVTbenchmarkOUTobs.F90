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
! !ROUTINE: readLVTbenchmarkOUTobs
! \label(readLVTbenchmarkOUTobs)
!
! !INTERFACE:
subroutine readLVTbenchmarkOUTobs(source)
! !USES:   
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod
  use LVTbenchmarkOUT_obsMod

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
  real             :: obsData(LVT_rc%lnc, LVT_rc%lnr)
  
  integer          :: t
  integer          :: ios
  integer          :: ftn
  integer          :: varid
  integer          :: c,r

  character*100      :: cdate, cdate1

   
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source), &
        LVT_rc%dda(source), LVT_rc%dhr(source),LVT_rc%dmn(source)

   write(unit=cdate1, fmt='(i4.4, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source)
     
   fname = trim(lvtbenchobs(source)%odir)//'/TRAINING/LVT_HIST_OUT_' &         
        //trim(cdate)//'00.nc'
     
  inquire(file=trim(fname),exist=file_exists)
  
  write(LVT_logunit,*) '[INFO] reading LVT output ',trim(fname)
  if(file_exists) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
     ios = nf90_open(path=fname,mode=NF90_NOWRITE, ncid=ftn)
     
     ios = nf90_inq_varid(ftn,trim(lvtbenchobs(source)%vname), varid)
     ios = nf90_get_var(ftn,varid,obsData)
     ios = nf90_close(ftn)
#endif

  else
     write(LVT_logunit,*) '[WARN] Warning: LVT file ',trim(fname),' does not exist'
     obsData = -9999.0
  endif
  if (LVT_MOC_QLE(source).ge.1) then 
     call LVT_logSingleDataStreamVar(LVT_MOC_QLE, source, obsData,vlevel=1,&
          units="W/m2")
  endif

end subroutine readLVTbenchmarkOUTobs

