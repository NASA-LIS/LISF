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
! !ROUTINE: readLIS6output
! \label(readLIS6output)
!
! !INTERFACE:
subroutine readLIS6output(source)
! !USES:   
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod
  use LIS6outputMod,    only : lis6output

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
! !NOTES: Currently limited to the LIS6 outputs produced for the 
! NOHRSC with Noah 3.2
! 
!EOP

  character*100    :: fname 
  logical          :: file_exists
  real             :: sca(LVT_rc%lnc, LVT_rc%lnr)
  real             :: snd(LVT_rc%lnc, LVT_rc%lnr)
  integer          :: scaid, sndid
  integer          :: t,iret
  integer          :: ftn
  integer          :: c,r

  character*100      :: cdate, cdate1

   
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') &
        LVT_rc%dyr(source), LVT_rc%dmo(source), &
        LVT_rc%dda(source)

   write(unit=cdate1, fmt='(i4.4)') &
        LVT_rc%dyr(source)
     
   fname = trim(lis6output(source)%odir)//'/' & 
        //trim(cdate1)//'/NOAH32.'&
        //trim(cdate)//'0000.d01.nc'
     
  inquire(file=trim(fname),exist=file_exists)
  
  if(file_exists) then
     write(LVT_logunit,*) '[INFO] reading LIS6 output ',trim(fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 
        
        call LVT_verify(nf90_inq_varid(ftn,"Snowcover",scaid),&
             'nf90_inq_varid failed for Snowcover')
        call LVT_verify(nf90_inq_varid(ftn,"SnowDepth",sndid),&
             'nf90_inq_varid failed for SnowDepth')
        
        call LVT_verify(nf90_get_var(ftn,scaid, sca),&
             'Error in nf90_get_var Snowcover')
        call LVT_verify(nf90_get_var(ftn,sndid,snd),&
             'Error in nf90_get_var SnowDepth')
     endif
     iret = nf90_close(ftn)
#endif
  else
     write(LVT_logunit,*) '[WARN] Warning: LIS6 file ',trim(fname),' does not exist'
     snd = -9999.0
     sca = -9999.0
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_snowcover, source,sca,vlevel=1,units="-")
  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH, source,snd,vlevel=1,units="m")
  
end subroutine readLIS6output

