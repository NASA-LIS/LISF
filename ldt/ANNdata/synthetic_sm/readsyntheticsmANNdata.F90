!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readsyntheticsmANNdata
! \label{readsyntheticsmANNdata}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readsyntheticsmANNdata(n,iomode)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_ANNMod
  use syntheticsm_ANNdataMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: iomode
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the synthetic
! soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  integer           :: ftn
  character(len=LDT_CONST_PATH_LEN) :: fname
  character*3       :: fnest
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  write(fnest,'(i3.3)') n
  alarmCheck = LDT_isAlarmRinging(LDT_rc,"Synthetic sm alarm "//trim(fnest))

  if(alarmCheck) then
     syntheticsmobs(n)%smobs = LDT_rc%udef
     smobs= LDT_rc%udef
        
     syntheticsmobs(n)%startmode = .false. 
     
     call create_ANNsyntheticsm_filename(syntheticsmobs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
     
     inquire(file=trim(fname),exist=file_exists)
     if(file_exists) then
        
        write(LDT_logunit,*) 'Reading ..',trim(fname)
        ftn = LDT_getNextUnitNumber()
        open(ftn,file=trim(fname),form='unformatted')
        read(ftn) smobs
        call LDT_releaseUnitNumber(ftn)
     endif
     
     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              syntheticsmobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
     
     call LDT_logSingleANNdata(n,&
          syntheticsmobs(n)%smobs,  &
          pindex=1, &
          iomode = iomode, &
          name = "SoilMoist",&
          units="m3/m3")
  endif

end subroutine readsyntheticsmANNdata

!BOP
! !ROUTINE: create_ANNsyntheticsm_filename
! \label{create_ANNsyntheticsm_filename}
! 
! !INTERFACE: 
subroutine create_ANNsyntheticsm_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the synthetic filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the synthetic soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated synthetic filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/SOILM_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'0000.bin'
  
end subroutine create_ANNsyntheticsm_filename
