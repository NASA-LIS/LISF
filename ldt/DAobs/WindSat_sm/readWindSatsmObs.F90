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
! !ROUTINE: readWindSatsmObs
! \label{readWindSatsmObs}
! 
! !REVISION HISTORY: 
!  28 Jan 2013: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readWindSatsmObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc, LDT_domain
  use LDT_timeMgrMod,   only : LDT_get_julss
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use WindSatsm_obsMod, only : WindSatsmobs
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the standard 
! LPRM soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  integer           :: iret
  integer           :: ftn,istat
  integer           :: currTime, refTime
  integer           :: yr,mo,da,hour,min,ss
  real              :: dt
  integer           :: fnd
  character(len=LDT_CONST_PATH_LEN)     :: smname,tsname,tmname,clsname
  logical*1         :: li(WindSatsmobs(n)%mi)
  logical*1         :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: sm(WindSatsmobs(n)%mi)
  real              :: ts(WindSatsmobs(n)%mi)
  real*8            :: tm(WindSatsmobs(n)%mi)
  real              :: tm1(WindSatsmobs(n)%mi)
  integer*2         :: cls(WindSatsmobs(n)%mi)
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LDT_rc%hr)*3600 + 60*LDT_rc%mn + LDT_rc%ss
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  WindSatsmobs(n)%smobs = LDT_rc%udef
  smobs = LDT_rc%udef

  fnd = 0 
  if(WindSatsmobs(n)%startmode.or.alarmCheck) then 
     
     WindSatsmobs(n)%startmode = .false. 

     call create_WindSatsm_filename(WindSatsmobs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, &
          smname,tmname,tsname,clsname)

     inquire(file=trim(smname),exist=file_exists)
     if(file_exists) then

        write(LDT_logunit,*) 'Reading ..',trim(smname)
        ftn = LDT_getNextUnitNumber()
        open(ftn,file=trim(smname),form='unformatted',status='old',&
             access='direct',recl=WindSatsmobs(n)%mi*4,iostat=istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) sm
        endif
        call LDT_releaseUnitNumber(ftn)
        
     endif

     inquire(file=trim(tmname),exist=file_exists)
     if(file_exists) then

        write(LDT_logunit,*) 'Reading ..',trim(tmname)
        ftn = LDT_getNextUnitNumber()
        open(ftn,file=trim(tmname),form='unformatted',status='old',&
             access='direct',recl=WindSatsmobs(n)%mi*8,iostat=istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) tm
        endif
        call LDT_releaseUnitNumber(ftn)
        
     endif

     inquire(file=trim(tsname),exist=file_exists)
     if(file_exists) then

        write(LDT_logunit,*) 'Reading ..',trim(tsname)
        ftn = LDT_getNextUnitNumber()
        open(ftn,file=trim(tsname),form='unformatted',status='old',&
             access='direct',recl=WindSatsmobs(n)%mi*4,iostat=istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) ts
        endif
        call LDT_releaseUnitNumber(ftn)
        
     endif
     inquire(file=trim(clsname), exist=file_exists)
     if(file_exists) then 
        write(LDT_logunit,*) 'Reading ..',trim(clsname)
        ftn = LDT_getNextUnitNumber()
        open(ftn,file=trim(clsname),form='unformatted',status='old',&
             access='direct',recl=WindSatsmobs(n)%mi*8,iostat=istat)
        if(istat.eq.0) then 
           read(ftn,rec=1) cls
        endif
        call LDT_releaseUnitNumber(ftn)
     endif
     
     li = .false.
     do c=1,WindSatsmobs(n)%mi
        if(sm(c).ne.-999.0) then 
           if(cls(c).eq.80.or.cls(c).eq.90.and.ts(c).gt.0) then 
              li(c) = .true.
           endif
        endif
     enddo
     
     call neighbor_interp(LDT_rc%gridDesc(n,:),li,sm,&
          lo,WindSatsmobs(n)%smobs,&
          WindSatsmobs(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          WindSatsmobs(n)%n11,LDT_rc%udef,iret)

     do c=1,WindSatsmobs(n)%mi
        tm1(c) = tm(c)
     enddo
     call neighbor_interp(LDT_rc%gridDesc(n,:),li,tm1,&
          lo,WindSatsmobs(n)%tmobs,&
          WindSatsmobs(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          WindSatsmobs(n)%n11,LDT_rc%udef,iret)

  endif

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(WindSatsmobs(n)%smobs(c,r).ne.LDT_rc%udef) then 
           if(LDT_rc%mo.eq.11.or.LDT_rc%mo.eq.12.or.LDT_rc%mo.eq.1.or.&
                LDT_rc%mo.eq.2) then 
              WindSatsmobs(n)%smobs = LDT_rc%udef
           else
              yr = 2000
              mo = 1
              da = 1
              hour = 12
              min = 0 
              ss = 0 
              
              call LDT_get_julss(yr, mo, da, hour, min, ss, reftime)
              call LDT_get_julss(LDT_rc%yr, LDT_rc%mo, LDT_rc%da, &
                   LDT_rc%hr, 0, 0, currTime)
              dt = WindSatsmobs(n)%tmobs(c,r) + reftime - & 
                   currTime
              if(.not.(WindSatsmobs(n)%tmobs(c,r).gt.0.0.and.&
                   dt .ge.0.and.dt.le.LDT_rc%ts)) then 
                 WindSatsmobs(n)%smobs(c,r) = LDT_rc%udef
              endif
           endif
           
        endif
        if(WindSatsmobs(n)%smobs(c,r).ne.LDT_rc%udef) then 
           fnd = 1
        endif
     enddo
  enddo

!  if(fnd.eq.1) then 
!     open(100,file='smdata.bin',form='unformatted')
!     write(100) WindSatsmobs(n)%smobs
!     close(100)
!     stop
!  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       WindSatsmobs(n)%smobs,vlevel=1)

end subroutine readWindSatsmObs




!BOP
! !ROUTINE: create_WindSatsm_filename
! \label{create_WindSatsm_filename}
! 
! !INTERFACE: 
subroutine create_WindSatsm_filename(ndir,yr, mo,da, smname,&
     tmname,tsname,clsname)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: smname,tmname,tsname,clsname
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the WindSat filenames based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the WindSat soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[smame] name of soil moisture data file
!  \item[tmname] name of the WindSat time file
!  \item[tsname] name of the WindSat surface temperature file
!  \item[clsname] name of the WindSat surface class file
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
 smname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.sm2'

  tmname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.tm2'

  tsname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.ts2'

  clsname = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/d'//trim(fyr)//trim(fmo)//trim(fda)//'GEZ25av_d.cls2'
  
end subroutine create_WindSatsm_filename
