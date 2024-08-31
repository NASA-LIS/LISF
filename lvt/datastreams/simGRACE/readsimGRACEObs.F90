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
! !ROUTINE: readsimGRACEObs
! \label{readsimGRACEObs}
!
! !INTERFACE: 
subroutine readsimGRACEObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use simGRACE_obsMod, only : simGRACEObs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The simGRACE output is available at monthly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP

  character*100           :: filename
  logical                 :: file_exists
  integer                 :: c,r,t,kk
  integer                 :: ftn
  integer                 :: iret
  real                    :: tws(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: obsData_inp(simGRACEobs(source)%nc,&
       simGRACEobs(source)%nr)
  real             :: obsData_inp_1d(simGRACEobs(source)%nc*simGRACEobs(source)%nr)
  real             :: obsData_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1        :: li(simGRACEobs(source)%nc*simGRACEobs(source)%nr)
  logical*1        :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real,  allocatable      :: tws_in(:)
  logical*1, allocatable  :: bitmap_in(:)
  logical*1               :: bitmap_out(LVT_rc%lnc*LVT_rc%lnr)
  type(ESMF_Time)         :: cTime
  type(ESMF_TimeInterval) :: tw
  integer                 :: yr,mo,da,hr,mn,ss
  integer                 :: status

  tws = LVT_rc%udef

  if(simGRACEobs(source)%useRawData.eq.1) then 
     
     if(simGRACEobs(source)%mo.ne.LVT_rc%d_nmo(source).and.&
          LVT_rc%dda(source).eq.1) then 
        !             LVT_rc%dda(source).eq.15) then 
        
        if(simGRACEobs(source)%startFlag) then 
           simGRACEobs(source)%yr = LVT_rc%dyr(source)
           simGRACEobs(source)%mo = LVT_rc%dmo(source)
           simGRACEobs(source)%startFlag = .false. 
        endif
        
        call create_simGRACE_Raw_filename(simGRACEobs(source)%odir, &
             simGRACEobs(source)%yr,&
             simGRACEobs(source)%mo, filename)
        
        inquire(file=trim(filename),exist=file_exists) 
        
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading simGRACE file ',trim(filename)
           
           allocate(tws_in(simGRACEObs(source)%nc*simGRACEobs(source)%nr))
           allocate(bitmap_in(simGRACEObs(source)%nc*simGRACEobs(source)%nr))
           
           ftn = LVT_getNextUnitNumber()
           open(ftn,file=trim(filename),form='unformatted')
           read(ftn) tws_in
           call LVT_releaseUnitNumber(ftn)
           
           bitmap_in = .true. 
           do t=1,simGRACEObs(source)%nc*simGRACEobs(source)%nr
              if(tws_in(t).eq.-9999.0) then 
                 bitmap_in(t) = .false.
              endif
           enddo
           
           call bilinear_interp(LVT_rc%gridDesc,bitmap_in,tws_in,    &
                bitmap_out,tws,simGRACEobs(source)%nc*simGRACEobs(source)%nr,&
                LVT_rc%lnc*LVT_rc%lnr,   &
                simGRACEobs(source)%rlat,simGRACEobs(source)%rlon, &
                simGRACEobs(source)%w11,simGRACEobs(source)%w12,   &
                simGRACEobs(source)%w21,simGRACEobs(source)%w22,   &
                simGRACEobs(source)%n11,simGRACEobs(source)%n12,   &
                simGRACEobs(source)%n21,simGRACEobs(source)%n22,   &
                LVT_rc%udef,iret)
           
        else
           tws = LVT_rc%udef
        endif
        
        simGRACEobs(source)%yr = LVT_rc%d_nyr(source)
        simGRACEobs(source)%mo = LVT_rc%d_nmo(source)
        
     else
        tws  = LVT_rc%udef
     endif
  else
     if(simGRACEobs(source)%config.eq."default".or.&
          simGRACEobs(source)%config.eq."follow-on") then 
        if(simGRACEobs(source)%mo.ne.LVT_rc%d_nmo(source).and.&
             LVT_rc%dda(source).eq.1) then 
!             LVT_rc%dda(source).eq.15) then 
           
           if(simGRACEobs(source)%startFlag) then 
              simGRACEobs(source)%yr = LVT_rc%dyr(source)
              simGRACEobs(source)%mo = LVT_rc%dmo(source)
              simGRACEobs(source)%startFlag = .false. 
           endif
           
           call create_simGRACE1_filename(simGRACEobs(source)%odir, &
                simGRACEobs(source)%yr,&
                simGRACEobs(source)%mo, filename)
           
           inquire(file=trim(filename),exist=file_exists) 
           
           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading simGRACE file ',trim(filename)
              
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='unformatted')
              read(ftn) obsData_inp
              call LVT_releaseUnitNumber(ftn)
              
           else
              obsData_inp = LVT_rc%udef
           endif
           
           simGRACEobs(source)%yr = LVT_rc%d_nyr(source)
           simGRACEobs(source)%mo = LVT_rc%d_nmo(source)
           
           obsData_inp_1d = LVT_rc%udef

           li = .false. 
           do r=1,simGRACEobs(source)%nr
              do c=1,simGRACEobs(source)%nc
                 if(obsData_inp(c,r).ne.LVT_rc%udef) then 
                    obsData_inp_1d(c+(r-1)*simGRACEobs(source)%nc) = &
                         obsData_inp(c,r)
                    li(c+(r-1)*simGRACEobs(source)%nc) = .true. 
                 endif
              enddo
           enddo

           call neighbor_interp(LVT_rc%gridDesc,li,obsData_inp_1d,&
                lo, obsData_out_1d, simGRACEobs(source)%nc*simGRACEobs(source)%nr, &
                LVT_rc%lnc*LVT_rc%lnr,&
                simGRACEobs(source)%rlat1, simGRACEobs(source)%rlon1, &
                simGRACEobs(source)%n111,LVT_rc%udef, status)
           
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 tws(c,r) = obsData_out_1d(c+(r-1)*LVT_rc%lnc)
              enddo
           enddo
           
        else
           tws  = LVT_rc%udef
        endif
     elseif(simGRACEobs(source)%config.eq."quick-look") then 
        call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h  = LVT_rc%dhr(source), &
             m  = LVT_rc%dmn(source), & 
             s  = LVT_rc%dss(source), &
             calendar = LVT_calendar, & 
             rc = status)
        call LVT_verify(status)
        !1 day
        call ESMF_TimeIntervalSet(tw,s=86400,rc=status)
        ctime = ctime - tw
        
        call ESMF_TimeGet(ctime,yy=yr,mm=mo,dd=da,calendar=LVT_calendar,&
             rc=status)
        call LVT_verify(status)
        
        if(simGRACEobs(source)%da.ne.da) then 

           if(simGRACEobs(source)%startFlag) then 
              simGRACEobs(source)%yr = yr
              simGRACEobs(source)%mo = mo
              simGRACEobs(source)%da = da
              simGRACEobs(source)%startFlag = .false. 
           endif
           
           call create_simGRACE3_filename(simGRACEobs(source)%odir, yr, &
                mo, da, filename)
           
           inquire(file=trim(filename),exist=file_exists) 
           
           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading simGRACE file ',trim(filename)
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='unformatted')
              read(ftn) obsData_inp
              call LVT_releaseUnitNumber(ftn)
              
           else
              obsData_inp = LVT_rc%udef
           endif
           
           simGRACEobs(source)%yr = LVT_rc%d_nyr(source)
           simGRACEobs(source)%mo = LVT_rc%d_nmo(source)
           
           obsData_inp_1d = LVT_rc%udef

           li = .false. 
           do r=1,simGRACEobs(source)%nr
              do c=1,simGRACEobs(source)%nc
                 if(obsData_inp(c,r).ne.LVT_rc%udef) then 
                    obsData_inp_1d(c+(r-1)*simGRACEobs(source)%nc) = &
                         obsData_inp(c,r)
                    li(c+(r-1)*simGRACEobs(source)%nc) = .true. 
                 endif
              enddo
           enddo

           call neighbor_interp(LVT_rc%gridDesc,li,obsData_inp_1d,&
                lo, obsData_out_1d, &
                simGRACEobs(source)%nc*simGRACEobs(source)%nr, &
                LVT_rc%lnc*LVT_rc%lnr,&
                simGRACEobs(source)%rlat1, simGRACEobs(source)%rlon1, &
                simGRACEobs(source)%n111,LVT_rc%udef, status)
           
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 tws(c,r) = obsData_out_1d(c+(r-1)*LVT_rc%lnc)
              enddo
           enddo
           
        else
           tws  = LVT_rc%udef
        endif              
        simGRACEobs(source)%yr = yr
        simGRACEobs(source)%mo = mo
        simGRACEobs(source)%da = da           
        
     else
        call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h  = LVT_rc%dhr(source), &
             m  = LVT_rc%dmn(source), & 
             s  = LVT_rc%dss(source), &
             calendar = LVT_calendar, & 
             rc = status)
        call LVT_verify(status)
        !7 day
        call ESMF_TimeIntervalSet(tw,s=604800,rc=status)
        ctime = ctime - tw
        
        call ESMF_TimeGet(ctime,yy=yr,mm=mo,dd=da,calendar=LVT_calendar,&
             rc=status)
        call LVT_verify(status)
        
        if(simGRACEobs(source)%da.ne.da) then 
           if(simGRACEobs(source)%startFlag) then 
              simGRACEobs(source)%yr = yr
              simGRACEobs(source)%mo = mo
              simGRACEobs(source)%da = da
              simGRACEobs(source)%startFlag = .false. 
           endif
           
           call create_simGRACE2_filename(simGRACEobs(source)%odir, yr, &
                mo, da, filename)
           
           inquire(file=trim(filename),exist=file_exists) 
           
           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading simGRACE file ',trim(filename)
              
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='unformatted')
              read(ftn) tws
              call LVT_releaseUnitNumber(ftn)
              
           else
              tws = LVT_rc%udef
           endif
           simGRACEobs(source)%yr = yr
           simGRACEobs(source)%mo = mo
           simGRACEobs(source)%da = da
           
        endif
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_TWS,source,tws,vlevel=1,units="mm")
  
end subroutine readsimGRACEObs

!BOP
! 
! !ROUTINE: create_simGRACE1_filename
! \label{create_simGRACE1_filename}
!
! !INTERFACE: 
subroutine create_simGRACE1_filename(odir,yr,mo,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for simGRACE_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      simGRACE base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the simGRACE_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr,mo
  character(len=*)             :: filename
!EOP

  character*6             :: fdate
  
  write(unit=fdate, fmt='(i4.4,i2.2)') yr,mo
  
!  filename = trim(odir)//'/simGRACE_obs_'//trim(fdate)//'.bin'
  filename = trim(odir)//'/GRACE_obs_'//trim(fdate)//'.bin'
  
end subroutine create_simGRACE1_filename


!BOP
! 
! !ROUTINE: create_simGRACE2_filename
! \label{create_simGRACE2_filename}
!
! !INTERFACE: 
subroutine create_simGRACE2_filename(odir,yr,mo,da, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for simGRACE_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      simGRACE base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the simGRACE_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr,mo,da
  character(len=*)             :: filename
!EOP

  character*8           :: fdate
  
  write(unit=fdate, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
  
  filename = trim(odir)//'/simGRACE_obs_'//trim(fdate)//'.bin'
  
end subroutine create_simGRACE2_filename

!BOP
! 
! !ROUTINE: create_simGRACE3_filename
! \label{create_simGRACE3_filename}
!
! !INTERFACE: 
subroutine create_simGRACE3_filename(odir,yr,mo,da, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for simGRACE_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      simGRACE base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the simGRACE_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr,mo,da
  character(len=*)             :: filename
!EOP

  character*8           :: fdate
  
  write(unit=fdate, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
  
  filename = trim(odir)//'/GRACE_obs_'//trim(fdate)//'.bin'
  
end subroutine create_simGRACE3_filename


!BOP
! 
! !ROUTINE: create_simGRACE_raw_filename
! \label{create_simGRACE_raw_filename}
!
! !INTERFACE: 
subroutine create_simGRACE_raw_filename(odir,yr,mo,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for simGRACE_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      simGRACE base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the simGRACE_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr,mo
  character(len=*)             :: filename
!EOP

  character*6             :: fdate
  
  write(unit=fdate, fmt='(i4.4,i2.2)') yr,mo
  
  filename = trim(odir)//'/simsimGRACE_'//trim(fdate)//'0112.bin'
  
end subroutine create_simGRACE_raw_filename
