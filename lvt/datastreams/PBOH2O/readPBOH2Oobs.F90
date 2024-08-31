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
! !ROUTINE: readPBOH2OObs
! \label{readPBOH2OObs}
!
! !INTERFACE: 
subroutine readPBOH2OObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use PBOH2O_obsMod,      only : pboh2oobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
!
  integer,    intent(in) :: source
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for PBOH2O station data. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 May 2014: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! !ARGUMENTS: 

!EOP
  real*8           :: lis_prevtime
  type(ESMF_Time)  :: pboh2otime1, pboh2otime2
  integer          :: k,t
  integer          :: ftn
  character*100    :: filename
  real             :: gmt
  real                :: time
  logical                 :: file_exists
  type(ESMF_TimeInterval) :: dayInterval
  type(ESMF_Time)         :: initTime,obstime1,obstime2,starttime,obstime
  character*100       :: line
  character*4         :: fyr
  integer             :: iloc,i,ios,ios1
  integer             :: obs_yr,obs_mo,obs_da,obs_hr,obs_mn,obs_ss
  real                :: snod(LVT_rc%lnc,LVT_rc%lnr),sm(LVT_rc%lnc,LVT_rc%lnr)
  real                :: swe(LVT_rc%lnc,LVT_rc%lnr)
  real                :: col,row
  real                :: snod_val,sm_val,swe_val
  integer             :: stn_col,stn_row
  integer             :: yr, mo, da, hr, mn, ss, doy
  integer             :: st,et,c,r
  integer             :: status

  snod = LVT_rc%udef
  swe =  LVT_rc%udef
  sm =  LVT_rc%udef


  if(pboh2oobs(source)%yr.ne.LVT_rc%dyr(source).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 

     call ESMF_TimeSet(pboh2oobs(source)%startTime,yy=LVT_rc%dyr(source),&
          mm = 1,dd=1,h = 0, m = 0, calendar=LVT_calendar,&
          rc=status)
     call LVT_verify(status, 'error in ESMF_TimeSet in readPboh2oobs(Source)')     

     pboh2oobs(source)%yr = LVT_rc%dyr(source)
     
     pboh2oobs(source)%sm  = LVT_rc%udef
     pboh2oobs(source)%snod = LVT_rc%udef
     pboh2oobs(source)%swe = LVT_rc%udef

     do k=1,pboh2oobs(source)%n_stns
        write(fyr,fmt='(i4.4)') LVT_rc%dyr(source)
        filename = trim(pboh2oobs(source)%odir)//'/'//&
             trim(fyr)//'/'//&
             trim(pboh2oobs(source)%stn_name(k))//'_'//&
             trim(fyr)//'_v1.csv'
        inquire(file=filename,exist=file_exists)
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading ',trim(filename)
           ftn = LVT_getNextUnitNumber()
           open(ftn,file=trim(filename),form='formatted')
           ios = 0 
           read(ftn,*)
           read(ftn,*)
           read(ftn,*)
           read(ftn,'(a)') line
           iloc = index(line,"#")
           line = line(iloc+1:len(line))
           read(line,*) pboh2oobs(source)%stnlat(k),pboh2oobs(source)%stnlon(k)
           pboh2oobs(source)%stnlon(k) = pboh2oobs(source)%stnlon(k) - 360.0
           read(ftn,*)
           read(ftn,*)
           
           do while(ios.eq.0)
              read(ftn,'(a)',iostat=ios) line
              if(ios.ne.0) exit
              iloc = index(line,"-")
              read(line(1:iloc-1),*) yr
              line = line(iloc+1:len(line))
              iloc = index(line,"-")
              read(line(1:iloc-1),*) mo
              line = line(iloc+1:len(line))
              iloc = index(line,"T")
              read(line(1:iloc-1),*) da
              line = line(iloc+1:len(line))
              iloc = index(line,":")
              read(line(1:iloc-1),*) hr
              line = line(iloc+1:len(line))
              iloc = index(line,":")
              read(line(1:iloc-1),*) mn
              line = line(iloc+1:len(line))
              iloc = index(line,"Z")
              read(line(1:iloc-1),*) ss
              line = line(iloc+1:len(line))
              
              do i=1,13
                 iloc = index(line,",")
                 line = line(iloc+1:len(line))          
              enddo
              
              iloc = index(line,",")
              read(line(1:iloc-1),*,iostat=ios1) snod_val
              if(ios1.ne.0) snod_val = LVT_rc%udef

              line = line(iloc+1:len(line))
              
              do i=1,3
                 iloc = index(line,",")
                 line = line(iloc+1:len(line))          
              enddo
              iloc = index(line,",")
              read(line(1:iloc-1),*,iostat=ios1) swe_val
              if(ios1.ne.0) swe_val = LVT_rc%udef
              line = line(iloc+1:len(line))
              
              iloc = index(line,",")
              line = line(iloc+1:len(line))          
              
              iloc = index(line,",")
              read(line(1:iloc-1),*,iostat=ios1) sm_val
              if(ios1.ne.0) sm_val = LVT_rc%udef
              
!              if(sm_val.ne.LVT_rc%udef) then 
!                 print*, pboh2oobs(source)%stnlat(k),pboh2oobs(source)%stnlon(k),sm_val
!              endif

              call ESMF_TimeSet(obsTime, &
                   yy=yr,mm=mo,dd=da,h=hr,m=mn,s=ss,&
                   calendar=LVT_calendar,rc=status)
              call LVT_verify(status, 'error in ESMF_TimeSet in readPboh2oobs(Source)') 
              
              t = nint((obsTime - pboh2oobs(source)%startTime)/&
                   pboh2oobs(source)%timestep) + 1                   
              
              pboh2oobs(source)%snod(k,t) = snod_val
              pboh2oobs(source)%swe(k,t) = swe_val
              pboh2oobs(source)%sm(k,t) = sm_val
              
           enddo
        endif
        call LVT_releaseUnitNumber(ftn)
     enddo
  end if

  call ESMF_TimeSet(obstime1,yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), dd= LVT_rc%dda(source), h=LVT_rc%dhr(source),&
       m = LVT_rc%dmn(source), s = LVT_rc%dss(source), calendar=LVT_calendar,&
       rc=status)
  call LVT_verify(status, 'ESMF_TimeSet failed in readPboh2oobs(Source)')
  
  t = nint((obstime1 - pboh2oobs(source)%startTime)/pboh2oobs(source)%timestep)+1

  do k=1,pboh2oobs(source)%n_stns
     
     call latlon_to_ij(LVT_domain%lvtproj, &
          pboh2oobs(source)%stnlat(k),pboh2oobs(source)%stnlon(k),&
          col, row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if((stn_col.ge.1.and.stn_col.le.LVT_rc%lnc).and.&
          (stn_row.ge.1.and.stn_row.le.LVT_rc%lnr)) then 
        if(pboh2oobs(source)%snod(k,t).ne.LVT_rc%udef) then 
           snod(stn_col,stn_row) = &
                pboh2oobs(source)%snod(k,t)
        endif
        
        if(pboh2oobs(source)%swe(k,t).ne.LVT_rc%udef) then 
           swe(stn_col,stn_row) =&
                pboh2oobs(source)%swe(k,t)
        endif
        
        if(pboh2oobs(source)%sm(k,t).ne.LVT_rc%udef) then 
           sm(stn_col,stn_row) = &
                pboh2oobs(source)%sm(k,t)
        endif
        
     endif
     
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH, source,snod,vlevel=1,units="m")
  call LVT_logSingleDataStreamVar(LVT_MOC_SWE, source, swe,vlevel=1,units="kg/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST, source, sm,vlevel=1,units="m3/m3")


end subroutine readPBOH2OObs

