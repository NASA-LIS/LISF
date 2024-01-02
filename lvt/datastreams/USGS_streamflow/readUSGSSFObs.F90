!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: readUSGSSFObs
! \label{readUSGSSFObs}
!
! !INTERFACE: 
  subroutine readUSGSSFObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,        only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,     only : LVT_calendar, LVT_tick
  use LVT_logMod,         only : LVT_verify, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_logunit
  use USGSSF_obsMod, only : USGSSFobs
  use map_utils

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
!  11 May 2011: Sujay Kumar, Initial Specification
! 
!EOP
!--------------------------------------------------------------------------



  character*100       :: filename
  character*20        :: datestring, datestring2, qstring
  integer             :: stn_col, stn_row, c,r
  real                :: col,row
  integer             :: tind
  real                :: offset
  real                :: discharge
  integer             :: i,t,kk
  integer             :: ftn
  integer             :: iloc
  integer             :: ios,ios1,ios2,status
  integer             :: yr, mo, da, hr, mn
  character*100       :: line
  logical             :: file_exists
  type(ESMF_Time)     :: obsTime,obsTime1
  real                :: q(LVT_rc%lnc,LVT_rc%lnr)

  q = LVT_rc%udef
!every new year read the data, for each station and store it in memory

  if((LVT_rc%dyr(source).ne.USGSSFobs(source)%yr).or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 
     USGSSFobs(source)%yr = LVT_rc%dyr(source)
     USGSSFobs(source)%q = LVT_rc%udef

     call ESMF_TimeSet(USGSSFobs(source)%startTime, yy=LVT_rc%dyr(source), &
          mm= LVT_rc%dmo(source), &
          dd = 1, &
          h = 0, &
          m = 0, & 
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, "error in setting USGS start time")

     do i = 1, USGSSFobs(source)%n_stns
        call create_USGSSFobs_filename(USGSSFobs(source)%odir, &
             USGSSFobs(source)%stn_name(i), &
             LVT_rc%dyr(source),filename)
        inquire(file=trim(filename), exist=file_exists) 
        if(file_exists) then            
           write(LVT_logunit,*) &
                '[INFO] reading USGS streamflow file ',trim(filename)
           if(USGSSFobs(source)%version.eq.1) then 
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='formatted')
              ios = 0
              do while(ios.eq.0)
                 read(ftn,'(a)',iostat=ios) line
                 if(ios.ne.0) exit
                 if(line(1:4).eq."USGS") then
                    do kk=1,len(line)
                       if(line(kk:kk)==achar(9)) line(kk:kk) =';' !check for tab
                    enddo
                    iloc = index(line,";")
                    line = line(iloc+1:len(line))
                    iloc = index(line,";")
                    line = line(iloc+1:len(line))
                    
                    iloc = index(line,";")
                    read(line(1:iloc-1),*) datestring
                    line = line(iloc+1:len(line))
                    read(datestring,'(i4.4,1X,i2.2,1X,i2.2)') yr, mo, da
                    iloc = index(line,";")
                    do while(iloc.eq.1)
                       iloc = index(line,";")
                       if(iloc.eq.1) then
                          line = line(iloc+1:len(line))
                       endif
                    enddo
                    read(line(1:iloc-1),*,iostat=ios1) qstring
                    discharge = -9999.0
                    if(ios1.eq.0) then
                       read(qstring,*,iostat=ios2) discharge
                       if(ios2.eq.0) then
                          !convert from cubic ft/s to m3/s
                          discharge = discharge*(0.3048)**3
                          
                          call ESMF_TimeSet(obstime,yy=yr,&
                               mm=mo, dd=da, h=0,calendar=LVT_calendar,&
                               rc=status)
                          call LVT_verify(status,'ESMF_TimeSet in readUSGSSFObs')
                          
                          t = nint((obstime-USGSSFobs(source)%starttime)/&
                               USGSSFobs(source)%timestep)+1
                          if(t.ge.1.and.t.le.366) then
                             USGSSFobs(source)%q(i,t) = discharge
                          endif
                       endif
                    else
                    endif
                 endif
              end do
              call LVT_releaseUnitNumber(ftn)
           else
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(filename),form='formatted')
              ios = 0 
              do while(ios.eq.0) 
                 read(ftn,'(a)',iostat=ios) line
                 if(ios.ne.0) exit
                 if(line(1:4).eq."USGS") then 
                    do kk=1,len(line)
                       if(line(kk:kk)==achar(9).or.line(kk:kk)==achar(32))&
                            line(kk:kk) =';' !check for tab and spaces
                    enddo
                    iloc = index(line,";")
                    line = line(iloc+1:len(line))
                    iloc = index(line,";")
                    line = line(iloc+1:len(line))
                    
                    iloc = index(line,";")
                    read(line(1:iloc-1),*) datestring
                    line = line(iloc+1:len(line))              
                    read(datestring,'(i4.4,1X,i2.2,1X,i2.2)') &
                         yr, mo, da
                    
                    iloc = index(line,";")
                    read(line(1:iloc-1),*) datestring2
                    line = line(iloc+1:len(line))              
                    read(datestring2,'(i2.2,1X,i2.2)') &
                         hr,mn
                    
                    iloc = index(line,";")
                    line = line(iloc+1:len(line))
                    iloc = index(line,";")
                    
                    do while(iloc.eq.1) 
                       iloc = index(line,";")
                       if(iloc.eq.1) then 
                          line = line(iloc+1:len(line))     
                       endif
                    enddo
                    read(line(1:iloc-1),*,iostat=ios1) qstring
                    discharge = -9999.0
                    if(ios1.eq.0) then 
                       read(qstring,*,iostat=ios2) discharge        
                       if(ios2.eq.0) then 
                          !convert from cubic ft/s to m3/s      
                          discharge = discharge*(0.3048)**3
                          
                          call ESMF_TimeSet(obstime,yy=yr,&
                               mm=mo, dd=da, h=hr,m=mn,&
                               calendar=LVT_calendar,&
                               rc=status)
                          call LVT_verify(status,'ESMF_TimeSet in readUSGSSFObs')
                          
                          t = nint((obstime-USGSSFobs(source)%starttime)/&
                               USGSSFobs(source)%timestep)+1
                          if(t.ge.1.and.t.le.35136) then 
                             USGSSFobs(source)%q(i,t) = discharge
                          endif
                       endif
                    else
                    endif
                 endif
              end do
              call LVT_releaseUnitNumber(ftn)
           endif
        end if
     enddo
  endif

  call ESMF_TimeSet(obstime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'obstime1 set failed')

  offset = (obstime1-USGSSFobs(source)%starttime)/USGSSFobs(source)%timestep

  if((nint(offset)-offset).eq.0) then !only when LIS time matches the observation
     
     tind = nint((obstime1-USGSSFobs(source)%starttime)/&
          USGSSFobs(source)%timestep)+1

     do i=1,USGSSFobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, &
             USGSSFobs(source)%stnlat(i),USGSSFobs(source)%stnlon(i),&
             col,row)
        stn_col = nint(col)
        stn_row = nint(row)

        if(stn_col.gt.0.and.stn_row.gt.0.and.tind.ge.0.and.&
             stn_col.le.LVT_rc%lnc.and.stn_row.le.LVT_rc%lnr) then 
           if(USGSSFobs(source)%q(i,tind).ne.LVT_rc%udef) then 
              q(stn_col, stn_row) = &
                   USGSSFobs(source)%q(i,tind)
           endif
        endif
     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_streamflow,source, q,vlevel=1,units="m3/s")

end subroutine readUSGSSFObs

!BOP
! 
! !ROUTINE: create_USGSSFobs_filename
! \label(create_USGSSFobs_filename)
!
! !INTERFACE:
subroutine create_USGSSFobs_filename(odir,stn_name,yr,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*),     intent(in) :: odir
  character(len=*),     intent(in) :: stn_name
  integer,              intent(in) :: yr
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

  
  character(len=*)                 :: filename

  character*4      :: fyr

  write(unit=fyr,fmt='(i4.4)') yr
  
  filename = trim(odir)//'/'//trim(fyr)//'/'//trim(stn_name)//&
       '.'//trim(fyr)//'.dat'

end subroutine create_USGSSFobs_filename
