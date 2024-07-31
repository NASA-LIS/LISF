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
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
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
!BOP
! 
! !ROUTINE: readSNODEPmetObs
! \label{readSNODEPmetobs}
! 
! !REVISION HISTORY: 
!  13 Dec 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readSNODEPmetObs(source)
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use SNODEP_metobsMod,    only : SNODEPmetobs

  implicit none
! !ARGUMENTS: 
  integer,   intent(in)   :: source
! 
! !DESCRIPTION: 
!   This subroutine provides the data reader for in-situ snow depth observations
!   of SNODEP. The measurements are assumed to be available at an hourly 
!   interval. 
!
!  !TODO: 
!   filter out bad values -- 996, etc? 
!EOP

  real                    :: time
  integer                 :: status
  type(ESMF_Time)         :: snodeptime1, snodeptime2
  type(ESMF_TimeInterval) :: dayInterval
  type(ESMF_Time)         :: initTime
  real                    :: snod(LVT_rc%lnc, LVT_rc%lnr)
  real                    :: gmt
  integer                 :: c,r,st,et,t
  integer                 :: yr, mo, da, hr, mn, ss, doy
  real*8                  :: lis_prevtime

  snod  = LVT_rc%udef 

  time = LVT_rc%dhr(source)*3600+LVT_rc%dmn(source)*60+LVT_rc%dss(source)
  if((mod(time,86400.0).eq.0.0).or.&
       (LVT_rc%dda(source).ne.SNODEPmetobs(source)%da).or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 

     call ESMF_TimeSet(SNODEPmetobs(source)%starttime, yy=LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readSNODEPmetobs')

     SNODEPmetobs(source)%da = LVT_rc%dda(source)        
     call read_snodepdata(source, SNODEPmetobs(source)%odir,&
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source))

  endif

  call ESMF_TimeSet(snodeptime1, yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readSNODEPmetobs')
  
  t = nint((snodeptime1 - SNODEPmetobs(source)%starttime)/snodepmetobs(source)%ts) + 1
  do r=1, LVT_rc%lnr
     do c=1, LVT_rc%lnc
        if(SNODEPmetobs(source)%snod(c,r,t).ne.LVT_rc%udef) then 
           snod(c,r) = SNODEPmetobs(source)%snod(c,r,t)
        endif
     enddo
  enddo
  do r=1, LVT_rc%lnr
     do c=1, LVT_rc%lnc
!discard all data between june and oct and some additional 
!(conservative) QCs.
        if(LVT_rc%dmo(source).ge.6.and.LVT_rc%dmo(source).le.9) then 
           snod(c,r) = LVT_rc%udef
        elseif(LVT_rc%dmo(source).eq.5.and.snod(c,r).gt.0.5) then
           snod(c,r) = LVT_rc%udef
        elseif(LVT_rc%dmo(source).eq.4.and.snod(c,r).gt.1.0) then 
           snod(c,r) = LVT_rc%udef
        elseif(snod(c,r).gt.2.0) then 
           snod(c,r) = LVT_rc%udef
        endif        
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_snowdepth, source, snod,vlevel=1,units="m")

end subroutine readSNODEPmetObs

!BOP
! 
! !ROUTINE: read_snodepdata
!  \label{read_snodepdata} 
!
! !INTERFACE:
subroutine read_snodepdata(source, odir, yr, mo, da)
! 
! !USES:  
  use ESMF
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_timeMgrMod, only : LVT_calendar
  use LVT_logMod,     only : LVT_logunit, LVT_verify, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber
  use SNODEP_metobsMod,  only : SNODEPmetobs
  use map_utils 

  implicit none
!
! !INPUT PARAMETERS: 
  integer,    intent(in) :: source
  integer,    intent(in) :: yr
  integer,    intent(in) :: mo
  integer,    intent(in) :: da

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine reads an individual snodep file (data for a 
!   particular day). This routine also computes the temporal 
!   offset of the data relative 0z of a particular day. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      base directory for SNODEP data
!   \item[yr]        year of SNODEP data 
!   \item[mo]        month of SNODEP data
!   \item[da]        day of SNODEP data
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)       :: odir
!EOP
  logical                 :: file_exists
  character*100           :: line
  integer                 :: iloc, t,ios, ios1, ios2
  integer                 :: yy,mm,dd,hr
  real                    :: lat,lon,snod
  character*4             :: fyr
  character*2             :: fmo, fda
  integer                 :: ftn 
  integer                 :: stn_col, stn_row
  real                    :: col,row
  integer                 :: status
  integer                 :: data_index
  character*100           :: filename
  real                    :: tempc, dewpc, avgc
  type(ESMF_Time)         :: datatime, starttime

  write(fyr,fmt='(i4.4)') yr
  write(fmo,fmt='(i2.2)') mo
  write(fda,fmt='(i2.2)') da
  
  ftn = LVT_getNextUnitNumber()
  filename = trim(odir)//'/'//trim(fyr)//'/snodep_'//&
       trim(fyr)//trim(fmo)//trim(fda)//'.dat'
  inquire(file=trim(filename), exist=file_exists) 

  if(file_exists) then
 
     SNODEPmetobs(source)%snod = LVT_rc%udef

     open(ftn,file=trim(filename))
     write(LVT_logunit, *) '[INFO] Reading SNODEP ',trim(filename)
     read(ftn,'(a)') line
     ios = 0 
     do while(ios.eq.0) 
        read(ftn,'(a)',iostat=ios) line

        if(ios.ne.0) exit
        
        iloc = index(line,",")
        line = line(iloc+1:len(line))
        !year
        iloc = index(line,",")
        read(line(1:iloc-1),*) yy
        line = line(iloc+1:len(line))
        !month
        iloc = index(line,",")
        read(line(1:iloc-1),*) mm
        line = line(iloc+1:len(line))
        !day
        iloc = index(line,",")
        read(line(1:iloc-1),*) dd
        line = line(iloc+1:len(line))
        !hour
        iloc = index(line,",")
        read(line(1:iloc-1),*) hr
        line = line(iloc+1:len(line))
        
        !lat
        iloc = index(line,",")
        read(line(1:iloc-1),*) lat
        line = line(iloc+1:len(line))
        
        !lon
        iloc = index(line,",")
        read(line(1:iloc-1),*) lon
        line = line(iloc+1:len(line))
        
        !skip 
        do t=1,5
           iloc = index(line,",")
           line = line(iloc+1:len(line))
        enddo
        !tempc
        iloc = index(line,",")
        read(line(1:iloc-1),*,iostat=ios1) tempc
        line = line(iloc+1:len(line))
        !dewpc
        iloc = index(line,",")
        read(line(1:iloc-1),*, iostat=ios2) dewpc
        if(ios1.ne.0.and.ios2.eq.0) then 
           avgc = dewpc
        elseif(ios1.eq.0.and.ios2.ne.0) then
           avgc = tempc
        else
           avgc = (tempc+dewpc)/2
        endif
        line = line(iloc+1:len(line))

        !snod
        read(line,*) snod
        !average of dew point and air temperature is used to 
        !filter snow values. 
        if(snod.gt.990.or.avgc.gt.0) then 
           snod = LVT_rc%udef
        else           
           snod = snod/100 ! convert from cm to m. 
        endif

        call ESMF_TimeSet(datatime, yy=yy, mm=mm, dd=dd, h=hr, &
             calendar=LVT_calendar, rc=status)
        call LVT_verify(status, 'error in timeset: readSNODEPmetobs')
        
        call ESMF_TimeSet(starttime, yy=yy, mm=mm, dd=dd, h=0, & 
             calendar= LVT_calendar, rc=status)
        call LVT_verify(status, 'error in timeset: readSNODEPmetobs')
        
        data_index = (datatime - starttime)/SNODEPmetobs(source)%ts + 1
        
        call latlon_to_ij(LVT_domain%lvtproj, lat, lon, col, row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
             stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
           SNODEPmetobs(source)%snod(stn_col, stn_row, data_index) = snod
        endif
     enddo
  else
     SNODEPmetobs(source)%snod(:, :, :) = LVT_rc%udef
  endif

  call LVT_releaseUnitNumber(ftn)

end subroutine read_snodepdata
