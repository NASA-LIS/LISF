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
! !ROUTINE: readSNOTELObs
! \label{readSNOTELObs}
!
! !INTERFACE: 
subroutine readSNOTELObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_verify, LVT_endrun
  use LVT_timeMgrMod, only : LVT_calendar, LVT_tick
  use SNOTEL_obsMod,    only : snotelobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for SNOTEL station data. 
! LVT expects the data to be organized per calendar year, with 
! each file containing a daily data. Each reported observation is
! assumed to be time averaged. 
! 
! SWE (PILL) is recorded in inches
!
! At the start of the each year, the whole year's data is read and
! stored. At other times, LVT simply indexes into the stored
! arrays for retrieving the required values. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 AUG 2009: Sujay Kumar, Initial Specification
! 
!EOP

  integer                :: i,j,t,c,r,jj
  integer                :: stn_col, stn_row
  real                   :: col,row
  character*100          :: snotelfilename
  logical                :: file_exists
  logical                :: readflag
  integer                :: ftn, ios,ios1
  integer                :: yr,mo,da
  integer                :: status
  type(ESMF_Time)        :: snoteltime, snoteltime1
  integer                :: stnindex,tind
  real                   :: offset
  real                   :: swe_data, prcp_data
  character*100          :: line
  character*10           :: datestring
  integer                :: kk, iloc
  real                   :: swe(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: prcp(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: nswe(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: nprcp(LVT_rc%lnc,LVT_rc%lnr)

  swe   = 0 
  prcp  = 0
  nswe   = 0 
  nprcp  = 0

  if((LVT_rc%dyr(source).ne.snotelobs(source)%yr).or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     snotelobs(source)%yr = LVT_rc%dyr(source)
     snotelobs(source)%swe = LVT_rc%udef
     snotelobs(source)%prcp = LVT_rc%udef

     call ESMF_TimeSet(snotelobs(source)%startTime,  yy=snotelobs(source)%yr, &
          mm = LVT_rc%dmo(source), &
          dd = 1, &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'error in setting snotel start time')

     do jj=1,2 !repeat twice -- once to read the current year and one to read
        ! the previous -This is done because the snotel data is provided
        ! in water years. 
        do i=1,snotelobs(source)%nstns
           if(jj.eq.1) then 
              call create_SNOTEL_filename(snotelobs(source)%odir, &
                   snotelobs(source)%statename(i), &
                   snotelobs(source)%stnid(i),LVT_rc%dyr(source), snotelfilename)
           elseif(jj.eq.2) then 
              call create_SNOTEL_filename(snotelobs(source)%odir, &
                   snotelobs(source)%statename(i), &
                   snotelobs(source)%stnid(i),LVT_rc%dyr(source)+1, snotelfilename)
           endif

           inquire(file=trim(snotelfilename),exist=file_exists)

           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading SNOTEL data ',trim(snotelfilename)
              ftn=LVT_getNextUnitNumber()
              open(ftn,file=trim(snotelfilename),form='formatted')
              read(ftn,*) ! skip the header line

              readflag = .true. 
              do while(readflag) 
                 !              read(ftn,fmt='(i2.2,i2.2,i2.2,F5.2)', iostat=ios) mo,da,yr,swe_data
                 !              read(ftn,fmt='(i2.2,i2.2,i2.2,a)', iostat=ios) mo,da,yr,cswe_data
                 !              read(cswe_data,*,iostat=ios1) swe_data
                 read(ftn,'(a)',iostat=ios) line
                 if(ios.ne.0) then 
                    readflag = .false. 
                    exit
                 endif
                 do kk=1,len(line)
                    if(line(kk:kk)==achar(9)) line(kk:kk) =';'
                 enddo

                 iloc = index(line,";")
                 read(line(1:iloc-1),*) datestring
                 line = line(iloc+1:len(line))

                 read(datestring,'(i2.2,i2.2,i2.2)') mo,da,yr

                 iloc = index(line,";")
                 read(line(1:iloc-1),*,iostat=ios1) swe_data
                 line = line(iloc+1:len(line))
                 if(ios1.ne.0) swe_data = -9999.0

                 do kk=1,4 !skip
                    iloc = index(line,";")
                    line = line(iloc+1:len(line))
                 enddo

                 read(line,*,iostat=ios1) prcp_data
                 if(ios1.ne.0) prcp_data = -9999.0

                 if(ios.ne.0) then 
                    readflag = .false. 
                 else
                    if(yr.gt.60) then 
                       yr = yr + 1900 
                    else
                       yr = yr + 2000
                    endif

                    call ESMF_TimeSet(snoteltime,yy=yr,&
                         mm=mo, dd=da, h=0,calendar=LVT_calendar,&
                         rc=status)
                    call LVT_verify(status,'ESMF_TimeSet in readSnotelobs(Source)')

                    t = nint((snoteltime-snotelobs(source)%starttime)/&
                         snotelobs(source)%timestep)+1
                    if(t.ge.1.and.t.le.366) then 
                       if(swe_data.ge.0) then 
                          snotelobs(source)%swe(i,t) = swe_data*0.0254 !convert from inches to m
                       else
                          snotelobs(source)%swe(i,t) = LVT_rc%udef
                       endif

                       if(snotelobs(source)%swe(i,t).gt.20) then 
                          write(LVT_logunit,*) '[ERR]',i,t,snotelobs(source)%swe(i,t)
                          call LVT_endrun()
                       endif
                    endif

                    if(t.ge.1.and.t.le.366) then 
                       if(prcp_data.ge.0) then 
                          snotelobs(source)%prcp(i,t) = prcp_data*0.0254 !convert from inches to m
                       else
                          snotelobs(source)%prcp(i,t) = LVT_rc%udef
                       endif

                       if(snotelobs(source)%prcp(i,t).gt.20) then 
                          write(LVT_logunit,*) i,t,snotelobs(source)%prcp(i,t)
                          call LVT_endrun()
                       endif
                    endif
                 endif
                 if(ios.ne.0) then 
                    readflag = .false. 
                 endif
              enddo
              call LVT_releaseUnitNumber(ftn)
           endif
        enddo
     enddo
  endif

  call ESMF_TimeSet(snoteltime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'snoteltime1 set failed')

  offset = (snoteltime1-snotelobs(source)%starttime)/snotelobs(source)%timestep

  if((nint(offset)-offset).eq.0) then !only when LIS time matches the observation

     tind = nint((snoteltime1-snotelobs(source)%starttime)/snotelobs(source)%timestep)+1

     do i=1,snotelobs(source)%nstns
        call latlon_to_ij(LVT_domain%lvtproj, &
             snotelobs(source)%stnlat(i),snotelobs(source)%stnlon(i),&
             col,row)
        stn_col = nint(col)
        stn_row = nint(row)

        if(stn_col.gt.0.and.stn_row.gt.0.and.tind.ge.0.and.&
             stn_col.le.LVT_rc%lnc.and.stn_row.le.LVT_rc%lnr) then 
           if(snotelobs(source)%swe(i,tind).ne.LVT_rc%udef) then 
              swe(stn_col, stn_row) = swe(stn_col, stn_row) + & 
                   snotelobs(source)%swe(i,tind)
              nswe(stn_col, stn_row) = nswe(stn_col, stn_row) + 1
           endif
           if(snotelobs(source)%prcp(i,tind).ne.LVT_rc%udef) then 
              prcp(stn_col, stn_row) = prcp(stn_col, stn_row) + & 
                   snotelobs(source)%prcp(i,tind)
              nprcp(stn_col, stn_row) = nprcp(stn_col, stn_row) + 1
           endif
        endif
     enddo
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nswe(c,r).gt.0) then 
           swe(c,r) = swe(c,r)/nswe(c,r)
        else
           swe(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source, swe,vlevel=1,units="m")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nswe(c,r).gt.0) then 
           swe(c,r) =swe(c,r)*1000.0 !to mm
        else
           swe(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source, swe,vlevel=1,units="kg/m2")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nprcp(c,r).gt.0) then 
           prcp(c,r) = prcp(c,r)*1000.0/nprcp(c,r) !mm
        else
           prcp(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source, prcp,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source, prcp,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source, prcp,vlevel=1,units='kg/m2')

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(nprcp(c,r).gt.0) then 
           prcp(c,r) = prcp(c,r)/86400.0 !kg/m2s
        else
           prcp(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source, prcp,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source, prcp,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source, prcp,vlevel=1,units='kg/m2s') 

end subroutine readSNOTELObs

!BOP
! 
! !ROUTINE: create_SNOTEL_filename
! \label(create_SNOTEL_filename)
!
! !INTERFACE:
subroutine create_SNOTEL_filename(odir, stateid, stnid, yr,snotelname)
! 
! !USES:   
  use LVT_String_Utility
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: stateid
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: snotelname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the SNOTEL station
! 
!  The arguments are: 
!  \begin{description}
!   \item[stnid] Station ID 
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr

  write(fyr, '(i4.4)' ) yr

  snotelname = trim(odir)//'/'//trim(stateid)//'/'//LVT_StrLowCase(trim(stnid))//'_'&
       //trim(fyr)//'.tab'
  
end subroutine create_SNOTEL_filename
