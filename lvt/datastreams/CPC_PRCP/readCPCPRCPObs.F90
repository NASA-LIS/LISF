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
! !ROUTINE: readCPCPRCPobs
! \label{readCPCPRCPobs}
!
! !INTERFACE: 
subroutine readCPCPRCPobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_verify, LVT_endrun
  use LVT_timeMgrMod, only : LVT_calendar, LVT_tick
  use CPCPRCP_obsMod,    only : cpcprcpobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,  intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! TODO : 
!  * filename change
!  * check global data.  
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 AUG 2009: Sujay Kumar, Initial Specification
! 
!EOP

  integer             :: ftn 
  character*100       :: filename
  logical             :: file_exists
  real                :: prcp_in(cpcprcpobs(source)%nc*cpcprcpobs(source)%nr)
  real                :: prcp_in1(cpcprcpobs(source)%nc*cpcprcpobs(source)%nr)
  logical*1           :: lb(cpcprcpobs(source)%nc*cpcprcpobs(source)%nr)
  logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                :: prcp(LVT_rc%lnc*LVT_rc%lnr)
  real                :: prcp_final(LVT_rc%lnc, LVT_rc%lnr)
  integer             :: t,c,r
  integer             :: iret
  real                :: currTime
  logical             :: alarmCheck
  integer                      :: yr1, mo1, da1, hr1, mn1, ss1
  integer                      :: yr2, mo2, da2, hr2, mn2, ss2
  type(ESMF_Time)              :: time1
  type(ESMF_Time)              :: time2
  type(ESMF_TimeInterval)      :: lis_ts
  integer :: status

  prcp = LVT_rc%udef
  prcp_final = LVT_rc%udef

  ! New logic to prevent reads for sub-daily time periods.  Only read in at 
  ! 00Z.
  currTime = float(LVT_rc%dhr(source))*3600+ &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmCheck = (mod(currtime,86400.0).eq.0)

  ! Read the previous day's file.  This will be 24-hour accumulation ending
  ! at current time.
  yr1 = LVT_rc%dyr(source)
  mo1 = LVT_rc%dmo(source)
  da1 = LVT_rc%dda(source)
  hr1 = LVT_rc%dhr(source)
  mn1 = 0
  ss1 = 0
  call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar, rc=status)
  call LVT_verify(status)
  
  call ESMF_TimeIntervalSet(lis_ts, s = 86400, &
       rc=status)
  call LVT_verify(status)  
  time2 = time1 - lis_ts

  call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
        h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, rc=status)
  call LVT_verify(status)

  if (alarmCheck) then

     call create_CPCPRCP_filename(cpcprcpobs(source)%odir, &
          yr2,mo2,da2, &
          cpcprcpobs(source)%domainType, &
          cpcprcpobs(source)%userealtime, &
          filename)
  
     inquire(file=trim(filename),exist=file_exists)
     
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading CPC data ',trim(filename)
        ftn = LVT_getNextUnitNumber()
        open(ftn,file=trim(filename),form='unformatted',access='direct',&
             convert='little_endian',&
             recl=cpcprcpobs(source)%nc*cpcprcpobs(source)%nr*4)
        read(ftn,rec=1) prcp_in1
        call LVT_releaseUnitNumber(ftn)


        if(cpcprcpobs(source)%domainType.eq."GLOBAL") then 
           
           do r=1,cpcprcpobs(source)%nr
              do c=1,cpcprcpobs(source)%nc
                 if(c.ge.360) then
                    prcp_in((c-360+1)+(r-1)*cpcprcpobs(source)%nc) = &
                         prcp_in1(c+(r-1)*cpcprcpobs(source)%nc)
                 else
                    prcp_in((c+360-1)+(r-1)*cpcprcpobs(source)%nc) = &
                         prcp_in1(c+(r-1)*cpcprcpobs(source)%nc)
                 endif
              enddo
           enddo
        else
           prcp_in = prcp_in1
        endif
        
        
        lb = .false.
        
        do t=1, cpcprcpobs(source)%nc*cpcprcpobs(source)%nr
           if(prcp_in(t).ge.0) then
              prcp_in(t) = prcp_in(t)*0.1 !mm/day
              lb(t) = .true. 
           endif
        enddo
           
        call bilinear_interp(LVT_rc%gridDesc,lb,prcp_in, &
             lo,prcp, cpcprcpobs(source)%nc*cpcprcpobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr, cpcprcpobs(source)%rlat, &
             cpcprcpobs(source)%rlon, cpcprcpobs(source)%w11,  cpcprcpobs(source)%w12, &
             cpcprcpobs(source)%w21,  cpcprcpobs(source)%w22,  cpcprcpobs(source)%n11, &
             cpcprcpobs(source)%n12,  cpcprcpobs(source)%n21,  cpcprcpobs(source)%n22, &
             LVT_rc%udef, iret)     
        
        write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)
     else
        prcp = LVT_rc%udef
     endif

     do r=1,LVT_rc%lnr
        do c=1, LVT_rc%lnc
           prcp_final(c,r) = prcp(c+(r-1)*LVT_rc%lnc)
        enddo
     enddo

  end if  

  ! Convert mm/day to kg/m2s (note that 1 mm = 1 kg/m2)
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(prcp_final(c,r).ge.0) then 
           prcp_final(c,r) = prcp_final(c,r)/86400.0 !kg/m2s
        else
           prcp_final(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source,prcp_final,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source,prcp_final,vlevel=1,units='kg/m2s')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,vlevel=1,units='kg/m2s')
     
  ! Now convert from kg/m2s to kg/m2
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(prcp_final(c,r).ge.0) then 
           prcp_final(c,r) = prcp_final(c,r)*86400.0 !kg/m2
        else
           prcp_final(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  
  call LVT_logSingleDataStreamVar(LVT_MOC_rainf,source,prcp_final,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_rainfforc,source,prcp_final,vlevel=1,units='kg/m2')
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,vlevel=1,units='kg/m2') 
     
end subroutine readCPCPRCPobs

!BOP
! 
! !ROUTINE: create_CPCPRCP_filename
! \label(create_CPCPRCP_filename)
!
! !INTERFACE:
subroutine create_CPCPRCP_filename(odir, yr,mo, da, domaintype, realtime, &
     filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer         , intent(in)  :: yr
  integer         , intent(in)  :: mo
  integer         , intent(in)  :: da
  character(len=*), intent(in)  :: domaintype
  integer         , intent(in)  :: realtime
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: filename
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the CPCPRCP station
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
  character*2             :: fmo
  character*2             :: fda

  write(fyr, '(i4.4)' ) yr
  write(fmo, '(i2.2)' ) mo
  write(fda, '(i2.2)' ) da

  if(domaintype.eq."CONUS") then 
     if(realtime.eq.1) then 
        filename = trim(odir)//'/'//trim(fyr)//&
             '/PRCP_CU_GAUGE_V1.0CONUS_0.25deg.lnx.'//&
             trim(fyr)//trim(fmo)//trim(fda)//'.RT'
     else
        filename = trim(odir)//'/'//trim(fyr)//&
             '/PRCP_CU_GAUGE_V1.0CONUS_0.25deg.lnx.'//&
             trim(fyr)//trim(fmo)//trim(fda)
     endif
  elseif(domaintype.eq."GLOBAL") then 
     if(realtime.eq.1) then 
        filename = trim(odir)//'/'//trim(fyr)//&
             '/PRCP_CU_GAUGE_V1.0GLB_0.50deg.lnx.'//&
             trim(fyr)//trim(fmo)//trim(fda)//'.RT'
     else
        filename = trim(odir)//'/'//trim(fyr)//&
             '/PRCP_CU_GAUGE_V1.0GLB_0.50deg.lnx.'//&
             trim(fyr)//trim(fmo)//trim(fda)
     endif
  endif
  
end subroutine create_CPCPRCP_filename


