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
! !ROUTINE: readCEOP
! \label{readCEOP}
!
! !INTERFACE: 
subroutine readCEOP(source)
! 
! !USES: 
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod
  use LVT_logMod,       only : LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, &
       LVT_logunit, LVT_endrun, LVT_verify
  use CEOP_obsMod,      only : CEOPobs
  use map_utils
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS:
  integer,   intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads the CEOP data (surface meteorology, 
!  flux, and soil moisture and soil temperature) data. 
!  At the start of the program, the entire data is read 
!  and stored based on the time location, and the station 
!  id it correspond to. At future times, the read routine simply 
!  indexes into the right location. Depending upon the frequency of 
!  computing output statistics, the routine also computes time
!  average (between different LIS output intervals). The routine also 
!  interpolates the soil moisture and soil temperature data to the model's 
!  vertical resolution. Note that CEOP soil layer structure is relatively 
!  fine compared to Noah's vertical structure. As a result, only the top 
!  soil layer is likely to contain relevant information that can be used 
!  for comparison.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! 
! 
! 
! !ARGUMENTS: 

! 
! 
!EOP
  type(ESMF_Time)           :: time1, time2, startTime
  type(ESMF_Time)           :: currTime
  type(ESMF_TimeInterval)   :: ts
  integer                   :: yr, mo,da,hr,mn,ss,doy
  real                      :: gmt
  real*8                    :: lis_prevtime
  integer                   :: status
  integer                   :: ftn
  integer                   :: nts,tindex
  real                      :: col,row,lat,lon
  integer                   :: stn_col, stn_row
  integer                   :: dim1ID, latId, lonId, timeId,tskinId
  integer                   :: tskinflagId
  real, allocatable             :: stn_tskin(:)
  character, allocatable      :: stn_tskin_flag(:)
  real, allocatable             :: time(:)
  integer                   :: t,k,st,et,c,r
  logical                   :: file_exists
  character*100             :: sfcfile
  real                      :: tskin(LVT_rc%lnc,LVT_rc%lnr)
  integer                   :: ntskin(LVT_rc%lnc,LVT_rc%lnr)

  if(ceopobs(source)%startFlag) then 
     ceopobs(source)%startFlag = .false. 
     
     call ESMF_TimeSet(startTime, &
          yy = 1970, mm = 1, dd=1,h=0,m=0,s=0,&
          calendar = LVT_calendar, & 
          rc = status)
     call LVT_verify(status)

     call ESMF_TimeSet(ceopobs(source)%startTime,&
          yy=LVT_rc%dyr(source), mm = LVT_rc%dmo(source), dd=LVT_rc%dda(source),&
          h=LVT_rc%dhr(source),m=LVT_rc%dmn(source),s=LVT_rc%dss(source),&
          calendar = LVT_calendar, & 
          rc = status)
     call LVT_verify(status)

     if(ceopobs(source)%readsfc.eq.1) then 

        do k=1,ceopobs(source)%nstns
           sfcfile = trim(ceopobs(source)%odir)//'/'//&
                trim(ceopobs(source)%campaign(k))//'_'//&
                trim(ceopobs(source)%stnname(k))//'_'//&
                trim(ceopobs(source)%locname(k))//'_sfc.nc'
           inquire(file=trim(sfcfile),exist=file_exists) 
           if(file_exists) then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)            
              write(LVT_logunit,*) '[INFO] Reading ',trim(sfcfile)
              call LVT_verify(nf90_open(path=sfcfile, &
                   mode=NF90_NOWRITE,ncid=ftn),&
                   'nf90_open failed for '//trim(sfcfile))
              
              call LVT_verify(nf90_inq_dimid(ftn,'time',&
                   dim1ID),'nf90_inq_dimid  failed for time')
              
              call LVT_verify(nf90_inquire_dimension(ftn,&
                   dim1Id,len=nts),&
                   'nf90_inquire_dimension failed for time')

              allocate(stn_tskin(nts))
!              allocate(stn_tskin_flag(nts))
              allocate(time(nts))

              call LVT_verify(nf90_inq_varid(ftn,&
                   'latitude',&
                   latId),'nf90_inq_varid failed for latitude')
              call LVT_verify(nf90_inq_varid(ftn,&
                   'longitude',&
                   lonId),'nf90_inq_varid failed for longitude')
              call LVT_verify(nf90_get_var(ftn,&
                   latId,lat), &
                   'nf90_inq_var failed for latitude')
              call LVT_verify(nf90_get_var(ftn,&
                   lonId,lon), &
                   'nf90_inq_var failed for longitude')

              call latlon_to_ij(LVT_domain%lvtproj, &
                   lat,lon, col,row)

              ceopobs(source)%stn(k)%col = nint(col)
              ceopobs(source)%stn(k)%row = nint(row)

              call LVT_verify(nf90_inq_varid(ftn,&
                   'time',&
                   timeId),'nf90_inq_varid failed for time')
              call LVT_verify(nf90_get_var(ftn,&
                   timeId,time),&
                   'nf90_inq_var failed for time')
              call LVT_verify(nf90_inq_varid(ftn,&
                   'surface_temperature',&
                   tskinId),'nf90_inq_varid failed for surface_temperature')
              call LVT_verify(nf90_get_var(ftn,&
                   tskinId,stn_tskin),&
                   'nf90_get_var failed for surface_temperature')
!              call LVT_verify(nf90_inq_varid(ftn,&
!                   'surface_temperature_flag',&
!                   tskinflagId),&
!                   'nf90_inq_varid failed for surface_temperature_flag')
!              call LVT_verify(nf90_get_var(ftn,&
!                   tskinflagId,stn_tskin_flag),&
!                   'nf90_get_var failed for surface_temperature_flag')

              do t=1,nts
                 call ESMF_TimeIntervalSet(ts,&
                      s = nint(time(t)),rc=status)
                 call LVT_verify(status)
                 
                 currTime = startTime + ts
                 
                 call ESMF_TimeGet(currTime,yy=yr,mm=mo,dd=da,&
                      h=hr,m=mn,s=ss) 
                 tindex = (currTime - ceopobs(source)%startTime)/ceopobs(source)%ts
                 if(tindex.gt.0.and.stn_tskin(t).gt.240) then 
                    ceopobs(source)%stn(k)%tskin(tindex) = stn_tskin(t)
                 endif
              enddo
              call LVT_verify(nf90_close(ftn))
              
              write(LVT_logunit,*) '[INFO] Finished processing ',trim(sfcfile)
              deallocate(time)
              deallocate(stn_tskin)
#endif
           endif
        enddo
     endif
  endif

  call ESMF_TimeSet(time1,&
       yy=LVT_rc%dyr(source), mm = LVT_rc%dmo(source), dd=LVT_rc%dda(source),&
       h=LVT_rc%dhr(source),m=LVT_rc%dmn(source),s=LVT_rc%dss(source),&
       calendar = LVT_calendar, & 
       rc = status)
  call LVT_verify(status)

  yr = LVT_rc%dyr(source)
  mo = LVT_rc%dmo(source)
  da = LVT_rc%dda(source)
  hr = LVT_rc%dhr(source)
  mn = LVT_rc%dmn(source)
  ss = LVT_rc%dss(source)
  
  call LVT_tick(lis_prevtime, doy, gmt, yr,mo,da,hr,mn,ss,(-1)*LVT_rc%ts)
  
  call ESMF_TimeSet(time2, yy=yr, &
       mm = mo, &
       dd = da, &
       h = hr, &
       m = mn, &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readARMobs')

  st = ((time2 - ceopobs(source)%startTime)/ceopobs(source)%ts) + 1
  et = ((time1 - ceopobs(source)%startTime)/ceopobs(source)%ts) + 1
  
  tskin = 0 
  ntskin = 0
  
  do k=1,ceopobs(source)%nstns
     if(st.gt.0.and.et.gt.0) then 
        do t=st,et
           if(ceopobs(source)%stn(k)%tskin(t).gt.0) then 
              tskin(ceopobs(source)%stn(k)%col,&
                   ceopobs(source)%stn(k)%row) = tskin(ceopobs(source)%stn(k)%col,&
                   ceopobs(source)%stn(k)%row) + ceopobs(source)%stn(k)%tskin(t)
              ntskin(ceopobs(source)%stn(k)%col,&
                   ceopobs(source)%stn(k)%row) = ntskin(ceopobs(source)%stn(k)%col,&
                   ceopobs(source)%stn(k)%row) +1
           endif
        enddo
     endif
  enddo

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(ntskin(c,r).gt.0) then 
           tskin(c,r) = tskin(c,r)/ntskin(c,r)
        else
           tskin(c,r) = -9999.0
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_avgsurft,source, tskin,vlevel=1,units="K")

end subroutine readCEOP

