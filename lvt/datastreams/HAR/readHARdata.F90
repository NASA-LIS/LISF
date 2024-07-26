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
! !INTERFACE: readHARdata
! \label{readHARdata}
! 
! !ROUTINE: 
subroutine readHARdata(source)
! !USES: 
   use ESMF
   use HAR_dataMod
   use LVT_coreMod
   use LVT_histDataMod
   use LVT_logMod
   use LVT_timeMgrMod
!
! !DESCRIPTION: 
! 
!   This subroutine reads the HAR data, processes it onto the LVT
!   grid and logs it. 
!
! !REVISION HISTORY: 
!  20 Apr 2018  Sujay Kumar  Initial Specification
! 
!EOP

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

   ! Defaults
   implicit none
   
   ! Arguments
   integer,intent(in) :: source

   ! Local variables
   integer             :: ftn 
   character*100       :: filename
   logical             :: file_exists
   real                :: prcp_in(hardata(source)%nc,hardata(source)%nr)
   real                :: prcp_in1(hardata(source)%nc*hardata(source)%nr)
   logical*1           :: lb(hardata(source)%nc*hardata(source)%nr)
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
   integer :: precipid
   integer :: start(3), count(3)
   integer :: jda2

   ! Initialize variables
   prcp(:) = LVT_rc%udef
   prcp_final(:,:) = LVT_rc%udef

   currTime = float(LVT_rc%dhr(source))*3600+ &
        60*LVT_rc%dmn(source) + LVT_rc%dss(source)

   alarmCheck = (mod(currtime,86400.0).eq.0)

   ! Read precipitation.  This will be 24-hour accumulation.
   ! We can divide this equally into sub-time periods. 
   yr1 = LVT_rc%dyr(source)
   mo1 = LVT_rc%dmo(source)
   da1 = LVT_rc%dda(source)
   hr1 = LVT_rc%dhr(source)
   mn1 = 0
   ss1 = 0
   call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
   call LVT_verify(status)

   ! If it is 00Z, use previous day's time level
   if (mod(currtime,86400.0).eq.0) then
      call ESMF_TimeIntervalSet(lis_ts, s = 86400, &
           rc=status)
      call LVT_verify(status)  
   else
      call ESMF_TimeIntervalSet(lis_ts, s = 0, &
           rc=status)
      call LVT_verify(status)  
   end if
   time2 = time1 - lis_ts

   call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
        h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, &
        dayOfYear=jda2, rc=status)
   call LVT_verify(status)

   if (alarmCheck) then
      
      call create_HAR_filename(hardata(source)%odir,yr2,filename)

      inquire(file=trim(filename),exist=file_exists)

      if(file_exists) then 
         write(LVT_logunit,*) '[INFO] Reading HAR data ',trim(filename)
         iret = nf90_open(path=trim(filename),mode=NF90_NOWRITE, &
              ncid=ftn)
         if (iret .eq. 0) then
            call LVT_verify(nf90_inq_varid(ftn,"prcp",precipid),&
                 'nf90_inq_varid failed for prcp')

            start(1) = 1
            start(2) = 1
            start(3) = jda2
            count(1) = hardata(source)%nc
            count(2) = hardata(source)%nr
            count(3) = 1
            call LVT_verify(nf90_get_var(ftn,precipid,prcp_in, start=start, &
                 count=count),&
                 'Error in nf90_get_var precip')
            call LVT_verify(nf90_close(ftn))

            prcp_in1(:) = 0
            lb(:) = .false.
            t = 1
            do r = 1, hardata(source)%nr
               do c = 1,hardata(source)%nc
                  if(.not.isNaN(prcp_in(c,r))) then 
                     prcp_in1(t) = prcp_in(c,r)/3600.0 !from mm/hr to mm/s
                     if (prcp_in1(t) .ge. 0) then
                        lb(t) = .true.
                     end if
                  endif
                  t = t + 1
               end do ! c
            end do ! r

            ! EMK...Use budget-bilinear interpolation if HAR is at 
            ! coarser resolution than the analysis grid; otherwise, use
            ! upscale averaging.
            if (LVT_isAtAFinerResolution(hardata(source)%datares)) then
               call conserv_interp(LVT_rc%gridDesc,lb,prcp_in1, &
                    lo,prcp, hardata(source)%nc*hardata(source)%nr, &
                    LVT_rc%lnc*LVT_rc%lnr, hardata(source)%rlat, &
                    hardata(source)%rlon, &
                    hardata(source)%w112, hardata(source)%w122, &
                    hardata(source)%w212, hardata(source)%w222, &
                    hardata(source)%n112, hardata(source)%n122, &
                    hardata(source)%n212, hardata(source)%n222, &
                    LVT_rc%udef, iret)     
            else
               call upscaleByAveraging(&
                    hardata(source)%nc*hardata(source)%nr, &
                    LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
                    hardata(source)%n11, lb, &
                    prcp_in1, lo, prcp)
            end if
            write(LVT_logunit,*) '[INFO] Finished processing ',trim(filename)
               
         else
            write(LVT_logunit,*)'[ERR] Read error with HAR file ', &
                 trim(filename)
            prcp = LVT_rc%udef
         endif
      else 
         write(LVT_logunit,*)'[ERR] Missing HAR file ', trim(filename)
         prcp = LVT_rc%udef
      end if ! file_exists
         
      do r=1,LVT_rc%lnr
         do c=1, LVT_rc%lnc
            prcp_final(c,r) = prcp(c+(r-1)*LVT_rc%lnc)
         end do ! c
      end do ! r

   end if ! alarmCheck

   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,&
        vlevel=1,units='kg/m2s')
   
   ! Now convert from kg/m2s to kg/m2
   do r=1,LVT_rc%lnr
      do c=1,LVT_rc%lnc
         if(prcp_final(c,r).ge.0) then 
            prcp_final(c,r) = prcp_final(c,r)*86400.0 !kg/m2
         else
            prcp_final(c,r) = LVT_rc%udef
         endif
      enddo ! c
   enddo ! r
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final,&
        vlevel=1,units='kg/m2') 
      
end subroutine readHARdata

!------------------------------------------------------------------------------

!BOP
! !ROUTINE: create_HAR_filename
! \label{create_HAR_filename}
!
! !INTERFACE: 
subroutine create_HAR_filename(odir,yr,filename)

   implicit none
! !ARGUMENTS:
   character(len=*), intent(in) :: odir
   integer, intent(in) :: yr
   character(len=*), intent(out) :: filename
! 
! !DESCRIPTION: 
! This subroutine creates the name of the HAR data files
!
!EOP
   ! Local variables
   character*4             :: fyr

   write(fyr, '(i4.4)' ) yr

   filename  = trim(odir)//'/har_d10km_d_2d_prcp_'//trim(fyr)//'.nc'

end subroutine create_HAR_filename
