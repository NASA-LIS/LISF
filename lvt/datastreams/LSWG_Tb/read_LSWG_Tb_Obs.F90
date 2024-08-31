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
! !ROUTINE: read_LSWG_Tb_Obs
! \label{read_LSWG_Tb_Obs}
!
! !INTERFACE: 
subroutine read_LSWG_Tb_Obs(k)
! 
! !USES:   
  use ESMF

#if 0   
  use LVT_coreMod,    only : LVT_rc, LVT_domain
  use LVT_RTMhistDataMod,  only : LVT_logRTMSingleVar
  use LVT_RTMhistDatamod, only : LVT_RTMhistData
  use LVT_RTMobsDataMod,   only : LVT_RTMobsData
  use LVT_logMod,     only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_verify, LVT_endrun
  use LVT_timeMgrMod, only : LVT_calendar, LVT_tick
  use LSWG_Tb_obsMod,    only : lswg_Tbobs
  use map_utils
#endif
  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: k   !why argument?
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for LSWGTb station data. 
! LVT expects the data to be organized per calendar year, with 
! each file containing a daily data. Each reported observation is
! assumed to be time averaged. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 AUG 2009: Sujay Kumar, Initial Specification
! 
!EOP

  integer             :: ftn
  integer             :: kk, rtm_kk !kk because of k argument of this fxn call?
  integer             :: i,j,c,r,t,t_prev,offset
  integer             :: ios
  logical             :: file_exists
  logical             :: readflag
  character*4         :: sname
  real                :: col,row
  integer             :: stn_col, stn_row
  integer             :: yr, mo, da, gmt_hr, gmt_mn, gmt_sc, hr
  real                :: lat, lon, z_angle, rain_rate
  real, allocatable   :: tb(:) !data channel space by record
  integer, dimension(13) :: cloud

  !gridded Tb running Sum in data channel space
  real,   allocatable :: Tb1(:,:,:,:)  
  !gridded num obs in data channel space
  integer,allocatable :: nTb1(:,:,:,:) 
  !gridded Tb in RTM channel space to which to apply cloud mask
  real,   allocatable :: Tb2(:,:,:,:)  

  character*100       :: maskfile
  character*10        :: ftime
  type(ESMF_Time)     :: Tbtime, masktime, masktime_prev, temptime
  type(ESMF_TimeInterval) :: timestep, mask_timestep,doy
  integer             :: status
  integer             :: myyy, myd, mymo,myh,mym,mys
  integer             :: counter
 
#if 0  
if(lswg_Tbobs%start) then 
     allocate(tb(lswg_Tbobs%numchannels))
     allocate(Tb1(LVT_rc%lnc,LVT_rc%lnr, lswg_Tbobs%nts, lswg_Tbobs%numchannels))
     allocate(nTb1(LVT_rc%lnc,LVT_rc%lnr, lswg_Tbobs%nts,lswg_Tbobs%numchannels))
     Tb1 = 0.0
     nTb1 = 0
     ftn = LVT_getNextUnitNumber()
     inquire(file=trim(lswg_Tbobs%filename), exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading ',trim(lswg_Tbobs%filename)
        open(ftn,file=trim(lswg_Tbobs%filename),form='formatted')
        readflag = .true. 
        counter=0
        if (lswg_Tbobs%dataformat.eq.0) then 
           do while(readflag) 
              read(ftn,200,iostat=ios) sname, yr, mo, da, gmt_hr, gmt_mn, &
                   lat, lon, z_angle, rain_rate, &
                   (tb(i), i=1,lswg_Tbobs%numchannels)
!if (counter .le. 10) print*,counter, sname, yr, mo, da, gmt_hr, gmt_mn, &
!                   lat, lon, z_angle, rain_rate, &
!                   tb
              if(ios.ne.0) then
                 readflag = .false. 
              else
                 if(yr.gt.0.and.(trim(sname).eq.trim(lswg_Tbobs%sname))) then 
                    call ESMF_TimeSet(Tbtime, yy=yr,&
                         mm=mo, dd=da,h=gmt_hr,m=0,&
                         calendar=LVT_calendar, rc=status)
                    call LVT_verify(status, 'ESMF_timeset error in readLSWGObs')
                    if(gmt_mn.le.30) then !approximate to the nearest last hour
                       Tbtime = Tbtime
                    else
                       call ESMF_TimeIntervalSet(timestep, s=3600, rc=status)
                       call LVT_verify(status, 'ESMF_timeset error in readLSWGObs')
                       Tbtime = Tbtime + timestep
                    endif
                    t = ((Tbtime - lswg_Tbobs%startTime)/lswg_Tbobs%timestep)+1
                    call latlon_to_ij(LVT_domain%lvtproj, lat,lon, col,row)
                    stn_col = nint(col)
                    stn_row = nint(row)                 
                    if((stn_row.gt.0.and.stn_col.gt.0).and.&
                         (stn_row.le.LVT_rc%lnr.and.stn_col.le.LVT_rc%lnc)) then 
                       !consider a 0.125 dx box around the LIS grid point
                       !   -->1/2 of grid cell-->within 0.25 of center of grid cell
                    counter=counter+1
if (counter .le. 10) print*,counter, sname, yr, mo, da, gmt_hr, gmt_mn, &
                   lat, lon, z_angle, rain_rate, &
                   tb
                       if((abs(stn_col - col).le.0.25).and.&
                            (abs(stn_row - row).le.0.25)) then 
                          do kk=1,lswg_Tbobs%numchannels
                             if(tb(kk).gt.0.and.t.ge.1) then
                                Tb1(stn_col, stn_row,t,kk) = &
                                     Tb1(stn_col, stn_row,t,kk) + tb(kk)
                                nTb1(stn_col, stn_row,t,kk) = nTb1(stn_col, stn_row,t,kk) + 1
                             endif
                          enddo
                       endif
                    endif
                 endif
              endif
           enddo
        elseif (lswg_Tbobs%dataformat.eq.1) then 
           do while(readflag) 
              read(ftn,300,iostat=ios) lat, lon, da, yr, gmt_hr, gmt_mn, gmt_sc, &
                   (tb(i), i=1,lswg_Tbobs%numchannels)
              if(ios.ne.0) then
                 readflag = .false. 
              else
                 if(yr.gt.0) then 
                    call ESMF_TimeSet(Tbtime, yy=yr,h=0, &
                         calendar=LVT_calendar, rc=status)
                    da=da-1
                    call ESMF_TimeIntervalSet(doy, d = da, h=gmt_hr, m=0,rc=status)
                    call LVT_verify(status, 'ESMF_timeset error in readLSWGObs')
                    Tbtime=Tbtime+doy
                    call ESMF_TimeIntervalSet(timestep, s=3600, rc=status)
                    call LVT_verify(status, 'ESMF_timeset error in readLSWGObs')
                    if(gmt_mn.le.30) then !approximate to the nearest last hour
                       Tbtime = Tbtime
                    else
                       Tbtime = Tbtime + timestep
                    endif
                    t = ((Tbtime - lswg_Tbobs%startTime)/lswg_Tbobs%timestep)+1
                    call latlon_to_ij(LVT_domain%lvtproj, lat,lon, col,row)
                    stn_col = nint(col)
                    stn_row = nint(row)                 
                    if((stn_row.gt.0.and.stn_col.gt.0).and.&
                         (stn_row.le.LVT_rc%lnr.and.stn_col.le.LVT_rc%lnc)) then 
                       !consider a 0.125 dx box around the LIS grid point
                       !   -->1/2 of grid cell-->within 0.25 of center of grid cell
                       if((abs(stn_col - col).le.0.25).and.&
                            (abs(stn_row - row).le.0.25)) then 
                          do kk=1,lswg_Tbobs%numchannels
                             if(tb(kk).gt.0.and.t.ge.1) then 
                                Tb1(stn_col, stn_row,t,kk) = &
                                     Tb1(stn_col, stn_row,t,kk) + tb(kk)
                                nTb1(stn_col, stn_row,t,kk) = nTb1(stn_col, stn_row,t,kk) + 1
                             endif
                          enddo
                       endif
                    endif
                 endif
              endif
           enddo
        endif
200     format(a4,i5,i2.2,i2.2,i3,i2,2f9.2,f6.1,f7.2,20f6.1)
300     format(2f8.2,5i5,12f10.2)

        call LVT_releaseUnitNumber(ftn)
     else
        write(LVT_logunit,*) '[ERR] LSWG file ',lswg_Tbobs%filename, 'does not exist'
        write(LVT_logunit,*) '[ERR] Program stopping ... '
        call LVT_endrun
     endif

     !Compute average gridded Tb;store in lswg data structure
     do kk=1,lswg_Tbobs%numchannels
        do t=1, lswg_Tbobs%nts
           do c=1,LVT_rc%lnc
              do r=1,LVT_rc%lnr
                 if(nTb1(c,r,t,kk).gt.0) then 
                    lswg_Tbobs%Tb(c,r,t,kk) = Tb1(c,r,t,kk)/nTb1(c,r,t,kk)
                 else
                    !shouldn't be necessary as intialized in obsMod
                    !lswg_Tbobs%Tb(c,r,t,kk) = LVT_rc%udef 
                 endif
              enddo
           enddo
        enddo
     enddo

     deallocate(Tb1)
     deallocate(nTb1)

!Fill cloud mask data structure
     if(lswg_Tbobs%mask_option.eq.1) then 
! assume that the cloud masking is specified at 3 hr intervals
! the cloud masking data is assumed to be back averaged. we apply the 
! same cloud mask for the 3 hour window prior to the specified time. 

        call ESMF_TimeIntervalSet(mask_timestep, s=10800, rc=status)	
        call LVT_verify(status, 'ESMF_timeset error in readLSWGObs')

        ftn = LVT_getNextUnitNumber()
        inquire(file=trim(lswg_Tbobs%maskfile), exist=file_exists)
        if(file_exists) then  
           write(LVT_logunit,*) '[INFO] Reading cloud mask ',trim(lswg_Tbobs%maskfile)
           open(ftn,file=trim(lswg_Tbobs%maskfile),form='formatted')
           read(ftn,*)
           ios = 0
           do while(ios.eq.0)
              read(ftn,*,iostat=ios) yr, mo, da, hr, cloud
              call ESMF_TimeSet(masktime, yy=yr,&
                   mm=mo, dd=da,h=hr,m=0,&
                   calendar = LVT_calendar, rc=status)
              call LVT_verify(status, 'ESMF_timeset error in readLSWGObs')
              t = ((masktime -lswg_Tbobs%startTime)/lswg_Tbobs%timestep)+1
	      
	      masktime_prev = masktime -mask_timestep
              t_prev = nint((masktime_prev -lswg_Tbobs%startTime)/lswg_Tbobs%timestep)+1

              if(t_prev.le.0) t_prev = 1
              if(t.le.0) t = 1
              if(t.gt.t_prev) then 
                 t_prev = t_prev + 1 !move into the 3hr window
              endif
              if(cloud(lswg_Tbobs%maskcol).gt.lswg_Tbobs%cloud_pct) then 
                 lswg_Tbobs%mask(:,:,t_prev:t) = LVT_rc%udef
              endif
           enddo
        endif
        call LVT_releaseUnitNumber(ftn)
     endif

     lswg_Tbobs%start = .false.
  endif

  allocate(Tb2(LVT_rc%lnc,LVT_rc%lnr,1,LVT_RTMhistData%Tb%vlevels))
  Tb2=LVT_rc%udef
!Apply the cloud mask
  do kk=1,lswg_Tbobs%numchannels
     rtm_kk = lswg_Tbobs%data2rtm_channelmap(kk)
!print*, rtm_kk
     call ESMF_TimeSet(Tbtime, yy=LVT_rc%dyr(source), &
          mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
          s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
     call LVT_verify(status, 'ESMF_timeset error in readLSWGObs')
  
! is a time offset
     offset = nint((Tbtime - lswg_Tbobs%starttime)/lswg_Tbobs%timestep)+1
     if(offset.gt.0) then 
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(lswg_Tbobs%mask(c,r,offset).ne.LVT_rc%udef) then 
                 Tb2(c,r,1,rtm_kk) = lswg_Tbobs%Tb(c,r,offset,kk)
              endif
           enddo
        enddo
     endif
! is LVT_RTMobsData%Tb_obs initialized to LVT_rc%udef?????????
     call LVT_logRTMSingleVar(LVT_RTMobsData%Tb_obs, source,&
          Tb2(:,:,1,rtm_kk), vlevel =rtm_kk)
  enddo
  deallocate(Tb2)
#endif
end subroutine read_LSWG_Tb_Obs

!!$subroutine create_LSWGTb_filename(odir, stateid, stnid, yr, snotelname)
!!$  
!!$  implicit none
!!$
!!$! !ARGUMRENTS: 
!!$  character(len=*), intent(in)  :: odir
!!$  character(len=*), intent(in)  :: stateid
!!$  character(len=*), intent(in)  :: stnid
!!$  integer,          intent(in)  :: yr
!!$  character(len=*), intent(out) :: snotelname
!!$! !DESCRIPTION: 
!!$! 
!!$! This routine creates a filename for the LSWGTb station
!!$! 
!!$!  The arguments are: 
!!$!  \begin{description}
!!$!   \item[stnid] Station ID 
!!$!  \end{description}
!!$!EOP
!!$  character*4             :: fyr
!!$  
!!$  write(fyr, '(i4.4)' ) yr
!!$
!!$  snotelname = trim(odir)//'/'//trim(stateid)//'/'//trim(stnid)//'_'&
!!$       //trim(fyr)//'.dat'
!!$  
!!$end subroutine create_LSWGTb_filename
