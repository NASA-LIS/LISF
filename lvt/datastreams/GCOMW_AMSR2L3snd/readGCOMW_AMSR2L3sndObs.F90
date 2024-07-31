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
! !ROUTINE: readGCOMW_AMSR2L3sndObs
! \label{readGCOMW_AMSR2L3sndObs}
!
! !INTERFACE: 
subroutine readGCOMW_AMSR2L3sndObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use GCOMW_AMSR2L3snd_obsMod, only : GCOMW_AMSR2L3sndobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!  11 Dec 2014: Sujay Kumar, Initial Specification
! 
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j,status
  character*100     :: fname_A, fname_D
  integer           :: yr, mo, da, hr, mn, ss
  real              :: sndobs(LVT_rc%lnc, LVT_rc%lnr)
  real              :: sndobs_t(LVT_rc%lnc*LVT_rc%lnr)
  real              :: sndobs_av(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: count_av(LVT_rc%lnc, LVT_rc%lnr)
  type(ESMF_Time)   :: cTime, stTime, enTime

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LVT_rc%dhr(source))*3600 + &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  sndobs= LVT_rc%udef
        
  if(GCOMW_AMSR2L3sndobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     GCOMW_AMSR2L3sndobs(source)%startmode = .false. 

     if(LVT_rc%smoothObs.eq.1) then 
        sndobs_av = 0.0
        count_av = 0
        call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h  = LVT_rc%dhr(source), &
             m  = LVT_rc%dmn(source), & 
             s  = LVT_rc%dss(source), &
             calendar = LVT_calendar, & 
             rc = status)
        call LVT_verify(status)
        
        stTime = cTime - LVT_obsSmTwL
        enTime = cTime + LVT_obsSmTwL

        cTime = stTime
        do while (cTime.le.enTime) 
           call ESMF_TimeGet(cTime, yy = yr, mm=mo, dd=da,&
                h = hr, m = mn, s =ss, calendar=LVT_calendar, &
                rc=status)
           call LVT_verify(status)

           call create_GCOMW_AMSR2L3snd_A_filename(GCOMW_AMSR2L3sndobs(source)%odir, &
                yr, mo, da, fname_A)
           
           call create_GCOMW_AMSR2L3snd_D_filename(GCOMW_AMSR2L3sndobs(source)%odir, &
                yr, mo, da, fname_D)

           sndobs_t = LVT_rc%udef
           call read_amsr2_snd_data(source, fname_A, fname_D,sndobs_t)
           
           cTime = cTime + LVT_obsSmTwI

           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(sndobs_t(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                    sndobs_av(c,r) = sndobs_av(c,r) + sndobs_t(c+(r-1)*LVT_rc%lnc)
                    count_av(c,r) = count_av(c,r) + 1
                 endif
              enddo
           enddo
        enddo

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(count_av(c,r).gt.0) then 
                 sndobs(c,r) = sndobs_av(c,r)/count_av(c,r)
              else
                 sndobs(c,r) = LVT_rc%udef
              endif
           enddo
        enddo

     else

        call create_GCOMW_AMSR2L3snd_A_filename(GCOMW_AMSR2L3sndobs(source)%odir, &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), fname_A)
        
        call create_GCOMW_AMSR2L3snd_D_filename(GCOMW_AMSR2L3sndobs(source)%odir, &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), fname_D)
        
        call read_amsr2_snd_data(source, fname_A, fname_D,sndobs_t)

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              sndobs(c,r) = sndobs_t(c+(r-1)*LVT_rc%lnc)
           enddo
        enddo
     endif

  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_snowdepth, source,&
       sndobs,vlevel=1,units="m")
 
end subroutine readGCOMW_AMSR2L3sndObs


!BOP
! 
! !ROUTINE: read_amsr2_snd_data
! \label(read_amsr2_snd_data)
!
! !INTERFACE:
subroutine read_amsr2_snd_data(source, fname_A, fname_D, sndobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LVT_coreMod
  use LVT_logMod
  use map_utils,    only : latlon_to_ij
  use GCOMW_AMSR2L3snd_obsMod, only : GCOMW_AMSR2L3sndobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source 
  character (len=*)             :: fname_A
  character (len=*)             :: fname_D
  real                          :: sndobs_ip(LVT_rc%lnc*LVT_rc%lnr)


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname_A]        name of the GCOMW_AMSR2 AMSR-E file
!  \item[sndobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer             :: snd(1,GCOMW_AMSR2L3sndobs(source)%amsr2nc,GCOMW_AMSR2L3sndobs(source)%amsr2nr)
  integer           :: time(GCOMW_AMSR2L3sndobs(source)%amsr2nc,GCOMW_AMSR2L3sndobs(source)%amsr2nr)
  integer           :: sndtime(GCOMW_AMSR2L3sndobs(source)%amsr2nc,GCOMW_AMSR2L3sndobs(source)%amsr2nr)
  real                        :: snd_combined(GCOMW_AMSR2L3sndobs(source)%amsr2nc,GCOMW_AMSR2L3sndobs(source)%amsr2nr)
  real                        :: snd_data(GCOMW_AMSR2L3sndobs(source)%amsr2nc*GCOMW_AMSR2L3sndobs(source)%amsr2nr)
  logical*1                   :: snd_data_b(GCOMW_AMSR2L3sndobs(source)%amsr2nc*GCOMW_AMSR2L3sndobs(source)%amsr2nr)
  logical*1                   :: sndobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  logical                     :: file_exists
  integer                     :: c1,c,r,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid
  integer                     :: sndid, timeid
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  sndtime = -1.0
  snd_combined = LVT_rc%udef

  inquire(file=fname_A, exist=file_exists) 
  if(file_exists) then 
     write (LVT_logunit, *) "[INFO] Reading "//trim(fname_A)
     ios = nf90_open(path=trim(fname_A),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error opening file '//trim(fname_A))
     
     ios = nf90_inq_varid(nid, "Geophysical Data",sndid)
     call LVT_verify(ios, 'Error nf90_inq_varid: snd')
     
     ios = nf90_inq_varid(nid, "Time Information",timeid)
     call LVT_verify(ios, 'Error nf90_inq_varid: time')
  
     !values
     ios = nf90_get_var(nid, sndid, snd)
     call LVT_verify(ios, 'Error nf90_get_var: snd')
     
     ios = nf90_get_var(nid, timeid,time)
     call LVT_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LVT_verify(ios,'Error closing file '//trim(fname_A))

     do r=1, GCOMW_AMSR2L3sndobs(source)%amsr2nr
        do c=1, GCOMW_AMSR2L3sndobs(source)%amsr2nc
           c1 = c+GCOMW_AMSR2L3sndobs(source)%amsr2nc/2-1
           if(c1.gt.GCOMW_AMSR2L3sndobs(source)%amsr2nc) then 
              c1 = c-GCOMW_AMSR2L3sndobs(source)%amsr2nc/2-1
           endif
           
           if(time(c,r).lt.0.or.snd(1,c,r).le.0.001) then
              snd_combined(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1) = LVT_rc%udef
           else
              snd_combined(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1) = &
                   (snd(1,c,r)*0.1)/100.0
              sndtime(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1) = time(c,r)
           endif
        enddo
     enddo

  endif

  inquire(file=fname_D, exist=file_exists) 
  if(file_exists) then 
     write (LVT_logunit, *) "[INFO] Reading "//trim(fname_D)
     ios = nf90_open(path=trim(fname_D),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error opening file '//trim(fname_d))
     
     ios = nf90_inq_varid(nid, "Geophysical Data",sndid)
     call LVT_verify(ios, 'Error nf90_inq_varid: snd')
     
     ios = nf90_inq_varid(nid, "Time Information",timeid)
     call LVT_verify(ios, 'Error nf90_inq_varid: time')
  
     !values
     ios = nf90_get_var(nid, sndid, snd)
     call LVT_verify(ios, 'Error nf90_get_var: snd')
     
     ios = nf90_get_var(nid, timeid,time)
     call LVT_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LVT_verify(ios,'Error closing file '//trim(fname_d))


     do r=1, GCOMW_AMSR2L3sndobs(source)%amsr2nr
        do c=1, GCOMW_AMSR2L3sndobs(source)%amsr2nc
           c1 = c+GCOMW_AMSR2L3sndobs(source)%amsr2nc/2-1
           if(c1.gt.GCOMW_AMSR2L3sndobs(source)%amsr2nc) then 
              c1 = c-GCOMW_AMSR2L3sndobs(source)%amsr2nc/2-1
           endif
           if(snd_combined(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1).gt.0) then 
              !if valid data already exists, then overwrite only 
              !if a newer data is available. 
              if((time(c,r).gt.0).and.&
                   (time(c,r).gt.sndtime(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1)).and.&
                   (snd(1,c,r).gt.0.0)) then 
                 snd_combined(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1) = &
                      (snd(1,c,r)*0.1)/100.0
              end if
           else
              if(time(c,r).lt.0.or.snd(1,c,r).le.0.0) then
                 snd_combined(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1) = LVT_rc%udef
              else
                 snd_combined(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1) = &
                      (snd(1,c,r)*0.1)/100.0
                 sndtime(c1,GCOMW_AMSR2L3sndobs(source)%amsr2nr-r+1) = time(c,r)
              endif
              
           endif
        enddo
     enddo

  endif

  do r=1, GCOMW_AMSR2L3sndobs(source)%amsr2nr
     do c=1, GCOMW_AMSR2L3sndobs(source)%amsr2nc
        snd_data(c+(r-1)*GCOMW_AMSR2L3sndobs(source)%amsr2nc) = snd_combined(c,r)
        if(snd_combined(c,r).ne.LVT_rc%udef) then 
           snd_data_b(c+(r-1)*GCOMW_AMSR2L3sndobs(source)%amsr2nc) = .true. 
        else
           snd_data_b(c+(r-1)*GCOMW_AMSR2L3sndobs(source)%amsr2nc) = .false.
        endif
     enddo
  enddo

  if(LVT_isatAfinerResolution(GCOMW_AMSR2L3sndobs(source)%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LVT_rc%gridDesc(:),&
          snd_data_b, snd_data, sndobs_b_ip, sndobs_ip, &
          GCOMW_AMSR2L3sndobs(source)%amsr2nc*GCOMW_AMSR2L3sndobs(source)%amsr2nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          GCOMW_AMSR2L3sndobs(source)%rlat, GCOMW_AMSR2L3sndobs(source)%rlon, &
          GCOMW_AMSR2L3sndobs(source)%w11, GCOMW_AMSR2L3sndobs(source)%w12, &
          GCOMW_AMSR2L3sndobs(source)%w21, GCOMW_AMSR2L3sndobs(source)%w22, &
          GCOMW_AMSR2L3sndobs(source)%n11, GCOMW_AMSR2L3sndobs(source)%n12, &
          GCOMW_AMSR2L3sndobs(source)%n21, GCOMW_AMSR2L3sndobs(source)%n22, &
          LVT_rc%udef, ios)
  else
     call upscaleByAveraging(&
          GCOMW_AMSR2L3sndobs(source)%amsr2nc*GCOMW_AMSR2L3sndobs(source)%amsr2nr,&
          LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef, &
          GCOMW_AMSR2L3sndobs(source)%n11,snd_data_b, snd_data, sndobs_b_ip,sndobs_ip)
     
  endif

#endif
  
end subroutine read_amsr2_snd_data

!BOP
! !ROUTINE: create_GCOMW_AMSR2L3snd_A_filename
! \label{create_GCOMW_AMSR2L3snd_A_filename}
! 
! !INTERFACE: 
subroutine create_GCOMW_AMSR2L3snd_A_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GCOMW_AMSR2 filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GCOMW_AMSR2 soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated GCOMW_AMSR2 filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  logical           :: new_name
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  new_name = .false. 
  if(yr.ge.2015) then 
     if(yr.eq.2015) then 
        if(mo.ge.3) then 
           new_name = .true.
        else
           new_name = .false.
        endif
     else
        new_name = .true.
     endif
  endif
     
  if(new_name) then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
          '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '_01D_EQMA_L3SGSNDHG2210210.h5'
  else
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
          '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '_01D_EQMA_L3SGSNDHA1100100.h5'
  endif
end subroutine create_GCOMW_AMSR2L3snd_A_filename


!BOP
! !ROUTINE: create_GCOMW_AMSR2L3snd_D_filename
! \label{create_GCOMW_AMSR2L3snd_D_filename}
! 
! !INTERFACE: 
subroutine create_GCOMW_AMSR2L3snd_D_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GCOMW_AMSR2 filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GCOMW_AMSR2 soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated GCOMW_AMSR2 filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  logical           :: new_name
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  new_name = .false. 
  if(yr.ge.2015) then 
     if(yr.eq.2015) then 
        if(mo.ge.3) then 
           new_name = .true.
        else
           new_name = .false.
        endif
     else
        new_name = .true.
     endif
  endif
     
  if(new_name) then 
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
          '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '_01D_EQMD_L3SGSNDHG2210210.h5'
  else
     filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
          '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
          '_01D_EQMD_L3SGSNDHA1100100.h5'
  endif
end subroutine create_GCOMW_AMSR2L3snd_D_filename
