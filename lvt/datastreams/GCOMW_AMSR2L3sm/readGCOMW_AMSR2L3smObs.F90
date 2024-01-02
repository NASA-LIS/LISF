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
! !ROUTINE: readGCOMW_AMSR2L3smObs
! \label{readGCOMW_AMSR2L3smObs}
!
! !INTERFACE: 
subroutine readGCOMW_AMSR2L3smObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use GCOMW_AMSR2L3sm_obsMod, only : GCOMW_AMSR2L3smobs

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
  real              :: smobs(LVT_rc%lnc, LVT_rc%lnr)
  real              :: smobs_t(LVT_rc%lnc*LVT_rc%lnr)
  real              :: smobs_av(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: count_av(LVT_rc%lnc, LVT_rc%lnr)
  type(ESMF_Time)   :: cTime, stTime, enTime

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LVT_rc%dhr(source))*3600 + &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  smobs= LVT_rc%udef
        
  if(GCOMW_AMSR2L3smobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     GCOMW_AMSR2L3smobs(source)%startmode = .false. 

     if(LVT_rc%smoothObs.eq.1) then 
        smobs_av = 0.0
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

           call create_GCOMW_AMSR2L3sm_A_filename(GCOMW_AMSR2L3smobs(source)%odir, &
                yr, mo, da, fname_A)
           
           call create_GCOMW_AMSR2L3sm_D_filename(GCOMW_AMSR2L3smobs(source)%odir, &
                yr, mo, da, fname_D)

           smobs_t = LVT_rc%udef
           call read_GCOMW_AMSR2_data(source, fname_A, fname_D,smobs_t)
           
           cTime = cTime + LVT_obsSmTwI

           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(smobs_t(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                    smobs_av(c,r) = smobs_av(c,r) + smobs_t(c+(r-1)*LVT_rc%lnc)
                    count_av(c,r) = count_av(c,r) + 1
                 endif
              enddo
           enddo
        enddo

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(count_av(c,r).gt.0) then 
                 smobs(c,r) = smobs_av(c,r)/count_av(c,r)
              else
                 smobs(c,r) = LVT_rc%udef
              endif
           enddo
        enddo

     else

        call create_GCOMW_AMSR2L3sm_A_filename(GCOMW_AMSR2L3smobs(source)%odir, &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), fname_A)
        
        call create_GCOMW_AMSR2L3sm_D_filename(GCOMW_AMSR2L3smobs(source)%odir, &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), fname_D)
        
        call read_GCOMW_AMSR2_data(source, fname_A, fname_D,smobs_t)

        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              smobs(c,r) = smobs_t(c+(r-1)*LVT_rc%lnc)
           enddo
        enddo
     endif

  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source,&
       smobs,vlevel=1,units="m3/m3")
 
end subroutine readGCOMW_AMSR2L3smObs


!BOP
! 
! !ROUTINE: read_GCOMW_AMSR2_data
! \label(read_GCOMW_AMSR2_data)
!
! !INTERFACE:
subroutine read_GCOMW_AMSR2_data(source, fname_A, fname_D, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LVT_coreMod
  use LVT_logMod
  use map_utils,    only : latlon_to_ij
  use GCOMW_AMSR2L3sm_obsMod, only : GCOMW_AMSR2L3smobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source 
  character (len=*)             :: fname_A
  character (len=*)             :: fname_D
  real                          :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname_A]        name of the GCOMW_AMSR2 AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer             :: sm(1,GCOMW_AMSR2L3smobs(source)%amsr2nc,GCOMW_AMSR2L3smobs(source)%amsr2nr)
  integer           :: time(GCOMW_AMSR2L3smobs(source)%amsr2nc,GCOMW_AMSR2L3smobs(source)%amsr2nr)
  integer           :: smtime(GCOMW_AMSR2L3smobs(source)%amsr2nc,GCOMW_AMSR2L3smobs(source)%amsr2nr)
  real                        :: sm_combined(GCOMW_AMSR2L3smobs(source)%amsr2nc,GCOMW_AMSR2L3smobs(source)%amsr2nr)
  real                        :: sm_data(GCOMW_AMSR2L3smobs(source)%amsr2nc*GCOMW_AMSR2L3smobs(source)%amsr2nr)
  logical*1                   :: sm_data_b(GCOMW_AMSR2L3smobs(source)%amsr2nc*GCOMW_AMSR2L3smobs(source)%amsr2nr)
  logical*1                   :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  logical                     :: file_exists
  integer                     :: c1,c,r,i,j
  real                        :: rlat,rlon,ri,rj
  integer                     :: nid
  integer                     :: smid, timeid
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  smtime = -1.0
  sm_combined = LVT_rc%udef

  inquire(file=fname_A, exist=file_exists) 
  if(file_exists) then 
     write (LVT_logunit, *) "[INFO] Reading "//trim(fname_A)
     ios = nf90_open(path=trim(fname_A),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error opening file '//trim(fname_A))
     
     ios = nf90_inq_varid(nid, "Geophysical Data",smid)
     call LVT_verify(ios, 'Error nf90_inq_varid: sm')
     
     ios = nf90_inq_varid(nid, "Time Information",timeid)
     call LVT_verify(ios, 'Error nf90_inq_varid: time')
  
     !values
     ios = nf90_get_var(nid, smid, sm)
     call LVT_verify(ios, 'Error nf90_get_var: sm')
     
     ios = nf90_get_var(nid, timeid,time)
     call LVT_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LVT_verify(ios,'Error closing file '//trim(fname_A))


     do r=1, GCOMW_AMSR2L3smobs(source)%amsr2nr
        do c=1, GCOMW_AMSR2L3smobs(source)%amsr2nc
           c1 = c+GCOMW_AMSR2L3smobs(source)%amsr2nc/2-1
           if(c1.gt.GCOMW_AMSR2L3smobs(source)%amsr2nc) then 
              c1 = c-GCOMW_AMSR2L3smobs(source)%amsr2nc/2-1
           endif
           
           if(time(c,r).lt.0.or.sm(1,c,r).le.0.001) then
              sm_combined(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1) = LVT_rc%udef
           else
              sm_combined(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1) = &
                   (sm(1,c,r)*0.1)/100.0
              smtime(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1) = time(c,r)
           endif
        enddo
     enddo

  endif

  inquire(file=fname_D, exist=file_exists) 
  if(file_exists) then 
     write (LVT_logunit, *) "[INFO] Reading "//trim(fname_D)
     ios = nf90_open(path=trim(fname_D),mode=NF90_NOWRITE,ncid=nid)
     call LVT_verify(ios,'Error opening file '//trim(fname_d))
     
     ios = nf90_inq_varid(nid, "Geophysical Data",smid)
     call LVT_verify(ios, 'Error nf90_inq_varid: sm')
     
     ios = nf90_inq_varid(nid, "Time Information",timeid)
     call LVT_verify(ios, 'Error nf90_inq_varid: time')
  
     !values
     ios = nf90_get_var(nid, smid, sm)
     call LVT_verify(ios, 'Error nf90_get_var: sm')
     
     ios = nf90_get_var(nid, timeid,time)
     call LVT_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LVT_verify(ios,'Error closing file '//trim(fname_d))


     do r=1, GCOMW_AMSR2L3smobs(source)%amsr2nr
        do c=1, GCOMW_AMSR2L3smobs(source)%amsr2nc
           c1 = c+GCOMW_AMSR2L3smobs(source)%amsr2nc/2-1
           if(c1.gt.GCOMW_AMSR2L3smobs(source)%amsr2nc) then 
              c1 = c-GCOMW_AMSR2L3smobs(source)%amsr2nc/2-1
           endif
           if(sm_combined(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1).gt.0) then 
              !if valid data already exists, then overwrite only 
              !if a newer data is available. 
              if((time(c,r).gt.0).and.&
                   (time(c,r).gt.smtime(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1)).and.&
                   (sm(1,c,r).gt.0.001)) then 
                 sm_combined(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1) = &
                      (sm(1,c,r)*0.1)/100.0
              end if
           else
              if(time(c,r).lt.0.or.sm(1,c,r).le.0.001) then
                 sm_combined(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1) = LVT_rc%udef
              else
                 sm_combined(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1) = &
                      (sm(1,c,r)*0.1)/100.0
                 smtime(c1,GCOMW_AMSR2L3smobs(source)%amsr2nr-r+1) = time(c,r)
              endif
              
           endif
        enddo
     enddo

  endif

  do r=1, GCOMW_AMSR2L3smobs(source)%amsr2nr
     do c=1, GCOMW_AMSR2L3smobs(source)%amsr2nc
        sm_data(c+(r-1)*GCOMW_AMSR2L3smobs(source)%amsr2nc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LVT_rc%udef) then 
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3smobs(source)%amsr2nc) = .true. 
        else
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3smobs(source)%amsr2nc) = .false.
        endif
        if(sm_combined(c,r).gt.1) then 
           sm_combined(c,r) = LVT_rc%udef
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3smobs(source)%amsr2nc) = .false.
        endif
     enddo
  enddo

  if(LVT_isatAfinerResolution(GCOMW_AMSR2L3smobs(source)%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
          GCOMW_AMSR2L3smobs(source)%amsr2nc*GCOMW_AMSR2L3smobs(source)%amsr2nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          GCOMW_AMSR2L3smobs(source)%rlat, GCOMW_AMSR2L3smobs(source)%rlon, &
          GCOMW_AMSR2L3smobs(source)%w11, GCOMW_AMSR2L3smobs(source)%w12, &
          GCOMW_AMSR2L3smobs(source)%w21, GCOMW_AMSR2L3smobs(source)%w22, &
          GCOMW_AMSR2L3smobs(source)%n11, GCOMW_AMSR2L3smobs(source)%n12, &
          GCOMW_AMSR2L3smobs(source)%n21, GCOMW_AMSR2L3smobs(source)%n22, &
          LVT_rc%udef, ios)
  else
     call upscaleByAveraging(&
          GCOMW_AMSR2L3smobs(source)%amsr2nc*GCOMW_AMSR2L3smobs(source)%amsr2nr,&
          LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef, &
          GCOMW_AMSR2L3smobs(source)%n11,sm_data_b, sm_data, smobs_b_ip,smobs_ip)
     
  endif

!  open(100,file='smobs.bin',form='unformatted')
!  write(100) smobs_ip
!  close(100)
!  stop
#endif
  
end subroutine read_GCOMW_AMSR2_data

!BOP
! !ROUTINE: create_GCOMW_AMSR2L3sm_A_filename
! \label{create_GCOMW_AMSR2L3sm_A_filename}
! 
! !INTERFACE: 
subroutine create_GCOMW_AMSR2L3sm_A_filename(ndir, yr, mo,da, filename)
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
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
       '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
       '_01D_EQMA_L3SGSMCHF2210210.h5'        
!  filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
!       '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
!       '_01D_EQMA_L3SGSMCHA1110100.h5'
 

end subroutine create_GCOMW_AMSR2L3sm_A_filename


!BOP
! !ROUTINE: create_GCOMW_AMSR2L3sm_D_filename
! \label{create_GCOMW_AMSR2L3sm_D_filename}
! 
! !INTERFACE: 
subroutine create_GCOMW_AMSR2L3sm_D_filename(ndir, yr, mo,da, filename)
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
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
!  filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
!       '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
!       '_01D_EQMD_L3SGSMCHA1110100.h5'         
  filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
       '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
       '_01D_EQMD_L3SGSMCHF2210210.h5'
end subroutine create_GCOMW_AMSR2L3sm_D_filename
