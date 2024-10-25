!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readGCOMW_AMSR2L3smObs
! \label{readGCOMW_AMSR2L3smObs}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readGCOMW_AMSR2L3smObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc
  use LDT_timeMgrMod,   only : LDT_get_julss
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod  
  use GCOMW_AMSR2L3sm_obsMod, only : GCOMW_AMSR2L3smobs
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the GCOMW_AMSR2
! soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character(len=LDT_CONST_PATH_LEN)     :: fname_A, fname_D
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  GCOMW_AMSR2L3smobs(n)%smobs = LDT_rc%udef
  smobs= LDT_rc%udef

  call create_GCOMW_AMSR2L3sm_A_filename(GCOMW_AMSR2L3smobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname_A)
  
  call create_GCOMW_AMSR2L3sm_D_filename(GCOMW_AMSR2L3smobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname_D)
  
  call read_GCOMW_AMSR2_data(n, fname_A, fname_D,smobs)
  
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(smobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
           GCOMW_AMSR2L3smobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
        endif
     enddo
  enddo

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       GCOMW_AMSR2L3smobs(n)%smobs,vlevel=1)

end subroutine readGCOMW_AMSR2L3smObs


!BOP
! 
! !ROUTINE: read_GCOMW_AMSR2_data
! \label(read_GCOMW_AMSR2_data)
!
! !INTERFACE:
subroutine read_GCOMW_AMSR2_data(n, fname_A, fname_D, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod
  use LDT_logMod
  use map_utils,    only : latlon_to_ij
  use GCOMW_AMSR2L3sm_obsMod, only : GCOMW_AMSR2L3smobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname_A
  character (len=*)             :: fname_D
  real                          :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname_A]        name of the GCOMW_AMSR2 AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LDT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer             :: sm(1,GCOMW_AMSR2L3smobs(n)%amsr2nc,GCOMW_AMSR2L3smobs(n)%amsr2nr)
  integer           :: time(GCOMW_AMSR2L3smobs(n)%amsr2nc,GCOMW_AMSR2L3smobs(n)%amsr2nr)
  integer           :: smtime(GCOMW_AMSR2L3smobs(n)%amsr2nc,GCOMW_AMSR2L3smobs(n)%amsr2nr)
  real                        :: sm_combined(GCOMW_AMSR2L3smobs(n)%amsr2nc,GCOMW_AMSR2L3smobs(n)%amsr2nr)
  real                        :: sm_data(GCOMW_AMSR2L3smobs(n)%amsr2nc*GCOMW_AMSR2L3smobs(n)%amsr2nr)
  logical*1                   :: sm_data_b(GCOMW_AMSR2L3smobs(n)%amsr2nc*GCOMW_AMSR2L3smobs(n)%amsr2nr)
  logical*1                   :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  logical                     :: file_exists
  integer                     :: c1,c,r,i,j
  real                        :: ri,rj
  integer                     :: nid
  integer                     :: smid, timeid
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  smtime = -1.0
  sm_combined = LDT_rc%udef

  inquire(file=fname_A, exist=file_exists) 
  if(file_exists) then 
     write (LDT_logunit, *) "Reading "//trim(fname_A)
     ios = nf90_open(path=trim(fname_A),mode=NF90_NOWRITE,ncid=nid)
     call LDT_verify(ios,'Error opening file '//trim(fname_A))
     
     ios = nf90_inq_varid(nid, "Geophysical Data",smid)
     call LDT_verify(ios, 'Error nf90_inq_varid: sm')
     
     ios = nf90_inq_varid(nid, "Time Information",timeid)
     call LDT_verify(ios, 'Error nf90_inq_varid: time')
  
     !values
     ios = nf90_get_var(nid, smid, sm)
     call LDT_verify(ios, 'Error nf90_get_var: sm')
     
     ios = nf90_get_var(nid, timeid,time)
     call LDT_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LDT_verify(ios,'Error closing file '//trim(fname_A))


     do r=1, GCOMW_AMSR2L3smobs(n)%amsr2nr
        do c=1, GCOMW_AMSR2L3smobs(n)%amsr2nc
           c1 = c+GCOMW_AMSR2L3smobs(n)%amsr2nc/2-1
           if(c1.gt.GCOMW_AMSR2L3smobs(n)%amsr2nc) then 
              c1 = c-GCOMW_AMSR2L3smobs(n)%amsr2nc/2-1
           endif
           
           if(time(c,r).lt.0.or.sm(1,c,r).le.0.001) then
              sm_combined(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1) = LDT_rc%udef
           else
              sm_combined(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1) = &
                   (sm(1,c,r)*0.1)/100.0
              smtime(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1) = time(c,r)
           endif
        enddo
     enddo

  endif

  inquire(file=fname_D, exist=file_exists) 
  if(file_exists) then 
     write (LDT_logunit, *) "Reading "//trim(fname_D)
     ios = nf90_open(path=trim(fname_D),mode=NF90_NOWRITE,ncid=nid)
     call LDT_verify(ios,'Error opening file '//trim(fname_d))
     
     ios = nf90_inq_varid(nid, "Geophysical Data",smid)
     call LDT_verify(ios, 'Error nf90_inq_varid: sm')
     
     ios = nf90_inq_varid(nid, "Time Information",timeid)
     call LDT_verify(ios, 'Error nf90_inq_varid: time')
  
     !values
     ios = nf90_get_var(nid, smid, sm)
     call LDT_verify(ios, 'Error nf90_get_var: sm')
     
     ios = nf90_get_var(nid, timeid,time)
     call LDT_verify(ios, 'Error nf90_get_var: time')
     
     ios = nf90_close(ncid=nid)
     call LDT_verify(ios,'Error closing file '//trim(fname_d))


     do r=1, GCOMW_AMSR2L3smobs(n)%amsr2nr
        do c=1, GCOMW_AMSR2L3smobs(n)%amsr2nc
           c1 = c+GCOMW_AMSR2L3smobs(n)%amsr2nc/2-1
           if(c1.gt.GCOMW_AMSR2L3smobs(n)%amsr2nc) then 
              c1 = c-GCOMW_AMSR2L3smobs(n)%amsr2nc/2-1
           endif
           if(sm_combined(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1).gt.0) then 
              !if valid data already exists, then overwrite only 
              !if a newer data is available. 
              if((time(c,r).gt.0).and.&
                   (time(c,r).gt.smtime(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1)).and.&
                   (sm(1,c,r).gt.0.001)) then 
                 sm_combined(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1) = &
                      (sm(1,c,r)*0.1)/100.0
              end if
           else
              if(time(c,r).lt.0.or.sm(1,c,r).le.0.001) then
                 sm_combined(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1) = LDT_rc%udef
              else
                 sm_combined(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1) = &
                      (sm(1,c,r)*0.1)/100.0
                 smtime(c1,GCOMW_AMSR2L3smobs(n)%amsr2nr-r+1) = time(c,r)
              endif
              
           endif
        enddo
     enddo

  endif

  do r=1, GCOMW_AMSR2L3smobs(n)%amsr2nr
     do c=1, GCOMW_AMSR2L3smobs(n)%amsr2nc
        sm_data(c+(r-1)*GCOMW_AMSR2L3smobs(n)%amsr2nc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LDT_rc%udef) then 
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3smobs(n)%amsr2nc) = .true. 
        else
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3smobs(n)%amsr2nc) = .false.
        endif
        if(sm_combined(c,r).gt.1) then 
           sm_combined(c,r) = LDT_rc%udef
           sm_data_b(c+(r-1)*GCOMW_AMSR2L3smobs(n)%amsr2nc) = .false.
        endif
     enddo
  enddo

  if(LDT_isLDTatAfinerResolution(n,GCOMW_AMSR2L3smobs(n)%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
          GCOMW_AMSR2L3smobs(n)%amsr2nc*GCOMW_AMSR2L3smobs(n)%amsr2nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          GCOMW_AMSR2L3smobs(n)%w11, GCOMW_AMSR2L3smobs(n)%w12, &
          GCOMW_AMSR2L3smobs(n)%w21, GCOMW_AMSR2L3smobs(n)%w22, &
          GCOMW_AMSR2L3smobs(n)%n11, GCOMW_AMSR2L3smobs(n)%n12, &
          GCOMW_AMSR2L3smobs(n)%n21, GCOMW_AMSR2L3smobs(n)%n22, &
          LDT_rc%udef, ios)
  else
     call upscaleByAveraging(&
          GCOMW_AMSR2L3smobs(n)%amsr2nc*GCOMW_AMSR2L3smobs(n)%amsr2nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
          GCOMW_AMSR2L3smobs(n)%n11,sm_data_b, sm_data, smobs_b_ip,smobs_ip)
     
  endif

!  open(100,file='test_ip.bin',form='unformatted')
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
       '_01D_EQMA_L3SGSMCHA1110100.h5'         
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
 
  filename = trim(ndir)//'/'//trim(fyr)//'.'//trim(fmo)//&
       '/GW1AM2_'//trim(fyr)//trim(fmo)//trim(fda)//&
       '_01D_EQMD_L3SGSMCHA1110100.h5'         
end subroutine create_GCOMW_AMSR2L3sm_D_filename
