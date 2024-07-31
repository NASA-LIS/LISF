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
! !ROUTINE: readESACCIsmObs
! \label{readESACCIsmObs}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readESACCIsmObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc
  use LDT_timeMgrMod,   only : LDT_get_julss
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use ESACCIsm_obsMod, only : ESACCIsmobs
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the ESACCI
! soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character(len=LDT_CONST_PATH_LEN)     :: fname
  real              :: smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  ESACCIsmobs(n)%smobs = LDT_rc%udef
  smobs= LDT_rc%udef

  call create_ESACCIsm_filename(ESACCIsmobs(n)%odir, &
       ESACCIsmobs(n)%version, ESACCIsmobs(n)%sensor,&
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
  
  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then
     
     write(LDT_logunit,*) '[INFO] Reading ..',trim(fname)
     call read_ESACCI_data(n, fname, ESACCIsmobs(n)%version, smobs) ! NT: include version in reading
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
              ESACCIsmobs(n)%smobs(c,r) = smobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       ESACCIsmobs(n)%smobs,vlevel=1)

end subroutine readESACCIsmObs


!BOP
! 
! !ROUTINE: read_ESACCI_data
! \label(read_ESACCI_data)
!
! !INTERFACE:
subroutine read_ESACCI_data(n, fname, version, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod,  only : LDT_rc, LDT_domain
  use LDT_logMod,   only : LDT_verify
  use map_utils,    only : latlon_to_ij
  use ESACCIsm_obsMod, only : ESACCIsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  character (len=*)             :: fname
  real                          :: smobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real	                        :: version

! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the ESACCI NETCDF file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  land surface temperature is below freezing, if rain is present, if 
!  RFI is present, if residual error is above 0.5 or if optical depth
!  is above 0.8. Finally the routine combines both the C-band and X-band
!  retrievals. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the ESACCI AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LDT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer           :: sm(ESACCIsmobs(n)%esaccinc,ESACCIsmobs(n)%esaccinr)
  real              :: sm1(ESACCIsmobs(n)%esaccinc,ESACCIsmobs(n)%esaccinr)
  integer           :: flag(ESACCIsmobs(n)%esaccinc,ESACCIsmobs(n)%esaccinr)

  real              :: sm_combined(ESACCIsmobs(n)%esaccinc,ESACCIsmobs(n)%esaccinr)
  real              :: sm_data(ESACCIsmobs(n)%esaccinc*ESACCIsmobs(n)%esaccinr)
  logical*1         :: sm_data_b(ESACCIsmobs(n)%esaccinc*ESACCIsmobs(n)%esaccinr)
  logical*1         :: smobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  integer           :: c,r,i,j
  real              :: ri,rj
  integer           :: nid
  integer           :: smid, flagid
  integer           :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, 'sm',smid)
  call LDT_verify(ios, 'Error nf90_inq_varid: sm')
  
  ios = nf90_inq_varid(nid, 'flag',flagid)
  call LDT_verify(ios, 'Error nf90_inq_varid: flag')
  
  !values
  if(version .lt. 3) then
      ios = nf90_get_var(nid, smid, sm)
  else
      ios = nf90_get_var(nid, smid, sm1)
  endif
  call LDT_verify(ios, 'Error nf90_get_var: sm')
 
  ios = nf90_get_var(nid, flagid,flag)
  call LDT_verify(ios, 'Error nf90_get_var: flag')
  
  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))

  do r=1, ESACCIsmobs(n)%esaccinr
     do c=1, ESACCIsmobs(n)%esaccinc
!------------------------------------------------------------------------
! All data flagged for snow coverage or temperature below zero (flag=1), 
! dense vegetation (flag=2) and no convergence in the ESACCI algorithm 
! (flag =3) and undefined values are masked out. 
!------------------------------------------------------------------------
        if(version .lt. 3) then
            if(flag(c,r).ne.0.or.sm(c,r).lt.0) then 
               sm_combined(c,ESACCIsmobs(n)%esaccinr-r+1) = LDT_rc%udef
            else
              sm_combined(c,ESACCIsmobs(n)%esaccinr-r+1) = sm(c,r)*0.0001
            endif
        else
            if(flag(c,r).ne.0.or.sm1(c,r).lt.0) then ! NT: checking sm and sm1 separately
               sm_combined(c,ESACCIsmobs(n)%esaccinr-r+1) = LDT_rc%udef
            else
               sm_combined(c,ESACCIsmobs(n)%esaccinr-r+1) = sm1(c,r)
            endif
        endif
     enddo
  enddo
 
  do r=1, ESACCIsmobs(n)%esaccinr
     do c=1, ESACCIsmobs(n)%esaccinc
        sm_data(c+(r-1)*ESACCIsmobs(n)%esaccinc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LDT_rc%udef) then 
           sm_data_b(c+(r-1)*ESACCIsmobs(n)%esaccinc) = .true. 
        else
           sm_data_b(c+(r-1)*ESACCIsmobs(n)%esaccinc) = .false.
        endif
        if(sm_combined(c,r).gt.0.5) then 
           sm_combined(c,r) = LDT_rc%udef
           sm_data_b(c+(r-1)*ESACCIsmobs(n)%esaccinc) = .false.
        endif
     enddo
  enddo
  
!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LDT_rc%gridDesc(n,:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       ESACCIsmobs(n)%esaccinc*ESACCIsmobs(n)%esaccinr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       ESACCIsmobs(n)%n11,&
       LDT_rc%udef, ios)

#endif
  
end subroutine read_ESACCI_data

!BOP
! !ROUTINE: create_ESACCIsm_filename
! \label{create_ESACCIsm_filename}
! 
! !INTERFACE: 
subroutine create_ESACCIsm_filename(ndir, version, sensor, yr, mo,da, filename)
! !USES:   
  use LDT_logMod

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  real              :: version
  character(len=*)  :: sensor
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the ESACCI filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the ESACCI soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated ESACCI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  character (len=3) :: cversion3
  character (len=4) :: cversion4
  character (len=8) :: sensortxt  !20220621 Pang  

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  if (version .lt. 10.0) then
    write(unit=cversion3, fmt='(f3.1)') version
  else
    ! For future version, version number > 10, e.g., 10.2
    write(unit=cversion4, fmt='(f4.1)') version
  endif
  
  if(sensor == 'passive') then
      sensortxt = 'PASSIVE'
  elseif(sensor == 'active') then
      sensortxt = 'ACTIVE'
  elseif(sensor == 'combined') then
      sensortxt = 'COMBINED'
  else
      write(LDT_logunit,*) "[ERR] Invalid ESA CCI soil moisture sensor type was chosen."
      write(LDT_logunit,*) "[ERR] Please choose either 'passive', 'active', or 'combined'."
      call LDT_endrun
  endif

  if(version.eq.1) then 
     filename = trim(ndir)//'/'//trim(fyr)//'/ESACCI-L3S_SOILMOISTURE-SSMV-MERGED-' &
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv00.1.nc'
  elseif(version.eq.2) then 
     filename = trim(ndir)//'/'//trim(fyr)//'/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-' &
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv02.0.nc'
  elseif(version.eq.2.2) then 
     filename = trim(ndir)//'/'//trim(fyr)//'/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-' &
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv02.2.nc'
  else
     ! NT: for versions after 2.2
     if (version .lt. 10.0) then
         filename = trim(ndir)//'/'//trim(fyr)//'/ESACCI-SOILMOISTURE-L3S-SSMV-'//trim(sensortxt)//'-' & !20220621 Pang
              //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv0'//cversion3//'.nc'
     else
         filename = trim(ndir)//'/'//trim(fyr)//'/ESACCI-SOILMOISTURE-L3S-SSMV-'//trim(sensortxt)//'-' & !20220621 Pang
              //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv'//cversion4//'.nc'
     endif
  endif
end subroutine create_ESACCIsm_filename
