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
! !ROUTINE: readESACCIsmObs
! \label{readESACCIsmObs}
!
! !INTERFACE: 
subroutine readESACCIsmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, & 
       LVT_releaseUnitNumber
  use LVT_timeMgrMod,   only : LVT_get_julss
  use ESACCIsm_obsMod, only : ESACCIsmobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine provides the data reader for the standard 
! LPRM soil moisture retrieval product. 
! 
! !NOTES: 
!  The mismatches between the LVT time and LPRM data time is not
!  handled. This is not an issue as long as LVT analysis is not
!  done for sub-daily intervals. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
!  25 May 2012: Sujay Kumar, Updated for LPRM version 5. 
! 
!EOP
  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character*100     :: fname
  real              :: smobs(LVT_rc%lnc*LVT_rc%lnr)
  real              :: lat,lon
  integer	    :: version

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  ESACCIsmobs(source)%smobs = LVT_rc%udef
  smobs= LVT_rc%udef
        
  if(ESACCIsmobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     ESACCIsmobs(source)%startmode = .false. 
     version = ESACCIsmobs(source)%version
     call create_ESACCIsm_filename(ESACCIsmobs(source)%odir, &
          ESACCIsmobs(source)%version, & 
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), fname)

     inquire(file=trim(fname),exist=file_exists)

     if(file_exists) then
        write(LVT_logunit,*) '[INFO] Reading ESA CCI file ... ',trim(fname)
        call read_ESACCI_data(source, fname,version,smobs)
     else
        write(LVT_logunit,*) '[WARN] Unable to read-in ESA CCI file ...',trim(fname)
        write(LVT_logunit,*) '[WARN]  Note: Check length of filepath, as a long '
        write(LVT_logunit,*) '[WARN]     pathname can lead to not reading the file.' 
     endif

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(smobs(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
              ESACCIsmobs(source)%smobs(c,r) = smobs(c+(r-1)*LVT_rc%lnc)
           endif
        enddo
     enddo
  endif
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       ESACCIsmobs(source)%smobs,vlevel=1,units="m3/m3")

end subroutine readESACCIsmObs


!BOP
! 
! !ROUTINE: read_ESACCI_data
! \label(read_ESACCI_data)
!
! !INTERFACE:
subroutine read_ESACCI_data(source, fname,version, smobs_ip)
! 
! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LVT_coreMod,  only : LVT_rc
  use LVT_logMod,   only : LVT_verify
  use map_utils,    only : latlon_to_ij
  use ESACCIsm_obsMod, only : ESACCIsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer              :: source
  character (len=*)    :: fname
  real                 :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)
  integer	       :: version
!
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
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP


  integer        :: sm(ESACCIsmobs(source)%esaccinc,ESACCIsmobs(source)%esaccinr)
  real           :: sm1(ESACCIsmobs(source)%esaccinc,ESACCIsmobs(source)%esaccinr)

  integer        :: flag(ESACCIsmobs(source)%esaccinc,ESACCIsmobs(source)%esaccinr)
  real           :: sm_combined(ESACCIsmobs(source)%esaccinc,ESACCIsmobs(source)%esaccinr)
  real           :: sm_data(ESACCIsmobs(source)%esaccinc*ESACCIsmobs(source)%esaccinr)
  logical*1      :: sm_data_b(ESACCIsmobs(source)%esaccinc*ESACCIsmobs(source)%esaccinr)
  logical*1      :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  integer        :: c,r,i,j
  real           :: rlat,rlon,ri,rj
  integer        :: nid
  integer        :: smid, flagid
  integer        :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LVT_verify(ios,'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(nid, 'sm',smid)
  call LVT_verify(ios, 'Error nf90_inq_varid: sm')
  
  ios = nf90_inq_varid(nid, 'flag',flagid)
  call LVT_verify(ios, 'Error nf90_inq_varid: flag')
  
  !values
  if(version .eq.1 .or. version .eq. 2) then
  ios = nf90_get_var(nid, smid, sm)
  call LVT_verify(ios, 'Error nf90_get_var: sm')
  endif

  if(version .eq.3) then
  ios = nf90_get_var(nid, smid, sm1)
  call LVT_verify(ios, 'Error nf90_get_var: sm')
  endif

  
  ios = nf90_get_var(nid, flagid,flag)
  call LVT_verify(ios, 'Error nf90_get_var: flag')
  
  ios = nf90_close(ncid=nid)
  call LVT_verify(ios,'Error closing file '//trim(fname))

  do r=1, ESACCIsmobs(source)%esaccinr
     do c=1, ESACCIsmobs(source)%esaccinc
!------------------------------------------------------------------------
! All data flagged for snow coverage or temperature below zero (flag=1), 
! dense vegetation (flag=2) and no convergence in the ESACCI algorithm 
! (flag =3) and undefined values are masked out. 
!------------------------------------------------------------------------
        if(flag(c,r).ne.0.or.sm(c,r).lt.0) then 
           sm_combined(c,ESACCIsmobs(source)%esaccinr-r+1) = LVT_rc%udef
        else
	 if(version .eq. 1 .or. version .eq. 2)then
           sm_combined(c,ESACCIsmobs(source)%esaccinr-r+1) = sm(c,r)*0.0001
	 endif
	 if(version .eq. 3)then
	   sm_combined(c,ESACCIsmobs(source)%esaccinr-r+1) = sm1(c,r)
	 endif
        endif
     enddo
  enddo
 
  do r=1, ESACCIsmobs(source)%esaccinr
     do c=1, ESACCIsmobs(source)%esaccinc
        sm_data(c+(r-1)*ESACCIsmobs(source)%esaccinc) = sm_combined(c,r)
        if(sm_combined(c,r).ne.LVT_rc%udef) then 
           sm_data_b(c+(r-1)*ESACCIsmobs(source)%esaccinc) = .true. 
        else
           sm_data_b(c+(r-1)*ESACCIsmobs(source)%esaccinc) = .false.
        endif
        if(sm_combined(c,r).gt.0.5) then 
           sm_combined(c,r) = LVT_rc%udef
           sm_data_b(c+(r-1)*ESACCIsmobs(source)%esaccinc) = .false.
        endif
     enddo
  enddo
  
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
  call bilinear_interp(LVT_rc%gridDesc(:),&
       sm_data_b, sm_data, smobs_b_ip, smobs_ip, &
       ESACCIsmobs(source)%esaccinc*ESACCIsmobs(source)%esaccinr, &
       LVT_rc%lnc*LVT_rc%lnr, &
       ESACCIsmobs(source)%rlat, ESACCIsmobs(source)%rlon, &
       ESACCIsmobs(source)%w11, ESACCIsmobs(source)%w12, &
       ESACCIsmobs(source)%w21, ESACCIsmobs(source)%w22, &
       ESACCIsmobs(source)%n11, ESACCIsmobs(source)%n12, &
       ESACCIsmobs(source)%n21, ESACCIsmobs(source)%n22, &
       LVT_rc%udef, ios)

#endif
  
end subroutine read_ESACCI_data

!BOP
! !ROUTINE: create_ESACCIsm_filename
! \label{create_ESACCIsm_filename}
! 
! !INTERFACE: 
subroutine create_ESACCIsm_filename(ndir, version, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: version
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the ESACCI filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the ESACCI soil moisture directory
!  \item[version] version of the ESACCI data (1 or 2)
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated ESACCI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  if(version.eq.1) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/ESACCI-L3S_SOILMOISTURE-SSMV-MERGED-' &
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv00.1.nc'       
  elseif(version.eq.2) then 
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-' & 
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv02.2.nc'           
  elseif(version.eq.3) then
     filename = trim(ndir)//'/'//trim(fyr)//&
          '/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-' &
          //trim(fyr)//trim(fmo)//trim(fda)//'000000-fv03.2.nc'
  endif

end subroutine create_ESACCIsm_filename
