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
! !ROUTINE: readLPRMvodObs
! \label{readLPRMvodObs}
! 
! !REVISION HISTORY: 
!  28 May 2019: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readLPRMvodObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use LPRMvod_obsMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the 
! LPRM vegetation optical depth retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,c1,r1
  integer           :: ios
  integer           :: ftn
  character(len=LDT_CONST_PATH_LEN)     :: fname
  real              :: vodobs(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real              :: vod_file(LPRMvodobs(n)%nc,LPRMvodobs(n)%nr)
  real              :: vod_inp(LPRMvodobs(n)%nc*LPRMvodobs(n)%nr)
  logical*1         :: vod_b_inp(LPRMvodobs(n)%nc*LPRMvodobs(n)%nr)
  real              :: vod_out(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1         :: vod_b_out(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: lat(LPRMvodobs(n)%nr)
  real              :: lon(LPRMvodobs(n)%nc)
  integer           :: vodid,flagid,latid,lonid


!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  vodobs= LDT_rc%udef

  call create_LPRMvod_filename(LPRMvodobs(n)%odir, &
       LPRMvodobs(n)%data_designation,&
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)
  
  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then
     
     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=ftn)
     call LDT_verify(ios,'Error opening file '//trim(fname))
     
     ios = nf90_inq_varid(ftn, 'lat',latid)
     call LDT_verify(ios, 'Error nf90_inq_varid: lat')
     
     ios = nf90_inq_varid(ftn, 'lon',lonid)
     call LDT_verify(ios, 'Error nf90_inq_varid: lon')
     
     ios = nf90_inq_varid(ftn, 'vod',vodid)
     call LDT_verify(ios, 'Error nf90_inq_varid: vod')
     
     ios = nf90_get_var(ftn, latid, lat)
     call LDT_verify(ios, 'Error nf90_get_var: lat')
     
     ios = nf90_get_var(ftn, lonid, lon)
     call LDT_verify(ios, 'Error nf90_get_var: lon')
     
     ios = nf90_get_var(ftn, vodid, vod_file)
     call LDT_verify(ios, 'Error nf90_get_var: vod')
     
     ios = nf90_close(ncid=ftn)
     call LDT_verify(ios,'Error closing file '//trim(fname))
#endif
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     vod_inp = LDT_rc%udef
     vod_b_inp  = .false. 


     do r=1,LPRMvodobs(n)%nr
        do c=1,LPRMvodobs(n)%nc
           r1 = nint((lat(r)+89.875)/0.25)+1
           c1 = nint((lon(c)+179.875)/0.25)+1
           
           if(vod_file(c,r).ne.-999999.0) then 
              vod_inp(c1+(r1-1)*LPRMvodobs(n)%nc) = & 
                      vod_file(c,r)
              vod_b_inp(c1+(r1-1)*LPRMvodobs(n)%nc) = & 
                   .true. 
           endif
        enddo
     enddo
     
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          vod_b_inp,&
          vod_inp,&
          vod_b_out,&
          vod_out,&
          LPRMvodobs(n)%nc*LPRMvodobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          LPRMvodobs(n)%w11,&
          LPRMvodobs(n)%w12,&
          LPRMvodobs(n)%w21,&
          LPRMvodobs(n)%w22,&
          LPRMvodobs(n)%n11,&
          LPRMvodobs(n)%n12,&
          LPRMvodobs(n)%n21,&
          LPRMvodobs(n)%n22,&
          LDT_rc%udef,ios)


     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           vodobs(c,r) = vod_out(c+(r-1)*LDT_rc%lnc(n))
        enddo
     enddo

  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%vod_obs,&
       vodobs,vlevel=1)

end subroutine readLPRMvodObs


!BOP
! !ROUTINE: create_LPRMvod_filename
! \label{create_LPRMvod_filename}
! 
! !INTERFACE: 
subroutine create_LPRMvod_filename(odir, designation,yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: designation
  integer           :: yr, mo, da
  character (len=*) :: odir
! 
! !DESCRIPTION: 
!  This subroutine creates the timestamped LPRM VOD filename 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the LPRM VOD directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated LPRM VOD filename
! \end{description}
!EOP

  character*4      :: yyyy
  character*2      :: mm,dd

  write(yyyy,'(i4.4)') yr
  write(mm,'(i2.2)') mo
  write(dd,'(i2.2)') da
  
  filename=trim(odir)//'/'//trim(designation)//'/'//trim(yyyy)//&
       '/vodca_v01-0_'//trim(designation)//'_'//trim(yyyy)//'-'//&
       trim(mm)//'-'//trim(dd)//'.nc'
         
end subroutine create_LPRMvod_filename
