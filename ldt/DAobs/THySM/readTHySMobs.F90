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
! !ROUTINE: readTHySMobs
! \label{readTHySMobs}
! 
! !REVISION HISTORY: 
!  29 Mar 2021: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readTHySMobs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use THySM_obsMod
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
! the THySM soil moisture data
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists1,file_exists2
  integer           :: c,r,c1,r1
  integer           :: ios
  integer           :: ftn
  character(len=LDT_CONST_PATH_LEN)     :: fname_AM, fname_PM
  real              :: smobs(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real              :: sm_file(THySMobs(n)%nc, THySMobs(n)%nr)
  real              :: sm_inp(THySMobs(n)%nc*THySMobs(n)%nr)
  logical*1         :: sm_b_inp(THySMobs(n)%nc*THySMobs(n)%nr)
  real              :: sm_out(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1         :: sm_b_out(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: lat(THySMobs(n)%nr)
  real              :: lon(THySMobs(n)%nc)
  integer           :: smid,flagid,latid,lonid


!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  smobs= LDT_rc%udef

  !AM is descending and PM is ascending
  call create_THySM_filename(THySMobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, 'AM',fname_AM)
  
  inquire(file=trim(fname_AM),exist=file_exists1)
  if(file_exists1) then

     write(LDT_logunit,*) '[INFO] Reading ',trim(fname_AM)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     ios = nf90_open(path=trim(fname_AM),mode=NF90_NOWRITE,ncid=ftn)
     call LDT_verify(ios,'Error opening file '//trim(fname_AM))
     
     ios = nf90_inq_varid(ftn, 'Band1',smid)
     call LDT_verify(ios, 'Error nf90_inq_varid: Band1')
     
     ios = nf90_get_var(ftn, smid, sm_file)
     call LDT_verify(ios, 'Error nf90_get_var: Band1')
     
     ios = nf90_close(ncid=ftn)
     call LDT_verify(ios,'Error closing file '//trim(fname_AM))
#endif
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname_AM)

     sm_b_inp  = .false. 
     sm_inp = LDT_rc%udef
     
     do r=1,THySMobs(n)%nr
        do c=1,THySMobs(n)%nc
           if(sm_file(c,r).gt.0) then 
              sm_b_inp(c+(r-1)*THySMobs(n)%nc) = & 
                   .true.
              sm_inp(c+(r-1)*THySMobs(n)%nc) = &
                   sm_file(c,r)
           endif
        enddo
     enddo

  endif

  call create_THySM_filename(THySMobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, 'PM',fname_PM)
  
  inquire(file=trim(fname_PM),exist=file_exists2)
  if(file_exists2) then

     write(LDT_logunit,*) '[INFO] Reading ',trim(fname_PM)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     ios = nf90_open(path=trim(fname_PM),mode=NF90_NOWRITE,ncid=ftn)
     call LDT_verify(ios,'Error opening file '//trim(fname_PM))
     
     ios = nf90_inq_varid(ftn, 'Band1',smid)
     call LDT_verify(ios, 'Error nf90_inq_varid: Band1')
     
     ios = nf90_get_var(ftn, smid, sm_file)
     call LDT_verify(ios, 'Error nf90_get_var: Band1')
     
     ios = nf90_close(ncid=ftn)
     call LDT_verify(ios,'Error closing file '//trim(fname_PM))
#endif
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname_PM)

     do r=1,THySMobs(n)%nr
        do c=1,THySMobs(n)%nc
           if(sm_file(c,r).gt.0.and.&
                sm_inp(c+(r-1)*THySMobs(n)%nc).eq.-9999.0) then 
              sm_b_inp(c+(r-1)*THySMobs(n)%nc) = & 
                   .true.
              sm_inp(c+(r-1)*THySMobs(n)%nc) = &
                   sm_file(c,r)
           endif
        enddo
     enddo

  endif

  if(file_exists1.or.file_exists2) then 
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          sm_b_inp,&
          sm_inp,&
          sm_b_out,&
          sm_out,&
          THySMobs(n)%nc*THySMobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          THySMobs(n)%w11,&
          THySMobs(n)%w12,&
          THySMobs(n)%w21,&
          THySMobs(n)%w22,&
          THySMobs(n)%n11,&
          THySMobs(n)%n12,&
          THySMobs(n)%n21,&
          THySMobs(n)%n22,&
          LDT_rc%udef,ios)


     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           smobs(c,r) = sm_out(c+(r-1)*LDT_rc%lnc(n))
        enddo
     enddo

  endif
  
  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       smobs,vlevel=1)

end subroutine readTHySMobs


!BOP
! !ROUTINE: create_THySM_filename
! \label{create_THySM_filename}
! 
! !INTERFACE: 
subroutine create_THySM_filename(odir, yr, mo,da, overpass, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: odir
  character (len=*) :: overpass
! 
! !DESCRIPTION: 
!  This subroutine creates the timestamped THySM filename 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the THySM directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated THySM filename
! \end{description}
!EOP

  character*4      :: yyyy
  character*2      :: mm,dd

  write(yyyy,'(i4.4)') yr
  write(mm,'(i2.2)') mo
  write(dd,'(i2.2)') da
  
  filename=trim(odir)//'/'//trim(yyyy)//&
       '/SMAP-HYB-1KM-DAILY_'//trim(yyyy)//'.'//&
       trim(mm)//'.'//trim(dd)//'_'//trim(overpass)//'.nc'
         
end subroutine create_THySM_filename
