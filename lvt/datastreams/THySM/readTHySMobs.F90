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
! !ROUTINE: readTHySMobs
! \label{readTHySMobs}
!
! !INTERFACE: 
subroutine readTHySMobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod,   only : LVT_get_julss
  use THySM_obsMod,     only : THySMobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine provides the data reader for the THySM
! product. 
! 
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 5 April 2021: Sujay Kumar, Initial Specification
! 
!EOP
  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists1,file_exists2
  integer           :: c,r
  integer           :: ios
  integer           :: ftn
  character*100     :: fname_AM, fname_PM
  real              :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  real              :: sm_file(THySMobs(source)%nc, THySMobs(source)%nr)
  real              :: sm_inp(THySMobs(source)%nc*THySMobs(source)%nr)
  logical*1         :: sm_b_inp(THySMobs(source)%nc*THySMobs(source)%nr)
  real              :: sm_out(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: sm_b_out(LVT_rc%lnc*LVT_rc%lnr)
  real              :: lat(THySMobs(source)%nr)
  real              :: lon(THySMobs(source)%nc)
  integer           :: smid,flagid,latid,lonid


  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  smobs= LVT_rc%udef
  sm_b_inp  = .false. 
  sm_inp = LVT_rc%udef
  
  if(alarmCheck) then

     !AM is descending and PM is ascending
     call create_THySM_filename(THySMobs(source)%odir, &
          LVT_rc%dyr(source), LVT_rc%dmo(source), &
          LVT_rc%dda(source), 'AM',fname_AM)
  
     inquire(file=trim(fname_AM),exist=file_exists1)
     if(file_exists1) then
        
        write(LVT_logunit,*) '[INFO] Reading ',trim(fname_AM)
        
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        ios = nf90_open(path=trim(fname_AM),mode=NF90_NOWRITE,ncid=ftn)
        call LVT_verify(ios,'Error opening file '//trim(fname_AM))
        
        ios = nf90_inq_varid(ftn, 'Band1',smid)
        call LVT_verify(ios, 'Error nf90_inq_varid: Band1')
        
        ios = nf90_get_var(ftn, smid, sm_file)
        call LVT_verify(ios, 'Error nf90_get_var: Band1')
        
        ios = nf90_close(ncid=ftn)
        call LVT_verify(ios,'Error closing file '//trim(fname_AM))
#endif
        write(LVT_logunit,*) '[INFO] Finished reading ',trim(fname_AM)
        
        sm_b_inp  = .false. 
        sm_inp = LVT_rc%udef
        
        do r=1,THySMobs(source)%nr
           do c=1,THySMobs(source)%nc
              if(sm_file(c,r).gt.0) then 
                 sm_b_inp(c+(r-1)*THySMobs(source)%nc) = & 
                      .true.
                 sm_inp(c+(r-1)*THySMobs(source)%nc) = &
                      sm_file(c,r)
              endif
           enddo
        enddo
        
     endif
     
     call create_THySM_filename(THySMobs(source)%odir, &
          LVT_rc%dyr(source), LVT_rc%dmo(source), &
          LVT_rc%dda(source), 'PM',fname_PM)
     
     inquire(file=trim(fname_PM),exist=file_exists2)
     if(file_exists2) then
        
        write(LVT_logunit,*) '[INFO] Reading ',trim(fname_PM)
        
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        ios = nf90_open(path=trim(fname_PM),mode=NF90_NOWRITE,ncid=ftn)
        call LVT_verify(ios,'Error opening file '//trim(fname_PM))
        
        ios = nf90_inq_varid(ftn, 'Band1',smid)
        call LVT_verify(ios, 'Error nf90_inq_varid: Band1')
        
        ios = nf90_get_var(ftn, smid, sm_file)
        call LVT_verify(ios, 'Error nf90_get_var: Band1')
        
        ios = nf90_close(ncid=ftn)
        call LVT_verify(ios,'Error closing file '//trim(fname_PM))
#endif
        write(LVT_logunit,*) '[INFO] Finished reading ',trim(fname_PM)
        
        do r=1,THySMobs(source)%nr
           do c=1,THySMobs(source)%nc
              if(sm_file(c,r).gt.0.and.&
                   sm_inp(c+(r-1)*THySMobs(source)%nc).eq.-9999.0) then 
                 sm_b_inp(c+(r-1)*THySMobs(source)%nc) = & 
                      .true.
                 sm_inp(c+(r-1)*THySMobs(source)%nc) = &
                      sm_file(c,r)
              endif
           enddo
        enddo
        
     endif

     if(file_exists1.or.file_exists2) then 
        call bilinear_interp(LVT_rc%gridDesc(:),&
             sm_b_inp,&
             sm_inp,&
             sm_b_out,&
             sm_out,&
             THySMobs(source)%nc*THySMobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr, &
             THySMobs(source)%rlat,&
             THySMobs(source)%rlon,&
             THySMobs(source)%w11,&
             THySMobs(source)%w12,&
             THySMobs(source)%w21,&
             THySMobs(source)%w22,&
             THySMobs(source)%n11,&
             THySMobs(source)%n12,&
             THySMobs(source)%n21,&
             THySMobs(source)%n22,&
             LVT_rc%udef,ios)
        
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              smobs(c,r) = sm_out(c+(r-1)*LVT_rc%lnc)
           enddo
        enddo
        
     endif
     
  endif
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       smobs,vlevel=1,units="m3/m3")

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
