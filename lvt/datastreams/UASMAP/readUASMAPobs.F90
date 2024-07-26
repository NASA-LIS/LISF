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
! !ROUTINE: readUASMAPobs
! \label{readUASMAPobs}
!
! !INTERFACE: 
subroutine readUASMAPobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod,   only : LVT_get_julss
  use UASMAP_obsMod,     only : UASMAPobs
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
! This subroutine provides the data reader for the UASMAP
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
  logical           :: file_exists
  integer           :: c,r,r1
  integer           :: ios
  integer           :: ftn
  character*100     :: fname
  real              :: smobs(LVT_rc%lnc,LVT_rc%lnr)
  real              :: sm_file(UASMAPobs(source)%nc, UASMAPobs(source)%nr)
  real              :: sm_inp(UASMAPobs(source)%nc*UASMAPobs(source)%nr)
  logical*1         :: sm_b_inp(UASMAPobs(source)%nc*UASMAPobs(source)%nr)
  real              :: sm_out(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: sm_b_out(LVT_rc%lnc*LVT_rc%lnr)
  real              :: lat(UASMAPobs(source)%nr)
  real              :: lon(UASMAPobs(source)%nc)
  integer           :: smid,flagid,latid,lonid


  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  smobs= LVT_rc%udef

  if(alarmCheck) then

     !AM is descending and PM is ascending
     call create_UASMAP_filename(UASMAPobs(source)%odir, &
          LVT_rc%dyr(source), LVT_rc%dmo(source), &
          LVT_rc%dda(source), fname)


     inquire(file=trim(fname),exist=file_exists)
     if(file_exists) then
        
        write(LVT_logunit,*) '[INFO] Reading ',trim(fname)
        
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=ftn)
        call LVT_verify(ios,'Error opening file '//trim(fname))
        
        ios = nf90_inq_varid(ftn, 'SPL3SMP_D',smid)
        call LVT_verify(ios, 'Error nf90_inq_varid: SPL3SMP_D')
        
        ios = nf90_get_var(ftn, smid, sm_file)
        call LVT_verify(ios, 'Error nf90_get_var: SPL3SMP_D')
        
        ios = nf90_close(ncid=ftn)
        call LVT_verify(ios,'Error closing file '//trim(fname))
#endif
        write(LVT_logunit,*) '[INFO] Finished reading ',trim(fname)
        
        sm_b_inp  = .false. 
        sm_inp = LVT_rc%udef
        
        do r=1,UASMAPobs(source)%nr
           do c=1,UASMAPobs(source)%nc
              if(.not.isNaN(sm_file(c,r)).and.sm_file(c,r).gt.0) then
                 r1=UASMAPobs(source)%nr-r+1
                 sm_b_inp(c+(r1-1)*UASMAPobs(source)%nc) = & 
                      .true.
                 sm_inp(c+(r1-1)*UASMAPobs(source)%nc) = &
                      sm_file(c,r)
              endif
           enddo
        enddo
        
     endif
     
     if(file_exists) then 
        call bilinear_interp(LVT_rc%gridDesc(:),&
             sm_b_inp,&
             sm_inp,&
             sm_b_out,&
             sm_out,&
             UASMAPobs(source)%nc*UASMAPobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr, &
             UASMAPobs(source)%rlat,&
             UASMAPobs(source)%rlon,&
             UASMAPobs(source)%w11,&
             UASMAPobs(source)%w12,&
             UASMAPobs(source)%w21,&
             UASMAPobs(source)%w22,&
             UASMAPobs(source)%n11,&
             UASMAPobs(source)%n12,&
             UASMAPobs(source)%n21,&
             UASMAPobs(source)%n22,&
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

end subroutine readUASMAPobs


!BOP
! !ROUTINE: create_UASMAP_filename
! \label{create_UASMAP_filename}
! 
! !INTERFACE: 
subroutine create_UASMAP_filename(odir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: odir

! 
! !DESCRIPTION: 
!  This subroutine creates the timestamped UASMAP filename 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the UASMAP directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated UASMAP filename
! \end{description}
!EOP

  character*4      :: yyyy
  character*2      :: mm,dd

  write(yyyy,'(i4.4)') yr
  write(mm,'(i2.2)') mo
  write(dd,'(i2.2)') da
  
  filename=trim(odir)//'/'//&
       '/SPL3SMP_D_f06_1km_'//trim(yyyy)//&
       trim(mm)//trim(dd)//'.nc'
         
end subroutine create_UASMAP_filename
