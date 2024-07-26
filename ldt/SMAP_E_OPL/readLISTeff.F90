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
! !ROUTINE: readLIS_Teff
! \label{readLISTeff}
! 
! !REVISION HISTORY: 
!  12 JAN 2022: Yonghwan Kwon, Initial Specification
! 
! !INTERFACE: 
subroutine readLIS_Teff(n,yyyymmdd,hh,Orbit,teff)
! !USES:
  use LDT_coreMod
  use LDT_logMod
  use LDT_smap_e_oplMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
  character*8         :: yyyymmdd
  character*2         :: hh
  character*1         :: Orbit
  real                :: teff(LDT_rc%lnc(n),LDT_rc%lnr(n))
  
!EOP 
  integer             :: c,r
  real                :: tsoil(LDT_rc%lnc(n),LDT_rc%lnr(n),4)
  character*100       :: fname
  logical             :: file_exists
  real                :: kk, cc_6am, cc_6pm    !parameters for calculating effective soil temperature
                                               !from soil layer temperature at 6 AM and 6 PM local time

  teff = LDT_rc%udef
  tsoil = LDT_rc%udef

  call create_LISsoilT_filename(SMAPeOPL%LISdir, &
       yyyymmdd, hh, fname)

  inquire(file=trim(fname),exist=file_exists)  
  if(file_exists) then
     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_LIStsoil_data(n,fname,tsoil)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     kk = 1.007
     cc_6am = 0.246;    !Descending
     cc_6pm = 1.000;    !Ascending

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)

           !calculate effective soil temperature
           if(Orbit.eq."D") then
              if (tsoil(c,r,1).gt.273.15.and.tsoil(c,r,2).gt.273.15) then
                 teff(c,r) = kk * (cc_6am * tsoil(c,r,1) + (1 - cc_6am) * tsoil(c,r,2))
              else
                 teff(c,r) = LDT_rc%udef
              endif
           elseif(Orbit.eq."A") then
              if (tsoil(c,r,1).gt.273.15.and.tsoil(c,r,2).gt.273.15) then
                 teff(c,r) = kk * (cc_6pm * tsoil(c,r,1) + (1 - cc_6pm) * tsoil(c,r,2))
              else
                 teff(c,r) = LDT_rc%udef
              endif
           else
              teff(c,r) = LDT_rc%udef
           endif
        enddo
     enddo
  endif

end subroutine readLIS_Teff

!BOP
! 
! !ROUTINE: read_LIStsoil_data
! \label{read_LIStsoil_data}
!
! !INTERFACE:
subroutine read_LIStsoil_data(n,fname,tsoil)
! 
! !USES:
  use LDT_logMod
  use LDT_coreMod
  use LDT_smap_e_oplMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer, intent(in)           :: n
  character (len=*)             :: fname
!EOP

  integer     :: ios, nid
  integer     :: tsoilid             
  real        :: tsoil(LDT_rc%lnc(n),LDT_rc%lnr(n),4)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))

  ios = nf90_inq_varid(nid,'SoilTemp_inst',tsoilid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SoilTemp_inst')

  !values 

  ios = nf90_get_var(nid, tsoilid, tsoil, &
        start=(/1,1,1/), &
        count=(/LDT_rc%lnc(n),LDT_rc%lnr(n),4/))
  call LDT_verify(ios, 'Error nf90_get_var: tsoil')

  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))

#endif

end subroutine read_LIStsoil_data

!BOP
! !ROUTINE: create_LISsoilT_filename
! \label{create_LISsoilT_filename}
! 
! !INTERFACE:
subroutine create_LISsoilT_filename(LISdir,yyyymmdd,hh,filename)
! !USES:

  implicit none
! !ARGUMENTS:
  character(len=*)  :: filename
  character (len=*) :: LISdir
  character*8       :: yyyymmdd
  character*6       :: yyyymm
  character*2       :: hh
!EOP

  yyyymm = trim(yyyymmdd(1:6))

  filename = trim(LISdir)//'/'//trim(yyyymm)//&
             '/LIS_HIST_'//trim(yyyymmdd)//trim(hh)//'00.d01.nc'

end subroutine create_LISsoilT_filename
