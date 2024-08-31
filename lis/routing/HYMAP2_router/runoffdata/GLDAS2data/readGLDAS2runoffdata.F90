!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! 
! !DESCRIPTION:
! This routine is for reading in GLDAS 2.0 data in native NetCDF format.
! 
! !REVISION HISTORY: 
! 30 Jan 2016: Hiroko Beaudoing, Initial implementation
! 31 Aug 2016: Augusto Getirana, Fix file name format
! 
! !USES: 
!
!
!EOP 
subroutine readGLDAS2runoffdata(n,surface_runoff, baseflow)

  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use GLDAS2runoffdataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n))
  integer                       :: nc,nr
  integer                       :: c,r,t
  !real,   allocatable           :: qs(:,:)
  !real,   allocatable           :: qsb(:,:)
  integer                       :: ios, nid,qsid,qsbid
  integer                       :: ftn
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer                       :: doy, yr, mo, da, hr, mn, ss, ts
  real*8                        :: time
  real                          :: gmt
  logical                       :: file_exists
  real                          :: undef


  yr =LIS_rc%yr    !Next Hour
  mo =LIS_rc%mo
  da =LIS_rc%da
!  hr =GLDAS2runoffdata_struc(n)%outInterval*((LIS_rc%hr)/&
!       GLDAS2runoffdata_struc(n)%outInterval)
  hr=LIS_rc%hr-imod(LIS_rc%hr,int(GLDAS2runoffdata_struc(n)%outInterval/3600.))
  print*,int(GLDAS2runoffdata_struc(n)%outInterval/3600.),LIS_rc%hr,hr
  mn =0
  ss =0
  ts =GLDAS2runoffdata_struc(n)%outInterval

  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_GLDAS2_filename(GLDAS2runoffdata_struc(n)%odir,&
       GLDAS2runoffdata_struc(n)%model_name,&
       GLDAS2runoffdata_struc(n)%datares,&
       yr, mo, da, doy, hr, filename)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  if(trim(GLDAS2runoffdata_struc(n)%previous_filename)/=trim(filename))then
  
    GLDAS2runoffdata_struc(n)%previous_filename=filename


  surface_runoff = 0.0
  baseflow       = 0.0

  nc = GLDAS2runoffdata_struc(n)%nc
  nr = GLDAS2runoffdata_struc(n)%nr

!  allocate(qs(nc,nr))
!  allocate(qsb(nc,nr))
  
!  qs = LIS_rc%udef
!  qsb = LIS_rc%udef


  inquire(file=filename, exist=file_exists)
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading '//trim(filename)
        
     ios = nf90_open(path=filename,&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in readGLDAS2runoffdata')

     ios = nf90_inq_varid(nid,'Qs_acc',qsid)
     call LIS_verify(ios,'failed to inquire Qs_acc field in readGLDAS2runoffdata')

     ios = nf90_inq_varid(nid,'Qsb_acc',qsbid)
     call LIS_verify(ios,'failed to inquire Qsb_acc field in readGLDAS2runoffdata')

     ios = nf90_get_var(nid,qsid,GLDAS2runoffdata_struc(n)%qs)
     call LIS_verify(ios, 'failed to read Qs_acc field in readGLDAS2runoffdata')

     ios = nf90_get_var(nid,qsbid,GLDAS2runoffdata_struc(n)%qsb)
     call LIS_verify(ios, 'failed to read Qsb_acc field in readGLDAS2runoffdata')

     ios = nf90_get_att(nid,qsid,'missing_value',undef)
     call LIS_verify(ios,'failed to read missing_value in readGLDAS2runoffdata')

     call LIS_verify(nf90_close(nid))

! convert units from kg m-2 3hr-1 to kg m-2 sec-1.
   if ( GLDAS2runoffdata_struc(n)%outInterval .eq. 10800 ) then
    do r=1, nr
     do c=1, nc
        if ( GLDAS2runoffdata_struc(n)%qs(c,r) .ne. undef ) then
          GLDAS2runoffdata_struc(n)%qs(c,r) = GLDAS2runoffdata_struc(n)%qs(c,r) / GLDAS2runoffdata_struc(n)%outInterval
          GLDAS2runoffdata_struc(n)%qsb(c,r) = GLDAS2runoffdata_struc(n)%qsb(c,r) / GLDAS2runoffdata_struc(n)%outInterval
        endif
     end do
    end do
   else
     write(LIS_logunit,*) 'Runoff input units conversion not supported for this output interval'
     call LIS_endrun()
   endif
  else
     write(LIS_logunit,*) 'Failed to find '//trim(filename)
     call LIS_endrun()
  endif

  call interp_GLDAS2runoffdata(n,nc,nr,GLDAS2runoffdata_struc(n)%qs,surface_runoff)
  call interp_GLDAS2runoffdata(n,nc,nr,GLDAS2runoffdata_struc(n)%qsb,baseflow)

!  deallocate(qs)
!  deallocate(qsb)

  endif  

#endif
    
end subroutine readGLDAS2runoffdata

subroutine readGLDAS2evapdata(n,evap, potevap)

  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use GLDAS2runoffdataMod
  use LIS_fileIOMod
  use HYMAP_routingMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer,          intent(in) :: n
  real                         :: evap(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: potevap(LIS_rc%gnc(n),LIS_rc%gnr(n))
  integer                       :: nc,nr
  integer                       :: c,r,t
  real,   allocatable           :: evp(:,:)
  real,   allocatable           :: potevp(:,:)
  integer                       :: ios, nid,evpid,potevpid
  integer                       :: ftn
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer                       :: doy, yr, mo, da, hr, mn, ss, ts
  real*8                        :: time
  real                          :: gmt
  logical                       :: file_exists
  real                          :: undef
  REAL, PARAMETER:: LIS_CONST_LATVAP = 2.501e6 ! Latent heat for evapo for water in Noah


  yr =LIS_rc%yr    !Next Hour
  mo =LIS_rc%mo
  da =LIS_rc%da
  hr =GLDAS2runoffdata_struc(n)%outInterval*((LIS_rc%hr)/&
       GLDAS2runoffdata_struc(n)%outInterval)
  mn =0
  ss =0
  ts =GLDAS2runoffdata_struc(n)%outInterval

  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_GLDAS2_filename(GLDAS2runoffdata_struc(n)%odir,&
       GLDAS2runoffdata_struc(n)%model_name,&
       GLDAS2runoffdata_struc(n)%datares,&
       yr, mo, da, doy, hr, filename)

  evap = 0.0
  potevap       = 0.0

  nc = GLDAS2runoffdata_struc(n)%nc
  nr = GLDAS2runoffdata_struc(n)%nr

  allocate(evp(nc,nr))
  allocate(potevp(nc,nr))
  
  evp = LIS_rc%udef
  potevp = LIS_rc%udef

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=filename, exist=file_exists)
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading '//trim(filename)
        
     ios = nf90_open(path=filename,&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in readGLDAS2evapdata')

     ios = nf90_inq_varid(nid,'Evap_tavg',evpid)
     call LIS_verify(ios,'failed to inquire Evap_tavg field in readGLDAS2evapdata')

     ios = nf90_inq_varid(nid,'PotEvap_tavg',potevpid)
     call LIS_verify(ios,'failed to inquire PotEvap_tavg field in readGLDAS2evapdata')

     ios = nf90_get_var(nid,evpid,evp)
     call LIS_verify(ios, 'failed to read Evap_tavg field in readGLDAS2evapdata')

     ios = nf90_get_var(nid,potevpid,potevp)
     call LIS_verify(ios, 'failed to read PotEvap_tavg field in readGLDAS2evapdata')

     ios = nf90_get_att(nid,evpid,'missing_value',undef)
     call LIS_verify(ios,'failed to read missing_value in readGLDAS2evapdata')

     call LIS_verify(nf90_close(nid))

! Convert units of Potevap from Wm-2 to kg m-2 sec-1.
! units of evap is already in kg m-2 sec-1.
    do r=1, nr
     do c=1, nc
        if ( potevp(c,r) .ne. undef ) then
          potevp(c,r) = potevp(c,r) / LIS_CONST_LATVAP
        endif
     end do
    end do
  else
     write(LIS_logunit,*) 'Failed to find '//trim(filename)
     call LIS_endrun()
  endif
#endif

  call interp_GLDAS2runoffdata(n,nc,nr,evp,evap)
  call interp_GLDAS2runoffdata(n,nc,nr,potevp,potevap)

  deallocate(evp)
  deallocate(potevp)

    
end subroutine readGLDAS2evapdata
!BOP
! 
! !ROUTINE: create_GLDAS2_filename
! \label{create_GLDAS2_filename}
!
! !INTERFACE: 
subroutine create_GLDAS2_filename(odir,model_name, datares,&
     yr,mo,da, doy,hr,filename)

  use LIS_logMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

! 
! !USES:   
  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  character(len=*)             :: model_name
  real                         :: datares
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  integer                      :: doy
  integer                      :: hr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the GLDAS2 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            GLDAS2 base directory
!   \item[model\_name]     name of the model used in the GLDAS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the GLDAS2 file
!  \end{description}
! 
!EOP
  
  integer                 :: ftn
  character*4             :: fyr
  character*3             :: fdoy
  character*2             :: fmo, fhr, fda
  integer                 :: ierr
  character(len=LIS_CONST_PATH_LEN) :: list_name

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  
  if(datares .eq. 0.25) then   
     !ag - 31Aug2016
     list_name = 'ls '//trim(odir)//'/'//trim(fyr)//trim(fmo)//&
          '/GLDAS_'//trim(model_name)//'025_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '*.020.nc4 > GLDAS2_file'
     filename=trim(odir)//'/'//trim(fyr)//trim(fmo)//&
          '/GLDAS_'//trim(model_name)//'025_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '00.020.nc4'
!     list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
!          '/GLDAS_'//trim(model_name)//'025_3H.A'//&
!          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
!          '*.020.nc4 > GLDAS2_file'
  elseif(datares.eq. 1.0) then 
     list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
          '/GLDAS_'//trim(model_name)//'10_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '*.020.nc4 > GLDAS2_file'
     filename=trim(odir)//'/'//trim(fyr)//trim(fmo)//&
          '/GLDAS_'//trim(model_name)//'025_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '00.020.nc4'

  endif
  
!  call system(trim(list_name))
!  ftn = LIS_getNextUnitNumber()
!  open(ftn,file='GLDAS2_file',status='old',iostat=ierr)
!  do while(ierr.eq.0) 
!     read(ftn,'(a)',iostat=ierr) filename
!     if(ierr.ne.0) then 
!        exit
!     endif
!  enddo
!  call LIS_releaseUnitNumber(ftn)

end subroutine create_GLDAS2_filename

!BOP
!
! !ROUTINE: interp_GLDAS2runoffdata
!  \label{interp_GLDAS2runoffdata}
!
! !INTERFACE:
  subroutine interp_GLDAS2runoffdata(n, nc,nr,var_input,var_output)
!
! !USES:
    use LIS_coreMod
    use GLDAS2runoffdataMod
      
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the GLDAS2 variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (GLDAS2) grid
!    \item[nr]      number of rows in the input (GLDAS2) grid
!    \item[var_input] input variable to be interpolated
!    \item[lb]        input bitmap (true//false)
!    \item[var_output] resulting interpolated field
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
!BOP
!
! !ARGUMENTS:
    integer            :: n
    integer            :: nc
    integer            :: nr
    real               :: var_input(nc*nr)
    logical*1          :: lb(nc*nr)
    real               :: var_output(LIS_rc%lnc(n), LIS_rc%lnr(n))
    !EOP
    integer            :: ios
    integer            :: c,r
    logical*1          :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    real               :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))

    var_output = LIS_rc%udef
    lb = .false.
    do r = 1,nr
       do c = 1,nc
          if (var_input(c+(r-1)*nc).ne.LIS_rc%udef) then
             lb(c+(r-1)*nc) = .true.
          endif
       enddo
    enddo

    if(LIS_isAtAfinerResolution(n,GLDAS2runoffdata_struc(n)%datares)) then
       call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
            GLDAS2runoffdata_struc(n)%n11,                         & 
            LIS_rc%udef,ios)
    else
       call upscaleByAveraging(&
            nc*nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%udef, &
            GLDAS2runoffdata_struc(n)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          var_output(c,r) = go(c+(r-1)*LIS_rc%lnc(n))
       enddo
    enddo

  end subroutine interp_GLDAS2runoffdata


