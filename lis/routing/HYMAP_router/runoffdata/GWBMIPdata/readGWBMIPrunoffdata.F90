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
! !REVISION HISTORY: 
! 7 Jan 2016: Sujay Kumar, Initial implementation
! 
! !USES: 
subroutine readGWBMIPrunoffdata(n,surface_runoff, baseflow)

  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use GWBMIPrunoffdataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


!
! !DESCRIPTION: 
!
!EOP 

  implicit none

  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n))
  integer                       :: nc,nr
  integer                       :: c,r,t
  real,   allocatable           :: qs(:,:)
  real,   allocatable           :: qsb(:,:)
  integer                       :: ios
  integer                       :: ftn
  integer                       :: qsid, qsbid
  integer                       :: tindex
  character(len=LIS_CONST_PATH_LEN) :: filename
  type(ESMF_Time)               :: currTime, startTime
  type(ESMF_TimeInterval)       :: deltaT
  integer                       :: da_elapsed, hr_elapsed
  integer                       :: doy, yr, mo, da, hr, mn, ss, ts
  real*8                        :: time
  real                          :: gmt
  logical                       :: file_exists


  yr =LIS_rc%yr    !Next Hour
  mo =LIS_rc%mo
  da =LIS_rc%da
  hr =LIS_rc%hr
  mn =LIS_rc%mn
  ss =LIS_rc%ss
  ts =GWBMIPrunoffdata_struc(n)%outInterval

  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_GWBMIP_filename(GWBMIPrunoffdata_struc(n)%odir,&
       GWBMIPrunoffdata_struc(n)%model_prefix,&
       yr, filename)

  surface_runoff = 0.0
  baseflow       = 0.0

  nc = GWBMIPrunoffdata_struc(n)%nc
  nr = GWBMIPrunoffdata_struc(n)%nr

  call ESMF_TimeSet(startTime,yy=yr, &
       mm = 1, &
       dd = 1, &
       h = 0, &
       m = 0, &
       calendar = LIS_calendar, &
       rc=ios)
  call LIS_verify(ios, 'Error in ESMF_TimeSet: GWBMIPrunoffdata_init')
  
  call ESMF_TimeSet(currTime,yy=yr, &
       mm = mo, &
       dd = da, &
       h = 0, &
       m = 0, &
       s = 0, &
       calendar = LIS_calendar, &
       rc=ios)
  call LIS_verify(ios, 'Error in ESMF_TimeSet: GWBMIPrunoffdata_init')

  deltaT = currTime - startTime
  call ESMF_TimeIntervalGet(deltaT, d=da_elapsed,h=hr_elapsed)

  tindex = da_elapsed + 1
  allocate(qs(nc,nr))
  allocate(qsb(nc,nr))
  
  qs = LIS_rc%udef
  qsb = LIS_rc%udef
  
  inquire(file=filename, exist=file_exists)
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading '//trim(filename)
     ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(ios.ne.0) then
        write(LIS_logunit,*) '[ERR] Opening GWBMIP file ',trim(filename),&
             'failed'
        call LIS_endrun()
     endif
     call LIS_verify(nf90_inq_varid(ftn,'Qs',qsid),&
          'nf90_inq_varid failed for Qs')
     call LIS_verify(nf90_inq_varid(ftn,'Qsb',qsbid),&
          'nf90_inq_varid failed for Qsb')
     
     
     call LIS_verify(nf90_get_var(ftn,qsid,qs,&
          start=(/1,1,tindex/),count=(/nc,nr,1/)),&
          'Error in nf90_get_var Qs')
     call LIS_verify(nf90_get_var(ftn,qsid,qsb,&
          start=(/1,1,tindex/),count=(/nc,nr,1/)),&
          'Error in nf90_get_var Qsb')
     
     call LIS_verify(nf90_close(ftn))

  endif
  
  call interp_GWBMIPrunoffdata(n,nc,nr,qs,surface_runoff)
  call interp_GWBMIPrunoffdata(n,nc,nr,qsb,baseflow)

  deallocate(qs)
  deallocate(qsb)

    
end subroutine readGWBMIPrunoffdata

!BOP
! 
! !ROUTINE: create_GWBMIP_filename
! \label{create_GWBMIP_filename}
!
! !INTERFACE: 
subroutine create_GWBMIP_filename(odir, model_prefix, yr,filename)

  use LIS_logMod

! 
! !USES:   
  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  character(len=*)             :: model_prefix
  integer                      :: yr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the GWBMIP data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            GWBMIP base directory
!   \item[model\_name]     name of the model used in the NLDAS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the GWBMIP file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/'//trim(model_prefix)//'_'//trim(fyr)//'.nc'

end subroutine create_GWBMIP_filename




!BOP
!
! !ROUTINE: interp_GWBMIPrunoffdata
!  \label{interp_GWBMIPrunoffdata}
!
! !INTERFACE:
  subroutine interp_GWBMIPrunoffdata(n, nc,nr,var_input,var_output)
!
! !USES:
    use LIS_coreMod
    use GWBMIPrunoffdataMod
      
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the GWBMIP variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (GWBMIP) grid
!    \item[nr]      number of rows in the input (GWBMIP) grid
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
          if (var_input(c+(r-1)*nc).ne.1e+20) then
             lb(c+(r-1)*nc) = .true.
          endif
       enddo
    enddo

    if(LIS_isAtAfinerResolution(n,1.0)) then
!       call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
!            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
!            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
!            GWBMIPrunoffdata_struc(n)%n11,                         & 
!            LIS_rc%udef,ios)

       call bilinear_interp(LIS_rc%gridDesc,lb,var_input,  &
            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
            GWBMIPrunoffdata_struc(n)%w11,         & 
            GWBMIPrunoffdata_struc(n)%w12,         & 
            GWBMIPrunoffdata_struc(n)%w21,         & 
            GWBMIPrunoffdata_struc(n)%w22,         & 
            GWBMIPrunoffdata_struc(n)%n11,         & 
            GWBMIPrunoffdata_struc(n)%n12,         & 
            GWBMIPrunoffdata_struc(n)%n21,         & 
            GWBMIPrunoffdata_struc(n)%n22,         & 
            LIS_rc%udef,ios)

    else
       call upscaleByAveraging(&
            nc*nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%udef, &
            GWBMIPrunoffdata_struc(n)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          var_output(c,r) = go(c+(r-1)*LIS_rc%lnc(n))
       enddo
    enddo
    
    open(100,file='test.bin',form='unformatted')
    write(100) var_output
    close(100)
    stop

  end subroutine interp_GWBMIPrunoffdata


