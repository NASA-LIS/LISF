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
subroutine readMERRA2runoffdata(n,surface_runoff, baseflow)

  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use MERRA2runoffdataMod
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
  integer                       :: doy, yr, mo, da, hr, mn, ss, ts
  real*8                        :: time
  real                          :: gmt
  logical                       :: file_exists


  yr =LIS_rc%yr    !Next Hour
  mo =LIS_rc%mo
  da =LIS_rc%da
  hr =MERRA2runoffdata_struc(n)%outInterval*((LIS_rc%hr)/&
       MERRA2runoffdata_struc(n)%outInterval)
  mn =0
  ss =0
  ts =MERRA2runoffdata_struc(n)%outInterval

  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_MERRA2_filename(MERRA2runoffdata_struc(n)%odir,&
       yr, mo, da, filename)

  surface_runoff = 0.0
  baseflow       = 0.0

  nc = MERRA2runoffdata_struc(n)%nc
  nr = MERRA2runoffdata_struc(n)%nr

  tindex = hr+1
  
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
        write(LIS_logunit,*) '[ERR] Opening MERRA2 file ',trim(filename),&
             'failed'
        call LIS_endrun()
     endif
     call LIS_verify(nf90_inq_varid(ftn,'RUNOFF',qsid),&
          'nf90_inq_varid failed for RUNOFF')
     call LIS_verify(nf90_inq_varid(ftn,'BASEFLOW',qsbid),&
          'nf90_inq_varid failed for BASEFLOW')
     
     
     call LIS_verify(nf90_get_var(ftn,qsid,qs,&
          start=(/1,1,tindex/),count=(/nc,nr,1/)),&
          'Error in nf90_get_var RUNOFF')
     call LIS_verify(nf90_get_var(ftn,qsid,qsb,&
          start=(/1,1,tindex/),count=(/nc,nr,1/)),&
          'Error in nf90_get_var BASEFLOW')
     
     call LIS_verify(nf90_close(ftn))

  endif
  
  call interp_MERRA2runoffdata(n,nc,nr,qs,surface_runoff)
  call interp_MERRA2runoffdata(n,nc,nr,qsb,baseflow)

  deallocate(qs)
  deallocate(qsb)

    
end subroutine readMERRA2runoffdata

!BOP
! 
! !ROUTINE: create_MERRA2_filename
! \label{create_MERRA2_filename}
!
! !INTERFACE: 
subroutine create_MERRA2_filename(odir, yr,mo,da,filename)

  use LIS_logMod

! 
! !USES:   
  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the MERRA2 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            MERRA2 base directory
!   \item[model\_name]     name of the model used in the NLDAS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the MERRA2 file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  character*50            :: prefix

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  if (yr==1979 .and. mo>=2) then
     prefix = 'MERRA2_100'
  elseif (yr>1979 .and. yr<=1991) then
     prefix = 'MERRA2_100'
  elseif ( yr >= 1992 .and. yr <= 2000 ) then ! Since 2000 is last full year
     prefix = 'MERRA2_200'
  elseif ( yr >= 2001 .and. yr <= 2009 ) then ! Since 2009 is last full year
     prefix = 'MERRA2_300'
  elseif ( yr >= 2010 ) then
     prefix = 'MERRA2_400'
  else
!     write(LVT_logunit,*) '[ERR] merra2files: date out of range'
!     write(LVT_logunit,*) '[ERR] Supported years are from 1979-2-1 through ...'
!     call LVT_endrun()	
     filename = "none"
  endif
      
  filename = trim(odir)//'/'//trim(prefix)//'/stage/Y'//trim(fyr)//&
       '/M'//trim(fmo)//'/'//trim(prefix)//'.tavg1_2d_lnd_Nx.'//&
       trim(fyr)//trim(fmo)//trim(fda)//'.nc4'

end subroutine create_MERRA2_filename




!BOP
!
! !ROUTINE: interp_MERRA2runoffdata
!  \label{interp_MERRA2runoffdata}
!
! !INTERFACE:
  subroutine interp_MERRA2runoffdata(n, nc,nr,var_input,var_output)
!
! !USES:
    use LIS_coreMod
    use MERRA2runoffdataMod
      
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the MERRA2 variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (MERRA2) grid
!    \item[nr]      number of rows in the input (MERRA2) grid
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

    if(LIS_isAtAfinerResolution(n,0.5)) then
       call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
            MERRA2runoffdata_struc(n)%n11,                         & 
            LIS_rc%udef,ios)
    else
       call upscaleByAveraging(&
            nc*nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%udef, &
            MERRA2runoffdata_struc(n)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          var_output(c,r) = go(c+(r-1)*LIS_rc%lnc(n))/3600.0
       enddo
    enddo

  end subroutine interp_MERRA2runoffdata


