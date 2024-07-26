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
subroutine readGLDAS1runoffdata(n,surface_runoff, baseflow)

  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use GLDAS1runoffdataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if ( defined USE_GRIBAPI)
  use grib_api
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
  integer                       :: index,igrib,nvars
  integer                       :: ftn
  integer, allocatable          :: pid(:),tid(:)
  integer                       :: qs_index, qsb_index
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer                       :: doy, yr, mo, da, hr, mn, ss, ts
  real*8                        :: time
  real                          :: gmt
  logical                       :: file_exists


  yr =LIS_rc%yr    !Next Hour
  mo =LIS_rc%mo
  da =LIS_rc%da
  hr =GLDAS1runoffdata_struc(n)%outInterval*((LIS_rc%hr)/&
       GLDAS1runoffdata_struc(n)%outInterval)
  mn =0
  ss =0
  ts =GLDAS1runoffdata_struc(n)%outInterval

  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_GLDAS1_filename(GLDAS1runoffdata_struc(n)%odir,&
       GLDAS1runoffdata_struc(n)%model_name,&
       GLDAS1runoffdata_struc(n)%datares,&
       yr, mo, doy, hr, filename)

  qs_index   = 235
  qsb_index  = 234
  
  surface_runoff = 0.0
  baseflow       = 0.0

  nc = GLDAS1runoffdata_struc(n)%nc
  nr = GLDAS1runoffdata_struc(n)%nr

  allocate(qs(nc,nr))
  allocate(qsb(nc,nr))
  
  qs = LIS_rc%udef
  qsb = LIS_rc%udef


#if ( defined USE_GRIBAPI)
  inquire(file=filename, exist=file_exists)
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading '//trim(filename)
        
     call grib_open_file(ftn,trim(filename),'r',ios)
     if(ios.ne.0) then 
        write(LIS_logunit,*) &
             '[ERR] file not opened ',trim(filename)
        call LIS_endrun()
     endif
    
     call grib_count_in_file(ftn,nvars)
     
     allocate(pid(nvars))
     allocate(tid(nvars))
     
     do index=1,nvars
        call grib_new_from_file(ftn,igrib,ios)
        call LIS_verify(ios, &
             'grib_new_from_file failed in readGLDAS1runoffdata')
        
        call grib_get(igrib,"indicatorOfParameter",pid(index),ios)
        call LIS_verify(ios,&
             'grib_get failed for indicatorOfParameter in readGLDAS1runoffdata')
        
        call grib_get(igrib, "timeRangeIndicator",tid(index), ios)
        call LIS_verify(ios, &
             'grib_get failed for timeRangeIndicator in readGLDAS1runoffdata')

        if(pid(index).eq.qs_index) then 
           call retrieve_GLDAS1data(igrib, nc,nr,nvars,index,qs)
        endif

        if(pid(index).eq.qs_index) then 
           call retrieve_GLDAS1data(igrib, nc,nr,nvars,index,qsb)
        endif
     enddo
     
     call grib_close_file(ftn,ios)

     deallocate(pid)
     deallocate(tid)

  endif

  call interp_GLDAS1runoffdata(n,nc,nr,qs,surface_runoff)
  call interp_GLDAS1runoffdata(n,nc,nr,qsb,baseflow)

  deallocate(qs)
  deallocate(qsb)

    
#endif

end subroutine readGLDAS1runoffdata

!BOP
! 
! !ROUTINE: create_GLDAS1_filename
! \label{create_GLDAS1_filename}
!
! !INTERFACE: 
subroutine create_GLDAS1_filename(odir,model_name, datares,&
     yr,mo,doy,hr,filename)

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
  integer                      :: doy
  integer                      :: hr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the GLDAS1 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            GLDAS1 base directory
!   \item[model\_name]     name of the model used in the GLDAS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the GLDAS1 file
!  \end{description}
! 
!EOP
  
  integer                 :: ftn
  character*4             :: fyr
  character*3             :: fdoy
  character*2             :: fmo, fhr
  integer                 :: ierr
  character(len=LIS_CONST_PATH_LEN) :: list_name

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fhr, fmt='(i2.2)') hr
  
  if(datares .eq. 0.25) then   
     list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
          '/GLDAS_'//trim(model_name)//'025SUBP_3H.A'//&
          trim(fyr)//trim(fdoy)//'.'//trim(fhr)//&
          '*grb > GLDAS1_file'
  elseif(datares.eq. 1.0) then 
     if(model_name.ne."VIC") then 
        list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
             '/GLDAS_'//trim(model_name)//'10SUBP_3H.A'//&
             trim(fyr)//trim(fdoy)//'.'//trim(fhr)//&
             '*grb > GLDAS1_file'
     else
        list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
             '/GLDAS_'//trim(model_name)//'10_3H.A'//&
             trim(fyr)//trim(fdoy)//'.'//trim(fhr)//&
             '*grb > GLDAS1_file'
     endif
  endif
  call system(trim(list_name))         
     
  ftn = LIS_getNextUnitNumber()
  open(ftn,file='GLDAS1_file',status='old',iostat=ierr)
  do while(ierr.eq.0) 
     read(ftn,'(a)',iostat=ierr) filename
     if(ierr.ne.0) then 
        exit
     endif
  enddo
  call LIS_releaseUnitNumber(ftn)

end subroutine create_GLDAS1_filename


!BOP
!
! !ROUTINE: retrieve_GLDAS1data
! \label{retrieve_GLDAS1data}
!
! !INTERFACE:
  subroutine retrieve_GLDAS1data(igrib,nc,nr,nvars,index,gldas_var)
! !USES:
#if ( defined USE_GRIBAPI)
    use grib_api
#endif
    use LIS_logMod, only : LIS_verify
    
    implicit none
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine sums and counts up for each GLDAS1 variable.
!
! !REVISION HISTORY:
!
!EOP
!BOP
!
! !ARGUMENTS:
    integer              :: c,r,ios
    integer, intent(in)  :: igrib,nc,nr
    integer, intent(in)  :: nvars,index
    real                 :: var(nvars,nc*nr)
    real,    intent(out) :: gldas_var(nc*nr)

#if ( defined USE_GRIBAPI)
    call grib_get(igrib,"values",var(index,:),ios)
    call LIS_verify(ios,                                            &
         'grib_get failed for values in readGLDAS1data')

    do r = 1,nr
       do c = 1,nc
          if (var(index,c+(r-1)*nc).ne.9999.0) then
             gldas_var(c+(r-1)*nc) = var(index,c+(r-1)*nc)
          endif
       enddo
    enddo
#endif

  end subroutine retrieve_GLDAS1data



!BOP
!
! !ROUTINE: interp_GLDAS1runoffdata
!  \label{interp_GLDAS1runoffdata}
!
! !INTERFACE:
  subroutine interp_GLDAS1runoffdata(n, nc,nr,var_input,var_output)
!
! !USES:
    use LIS_coreMod
    use GLDAS1runoffdataMod
      
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the GLDAS1 variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (GLDAS1) grid
!    \item[nr]      number of rows in the input (GLDAS1) grid
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

    if(LIS_isAtAfinerResolution(n,GLDAS1runoffdata_struc(n)%datares)) then
       call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
            GLDAS1runoffdata_struc(n)%n11,                         & 
            LIS_rc%udef,ios)
    else
       call upscaleByAveraging(&
            nc*nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%udef, &
            GLDAS1runoffdata_struc(n)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          var_output(c,r) = go(c+(r-1)*LIS_rc%lnc(n))
       enddo
    enddo

  end subroutine interp_GLDAS1runoffdata


