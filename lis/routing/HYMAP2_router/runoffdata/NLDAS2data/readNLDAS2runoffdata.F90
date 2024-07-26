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
subroutine readNLDAS2runoffdata(n,surface_runoff, baseflow)
! !USES: 
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use NLDAS2runoffdataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if ( defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none

  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n))
!
!
! !DESCRIPTION: 
!
!  This subroutine reads the runoff fields from the NLDAS2 data, 
!  subsets to the current model time and conducts the geospatial transform
!  to the LIS grid. 
!  
!  The arguments are: 
!  \begin{description}
!   \item[n]               index of the nest
!   \item[surface_runoff]  surface runoff field generated from the NLDAS2 data
!   \item[baseflow]        baseflow field generated from the NLDAS2 data
!  \end{description}
!EOP 



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
  hr =NLDAS2runoffdata_struc(n)%outInterval*((LIS_rc%hr)/&
       NLDAS2runoffdata_struc(n)%outInterval)
  mn =0
  ss =0
  ts =NLDAS2runoffdata_struc(n)%outInterval

  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_NLDAS2_filename(NLDAS2runoffdata_struc(n)%odir,&
       NLDAS2runoffdata_struc(n)%model_name,&
       yr, mo, da, doy, hr, filename)

  qs_index   = 235
  qsb_index  = 234
  
  surface_runoff = 0.0
  baseflow       = 0.0

  nc = NLDAS2runoffdata_struc(n)%nc
  nr = NLDAS2runoffdata_struc(n)%nr

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
             'grib_new_from_file failed in readNLDAS2runoffdata')
        
        call grib_get(igrib,"indicatorOfParameter",pid(index),ios)
        call LIS_verify(ios,&
             'grib_get failed for indicatorOfParameter in readNLDAS2runoffdata')
        
        call grib_get(igrib, "timeRangeIndicator",tid(index), ios)
        call LIS_verify(ios, &
             'grib_get failed for timeRangeIndicator in readNLDAS2runoffdata')

        if(pid(index).eq.qs_index) then 
           call retrieve_NLDAS2data(igrib, nc,nr,nvars,index,qs)
        endif

        if(pid(index).eq.qsb_index) then 
           call retrieve_NLDAS2data(igrib, nc,nr,nvars,index,qsb)
        endif
     enddo
     
     call grib_close_file(ftn,ios)

     deallocate(pid)
     deallocate(tid)

  endif
  
  call interp_NLDAS2runoffdata(n,nc,nr,qs,surface_runoff)
  call interp_NLDAS2runoffdata(n,nc,nr,qsb,baseflow)

  deallocate(qs)
  deallocate(qsb)


#endif

end subroutine readNLDAS2runoffdata

!BOP
! 
! !ROUTINE: create_NLDAS2_filename
! \label{create_NLDAS2_filename}
!
! !INTERFACE: 
subroutine create_NLDAS2_filename(odir,model_name, &
     yr,mo,da,doy,hr,filename)
! !USES:   
  use LIS_logMod

  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  character(len=*)             :: model_name
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  integer                      :: doy
  integer                      :: hr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the NLDAS2 data
! based on the given date. 
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            NLDAS2 base directory
!   \item[model\_name]     name of the model used in the NLDAS2 run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[da]              month of data
!   \item[doy]             day of the year of data
!   \item[hr]              hour of data
!   \item[filename]        Name of the NLDAS2 file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*3             :: fdoy
  character*2             :: fmo, fda,fhr

 
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fdoy,fmt='(i3.3)') doy
  
  if ((model_name.eq."NOAH").or.&
       (model_name.eq."MOS").or.&
       (model_name.eq."VIC")) then
     filename = trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//      &
          '/NLDAS_'//trim(model_name)//'0125_H.A'//trim(fyr)//      &
          trim(fmo)//trim(fda)//'.'//trim(fhr)//'00.002.grb'
  elseif (model_name.eq."SAC") then
     filename = trim(odir)//'/'//trim(fyr)//'/'//                  &
          trim(fyr)//trim(fmo)//trim(fda)//'/'//trim(fyr)//  &
          trim(fmo)//trim(fda)//trim(fhr)//'.SAC.grb'
  elseif (model_name.eq."SACSM") then
     filename = trim(odir)//'/SM/'//trim(fyr)//'/'//trim(fyr)//    &
          trim(fmo)//trim(fda)//trim(fhr)//'.SAC.gdat'
  endif
    
end subroutine create_NLDAS2_filename


!BOP
!
! !ROUTINE: retrieve_NLDAS2data
! \label{retrieve_NLDAS2data}
!
! !INTERFACE:
  subroutine retrieve_NLDAS2data(igrib,nc,nr,nvars,index,nldas_var)
! !USES:
    use grib_api
    use LIS_logMod, only : LIS_verify
    
    implicit none
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This subroutine sums and counts up for each NLDAS2 variable.
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
    real,    intent(out) :: nldas_var(nc*nr)

    call grib_get(igrib,"values",var(index,:),ios)
    call LIS_verify(ios,                                            &
         'grib_get failed for values in readNLDAS2data')

    do r = 1,nr
       do c = 1,nc
          if (var(index,c+(r-1)*nc).ne.9999.0) then
             nldas_var(c+(r-1)*nc) = var(index,c+(r-1)*nc)
          endif
       enddo
    enddo

  end subroutine retrieve_NLDAS2data



!BOP
!
! !ROUTINE: interp_NLDAS2runoffdata
!  \label{interp_NLDAS2runoffdata}
!
! !INTERFACE:
  subroutine interp_NLDAS2runoffdata(n, nc,nr,var_input,var_output)
!
! !USES:
    use LIS_coreMod
    use NLDAS2runoffdataMod
      
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the NLDAS2 variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (NLDAS2) grid
!    \item[nr]      number of rows in the input (NLDAS2) grid
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

    if(LIS_isAtAfinerResolution(n,0.125)) then
       call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
            NLDAS2runoffdata_struc(n)%n11,                         & 
            LIS_rc%udef,ios)
    else
       call upscaleByAveraging(&
            nc*nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%udef, &
            NLDAS2runoffdata_struc(n)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          if(go(c+(r-1)*LIS_rc%lnc(n)).ne.LIS_rc%udef) then 
             var_output(c,r) = go(c+(r-1)*LIS_rc%lnc(n))/3600.0
          else
             var_output(c,r) = LIS_rc%udef
          endif
       enddo
    enddo

  end subroutine interp_NLDAS2runoffdata


