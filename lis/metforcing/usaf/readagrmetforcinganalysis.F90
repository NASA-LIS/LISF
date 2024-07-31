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
! !ROUTINE: readagrmetforcinganalysis
! \label{readagrmetforcinganalysis}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
! 21DEC2007, Marv Freimund, Simplify filename creation code
! 13APR2009, Sujay Kumar, Modified for a retrospective mode from AGRMET output
!  2SEP2009, Chris Franks, Modified level and time range for inputs
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
! 13APR2012, Chris Franks, Use instantaneous values and skin temp vs. 2m
!
! !INTERFACE:
subroutine readagrmetforcinganalysis(n,findex, order, agrfile, month)
! !USES:
  use LIS_coreMod, only         : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_timeMgrMod,only       : LIS_get_julhr,LIS_tick,LIS_time2date
  use LIS_logMod, only          : LIS_logunit, LIS_verify, LIS_warning,&
                                  LIS_endrun, LIS_abort, LIS_alert
  use LIS_spatialDownscalingMod
  use AGRMET_forcingMod,   only : agrmet_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,          intent(in) :: n
  integer,          intent(in) :: findex
  integer,          intent(in) :: order
  character(len=*), intent(in) :: agrfile
  integer,          intent(in) :: month
!
! !DESCRIPTION:
!  This routine reads the previously generated surface fields of temperature,
!  relative humidity, wind, and pressure, radiation, and precipitation
!
!  The indices of metdata1 and metdata2 correspond to the following
!  variables:
!  \begin{verbatim}
!   1 - 2m air temp
!   2 - 2m relative humidity
!   3 - shortwave
!   4 - longwave
!   5 - uwind
!   6 - vwind
!   7 - surface pressure
!   8 - precip
!  \end{verbatim}
!
!EOP
  integer                :: ftn
  logical                :: file_exists
  integer                :: iv, iv_total
  integer                :: c,r,t,kk
  integer                :: nagrmet,nvars
  integer                :: count1
  integer                :: igrib
  integer                :: var_index, vid
  logical                :: var_found
  logical                :: var_status(7)
  real                   :: missingValue
  real,      allocatable :: f(:)
  logical*1, allocatable :: lb(:)
  real                   :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical*1              :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                   :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                :: pds5(7), pds7(7), pds6(7),pds16(7)
  integer                :: pds5_val, pds7_val, pds6_val, pds16_val
  integer                :: rc,status,iret
  integer, save          :: step_count=1
  character(len=10)       :: str_count
  character(len=128)     :: dump_name
  character(len=255) :: message(20)
  integer :: mo
#if(defined USE_GRIBAPI) 
!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
!orig
!  pds5 = (/11, 51, 145, 144, 209, 1, 61/)
!  pds6 = (/1,   1,   1,   1,   105, 1,  1/)
!  pds7 = (/0,   0,   0,   0,   10, 0,  0/)
!  pds16 =(/ 1,   1,   1,   1,   1,   1, 133/)  !time avg/instantaneous field

  pds5 = (/ 11,  51, 145, 144, 209,   1, 61/)
  pds6 = (/105, 105,   1,   1, 105,   1,  1/)
  pds7 = (/  2,   2,   0,   0,  10,   0,  0/)
  pds16 =(/  1,   1,   1,   1,   1,   1, 133/)  !time avg/instantaneous field

  nagrmet = (agrmet_struc(n)%ncol*agrmet_struc(n)%nrow)

  var_status(:) = .false. ! EMK Initialize

  iv_total = 7
  inquire (file=agrfile, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(agrfile),'r',iret)
     if(iret.ne.0) then 
        write(LIS_logunit,*) &
             '[ERR] Could not open file: ',trim(agrfile)
        call LIS_endrun
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LIS_verify(iret, 'error in grib_count_in_file in read_agrmet')

     allocate(lb(agrmet_struc(n)%ncol*agrmet_struc(n)%nrow))
     allocate(f(agrmet_struc(n)%ncol*agrmet_struc(n)%nrow))
     
     do kk=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LIS_warning(iret, 'error in grib_new_from_file in read_agrmet')
        if(iret.ne.0) then 
           write(LIS_logunit,*) &
                '[ERR] Error code: ',iret
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(agrfile)
           deallocate(lb)
           deallocate(f)
           call LIS_endrun
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_agrmet')

        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in read_agrmet')

        call grib_get(igrib,'indicatorOfTypeOfLevel',pds6_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfTypeOfLevel in read_agrmet')

        call grib_get(igrib,'timeRangeIndicator',pds16_val,rc)
        call LIS_verify(rc, 'error in grib_get: timeRangeIndicator in read_agrmet')
        
        var_found = .false. 
        do iv=1,iv_total
           if((pds5_val.eq.pds5(iv)).and.&
                (pds7_val.eq.pds7(iv)).and.&
                (pds6_val.eq.pds6(iv)).and.&
                (pds16_val.eq.pds16(iv))) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo
        
        f = -9999.0        
        if(var_found) then 
           call grib_get(igrib,'values',f,rc)
           call LIS_warning(rc, 'error in grib_get:values in read_agrmet')
        endif

        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                '[ERR] Error code: ',rc
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(agrfile)
           write(LIS_logunit,*) 'for variables ',kk
           deallocate(lb)
           deallocate(f)
           call LIS_endrun
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in read_agrmet')

        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in read_agrmet')
        
        if(var_found) then 

           lb = .false.
           do t=1,nagrmet
              if(f(t).ne.missingValue) lb(t) = .true. 
              ! EMK...Bug fix...Should be in loop, and should check var_index,
              ! not kk
              if(var_index.eq.1) then !temperature
                 if(f(t).gt.400.or.f(t).lt.200) then
                    lb(t) = .false.
                 endif
              endif
           enddo
!EMK...This should not be outside of the t do loop, and check var_index!
!           if(kk.eq.1) then !temperature
!              if(f(t).gt.400.or.f(t).lt.200) then
!                 lb(t) = .false.
!              endif
!           endif

           if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 


              if(var_index.eq.7.and.LIS_rc%pcp_downscale(findex).ne.0) then 
                 call LIS_generatePcpClimoRatioField(n,findex,"AGRMET",&
                      month, agrmet_struc(n)%ncol*agrmet_struc(n)%nrow,&
                      f, lb)
              endif

              varfield = -9999.0

              call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,go,&
                   agrmet_struc(n)%ncol*agrmet_struc(n)%nrow, &
                   LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                   LIS_domain(n)%lat, LIS_domain(n)%lon,&
                   agrmet_struc(n)%w11_anl, agrmet_struc(n)%w12_anl, &
                   agrmet_struc(n)%w21_anl, agrmet_struc(n)%w22_anl, &
                   agrmet_struc(n)%n11_anl, agrmet_struc(n)%n12_anl, &
                   agrmet_struc(n)%n21_anl, agrmet_struc(n)%n22_anl, &
                   LIS_rc%udef, iret)

              if(var_index.eq.7.and.LIS_rc%pcp_downscale(findex).ne.0) then 
                 call LIS_pcpClimoDownscaling(n,findex,month, &
                      LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                      go,lo)
                 
              endif
              

!              if(pds5_val .eq. 144 ) then
!                 open(unit=1001, file="lwdn.0.bin_"//trim(str_count), form="unformatted")
!                 write(unit=1001) f
!                 close(unit=1001)
!                 open(unit=1001, file="lwdn.1.bin"//trim(str_count), form="unformatted")
!                 write(unit=1001) go
!                 close(unit=1001)
!              endif  
              count1 = 0 
              do r=1, LIS_rc%lnr(n)
                 do c=1,LIS_rc%lnc(n)
                    varfield(c,r) = go(c+count1)
                 enddo
                 count1 = count1 + LIS_rc%lnc(n)
              enddo
              
              call AGRMET_fillgaps(n,1,varfield)

              if(var_index.eq.1) vid = 1
              if(var_index.eq.2) vid = 2
              if(var_index.eq.3) vid = 3
              if(var_index.eq.4) vid = 4
              if(var_index.eq.5) vid = 5
              if(var_index.eq.6) vid = 6  !pressure
              if(var_index.eq.7) vid = 7  !precip
                            
              do r=1, LIS_rc%lnr(n)
                 do c=1,LIS_rc%lnc(n)
                    if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                       if(order.eq.1) then 

                          agrmet_struc(n)%metdata1(vid,&
                               LIS_domain(n)%gindex(c,r)) = &
                               varfield(c,r)
                       else
                          agrmet_struc(n)%metdata2(vid,&
                               LIS_domain(n)%gindex(c,r)) = &
                               varfield(c,r)
                       endif
                    endif
                 enddo
              enddo
           else if (trim(LIS_rc%met_interp(findex)) .eq. "average") then
              varfield = -9999.0
              mo = LIS_rc%lnc(n) * LIS_rc%lnr(n)
              call upscaleByAveraging(agrmet_struc(n)%mi111, mo, LIS_rc%udef, &
                   agrmet_struc(n)%n111_anl, lb, f, lo, go)
              count1 = 0
              do r = 1, LIS_rc%lnr(n)
                 do c = 1, LIS_rc%lnc(n)
                    varfield(c,r) = go(c+count1)
                 enddo
                 count1 = count1 + LIS_rc%lnc(n)
              enddo
              !call AGRMET_fillgaps(n, 1, varfield)

              if(var_index.eq.1) vid = 1
              if(var_index.eq.2) vid = 2
              if(var_index.eq.3) vid = 3
              if(var_index.eq.4) vid = 4
              if(var_index.eq.5) vid = 5
              if(var_index.eq.6) vid = 6  !pressure
              if(var_index.eq.7) vid = 7  !precip

              do r=1, LIS_rc%lnr(n)
                 do c=1,LIS_rc%lnc(n)
                    if(LIS_domain(n)%gindex(c,r).ne.-1) then
                       if(order.eq.1) then
                          agrmet_struc(n)%metdata1(vid,&
                               LIS_domain(n)%gindex(c,r)) = &
                               varfield(c,r)
                       else
                          agrmet_struc(n)%metdata2(vid,&
                               LIS_domain(n)%gindex(c,r)) = &
                               varfield(c,r)
                       endif
                    endif
                 enddo
              enddo

           else ! EMK Handle error
              write(lis_logunit,*) &
                   '[ERR] Unsupported AGRMET interpolation option ', &
                   trim(LIS_rc%met_interp(findex))
              write(lis_logunit,*)'[ERR] Aborting...'
              message(:) = ''
              message(1) = '[ERR] Program: LIS'
              message(2) = '  Routine: readagrmetforcinganalysis.'
              message(3) = '  Invalid interpolation selected'
              if (LIS_masterproc) then
                 call LIS_alert('LIS.readagrmetforcinganalysis', 1, &
                      message)
                 call LIS_abort(message)
              end if
           endif
        endif
        
     enddo
     call grib_close_file(ftn)

     deallocate(lb)
     deallocate(f)     
         
     do kk=1,iv_total
        if(.not.var_status(kk)) then 
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(agrfile)
           write(LIS_logunit,*) &
                '[ERR] kk,var_status = ',kk,var_status(kk)
           call LIS_endrun
        endif
     enddo
  else
     write(LIS_logunit,*) &
          '[ERR] Could not find file: ',trim(agrfile)
     call LIS_endrun
  endif

#endif
  step_count = step_count+1
end subroutine readagrmetforcinganalysis

