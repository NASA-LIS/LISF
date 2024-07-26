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
! !ROUTINE: read_nldas2b
!  \label{read_nldas2b}
!
! !REVISION HISTORY:
!  11 Apr 2000: Brian Cosgrove; changed code to use Forcing Mask (With
!               inland water filled in).  Deleted unused variables.
!  27 Apr 2000: Brian Cosgrove; changed code to use the original 
!               mask again since that is the
!               mask which  NCEP has already applied to the forcing data
!               by the time NASA gets it......not possible to use the 
!               expanded NASA forcing mask
!  1  May 2000: Brian Cosgrove; changed code so that if parameter 11 (sw)
!               is not found in hourly ncep data, it will just use
!               edas-based shortwave from the hourly ncep files
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!               not a hard-wired undefined value of -999.9 and -999.0
!  18 Aug 2000: Brian Cosgrove; changed code so that FMASK and not MASK
!               is used when ungribbing.  NCEP data already has a mask applied
!               to it and so may not be able to supply forcing data to
!               all LDAS land forcing points.  In areas where LDAS
!               forcing mask states that land exists, but where NCEP forcing
!               data is non-existant, assign undefined value to forcing data.
!  22 Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
!  05 Sep 2001: Brian Cosgrove; Removed dirnom and infile variables, changed
!               call to ungribncep to match removal.  Added code to make use
!               of precip weighting mask
!  02 Feb 2004: Sujay Kumar; Initial Specification in LIS
!  24 Aug 2007: Chuck Alonge; Modified for use with NLDAS2 data
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  14 Mar 2014: David Mocko: Fixed the use of Forcing "B" file data
!                            of NARR meteorology, radiation, precip,
!                            and/or forcing height fields
! !INTERFACE:
subroutine read_nldas2b(n, kk, findex, order, name,ferror)
! !USES:
  use LIS_coreMod, only        : LIS_rc, LIS_domain
  use LIS_logMod, only         : LIS_logunit, LIS_verify, LIS_warning
  use LIS_metforcingMod, only  : LIS_forc
  use nldas2_forcingMod,only : nldas2_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)       :: n
  integer, intent(in)       :: kk      ! Forecast member index
  integer, intent(in)       :: findex  ! Forcing index
  integer, intent(in)       :: order
  character(len=*), intent(in) :: name
  integer, intent(out)      :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  NLDAS2 data, transforms into 10 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the hourly NLDAS2 forecast file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_nldas2](\ref{interp_nldas2}) \newline
!    spatially interpolates a NLDAS2 variable
!  \end{description}
!EOP

  integer                   :: iv, ftn
  integer                   :: nldas2
  integer                   :: k,t,c,r,iret,rc
  integer                   :: var_index
  real                      :: missingValue 
  integer                   :: nvars
  integer                   :: igrib
  logical                   :: pcp_flag, var_found
  logical                   :: var_status(10)
  logical                   :: file_exists
  integer                   :: pds5(10), pds7(10),pds2(10)
  integer                   :: pds5_val, pds7_val, pds2_val
  logical*1, allocatable    :: lb(:)
  real, allocatable         :: f(:)
  real                      :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))

  ferror = 1
  iv = 0
#if(defined USE_GRIBAPI) 
  nldas2 = (nldas2_struc(n)%ncold*nldas2_struc(n)%nrold)
! Order of variables to be assigned in the NLDAS-2 Forcing "B" files.
! Note that this is NOT the order of the fields in the actual files,
! but the order in which they are assigned to "metdata". - dmm
  pds5 = (/ 011,051,204,179,033,034,001,061,063,007/) !parameter
  pds7 = (/ 001,001,000,000,001,001,001,000,000,001/) !htlev2
  pds2 = (/  84, 84, 84, 84, 84, 84, 84, 84, 84, 84/)

  varfield = 0 
  ferror = 1
  var_status = .false. 

  inquire (file=name, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(name),'r',iret)
     if(iret.ne.0) then 
        write(LIS_logunit,*) &
             '[ERR] Could not open file: ',trim(name)
        ferror = 0
        return
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LIS_verify(iret, 'error in grib_count_in_file in read_nldas2b')

     allocate(lb(nldas2_struc(n)%ncold*nldas2_struc(n)%nrold))
     allocate(f(nldas2_struc(n)%ncold*nldas2_struc(n)%nrold))
     
     do k=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LIS_warning(iret, 'error in grib_new_from_file in read_nldas2b')
        if(iret.ne.0) then 
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(name)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_nldas2b')

        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in read_nldas2b')

        call grib_get(igrib,'generatingProcessIdentifier',pds2_val,rc)
        call LIS_verify(rc, 'error in grib_get: generatingProcessIdentifier in read_nldas2b')

        var_found = .false. 
        do iv=1,nvars
           if((pds5_val.eq.pds5(iv)).and.&
                (pds7_val.eq.pds7(iv)).and.&
                (pds2_val.eq.pds2(iv))) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo

        f = LIS_rc%udef
        call grib_get(igrib,'values',f,rc)
        call LIS_warning(rc, 'error in grib_get:values in read_nldas2b')

        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(name)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in read_nldas2b')

        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in read_nldas2b')
        
        if(var_found) then 
           lb = .false.
           do t=1,nldas2
              if(f(t).ne.missingValue) lb(t) = .true. 
           enddo
           
           pcp_flag = .false. 
           if(var_index.eq.8.or.var_index.eq.9) pcp_flag = .true. 
           
           call interp_nldas2(n, findex,LIS_rc%mo, pcp_flag,nldas2,f,&
                lb,LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield)

           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    
                    ! MODEL LEVEL TAIR CASE
                    if (var_index == 1 ) then 
                       if ( nldas2_struc(n)%model_level_data .gt. 0 ) then 
                          if(order.eq.1) then
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model level data ?
                    end if ! var_index=1
                    
                    ! MODEL LEVEL SPFH CASE
                    if (var_index == 2 ) then 
                       if ( nldas2_struc(n)%model_level_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model level data ?
                    end if ! var_index=2
                    
                    ! USE MODEL BASED DSWRF CASE
                    if (var_index == 3 ) then 
                       if ( nldas2_struc(n)%model_dswrf_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model based dswrf data ?
                    end if ! var_index=3
                    
                    ! MODEL LEVEL USE AERODYNAMIC CONDUCTANCE
                    if (var_index == 4 ) then 
                       if ( nldas2_struc(n)%model_level_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,13,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,13,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model level data ?
                    end if ! var_index=4
                    
                    ! MODEL LEVEL UWIND CASE
                    if (var_index == 5 ) then 
                       if ( nldas2_struc(n)%model_level_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model level data ?
                    end if ! var_index=5
                    
                    ! MODEL LEVEL VWIND CASE
                    if (var_index == 6 ) then 
                       if ( nldas2_struc(n)%model_level_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model level data ?
                    end if ! var_index=6
                    
                    ! MODEL LEVEL PRESSURE CASE
                    if (var_index == 7 ) then 
                       if ( nldas2_struc(n)%model_level_press .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model level pressure ?
                    end if ! var_index=7
                    
                    ! MODEL BASED PRECIP CASE
                    if (var_index == 8 ) then 
                       if ( nldas2_struc(n)%model_pcp_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model pcp data ?
                    end if ! var_index=8
                    
                    ! MODEL BASED CONV. PRECIP CASE
                    if (var_index == 9 ) then 
                       if ( nldas2_struc(n)%model_pcp_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,var_index,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model pcp data ?
                    end if ! var_index=9
                    
                    ! MODEL FORCING HEIGHT CASE
                    if (var_index == 10 ) then 
                       if ( nldas2_struc(n)%model_level_data .gt. 0 ) then 
                          if(order.eq.1) then 
                             nldas2_struc(n)%metdata1(kk,12,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          elseif(order.eq.2) then 
                             nldas2_struc(n)%metdata2(kk,12,&
                                  LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                          endif
                       end if ! model level data?
                    end if ! var_index=10
                    
                 endif
              end do
           enddo

        endif
     enddo

     call grib_close_file(ftn)
     deallocate(f)     
     
     do k=1,nvars
        if(.not.var_status(k)) then 
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(name)
           ferror = 0
           
           return
        endif
     enddo
     deallocate(lb)
  else
     write(LIS_logunit,*) &
          '[ERR] Could not find file: ',trim(name)
     ferror = 0
  endif
#endif     

end subroutine read_nldas2b
