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
! !ROUTINE: read_gdasT1534
! \label{read_gdasT1534}
!
! !REVISION HISTORY:
!  20 June 2014: Sujay Kumar; initial implementation
! !INTERFACE:
subroutine read_gdasT1534( order, n, findex, &
     name, ferror,try )
! !USES:  
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_get_nstep, LIS_date2time
  use LIS_metforcingMod,  only : LIS_forc
  use gdasT1534_forcingMod,    only : gdasT1534_struc
  use LIS_logMod,         only : LIS_logunit

  implicit none
! !ARGUMENTS:
  integer, intent(in)          :: order    
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: name
  integer, intent(out)         :: ferror 
  integer, intent(inout)       :: try
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GDAST1534 forecast datasets, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read for the previous 
!    3hr bookend, order=2, read for the next 3 hr bookend)
!  \item[n]
!    index of the nest
!  \item[name00]
!    name of the instantaneous forecast file
!  \item[name03]
!    name of the 3 hour GDAST1534 forecast file
!  \item[name06]
!    name of the 6 hour GDAST1534 forecast file
!  \item[F06flag]
!    flag to indicate if 6hr forecast data is required for this interval
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \item[try]
!    index of the tries (in case of missing data)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gdasT1534](\ref{interp_gdasT1534}) \newline
!    spatially interpolates a GDAST1534 variable
!  \end{description}
!EOP
!==== Local Variables=======================
  
  integer :: iv, c,r,t
  integer :: ferror1, ferror2, ferror3
  integer :: ngdasT1534
  real    :: glbdata(gdasT1534_struc(n)%nmif,LIS_rc%ngrid(n))
  integer :: nstep

!=== End Variable Definition =======================

  glbdata = LIS_rc%udef
  ngdasT1534 = (gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold)
!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
  nstep = LIS_get_nstep(LIS_rc,n)

!--------------------------------------------------------------------------
! if there's a problem then ferror is set to zero
! read instantaneous fields
!--------------------------------------------------------------------------

  iv = 0

!--------------------------------------------------------------------------
! Set up to open file and retrieve specified field 
!--------------------------------------------------------------------------
  call retrieve_gdasT1534_variables(n, findex, name,glbdata, ferror)

!--------------------------------------------------------------------------
! Place the interpolated data into the LIS arrays
!--------------------------------------------------------------------------

  do iv=1,gdasT1534_struc(n)%nmif
     do t=1,LIS_rc%ngrid(n)
        if(order.eq.1) then 
           gdasT1534_struc(n)%metdata1(iv,t) = glbdata(iv,t)    
        else
           gdasT1534_struc(n)%metdata2(iv,t) = glbdata(iv,t)    
        endif
     enddo
  enddo
  return

end subroutine read_gdasT1534

!BOP
! 
! !ROUTINE: retrieve_gdasT1534_variables
! \label{retrieve_gdasT1534_variables}
! 
! !INTERFACE: 
subroutine retrieve_gdasT1534_variables(n, findex, fname, glbdata, errorcode)
! !USES: 
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit,LIS_getNextUnitNumber,& 
       LIS_releaseUnitNumber, LIS_verify, LIS_warning
  use gdasT1534_forcingMod,    only : gdasT1534_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer               :: n 
  integer               :: findex
  character(len=*)      :: fname
  real                  :: glbdata(gdasT1534_struc(n)%nmif,LIS_rc%ngrid(n))
  integer               :: errorcode
! 
! !DESCRIPTION: 
!   This subroutine retrieves GDAST1534 forcing variables, and interpolates
!  them to the LIS grid. 
! 
!EOP

  integer               :: ngdasT1534
  real, allocatable :: f(:)
  real, dimension(LIS_rc%lnc(n), LIS_rc%lnr(n)) :: varfield
  integer :: igrib
  integer :: iv,c,r,t
  real    :: missingValue 
  integer :: iret
  integer :: ftn 
  integer, dimension(gdasT1534_struc(n)%nmif) :: pds5, pds7, pds6,pds16
  integer :: pds5_val, pds7_val, pds16_val, pds6_val
  logical*1, allocatable :: lb(:)
  logical :: file_exists
  integer :: kk
  integer :: var_index
  integer :: nvars
  integer :: rc
  logical :: var_status(gdasT1534_struc(n)%nmif)
  logical :: pcp_flag, var_found
  integer :: grid_size

#if(defined USE_GRIBAPI) 
  pds5 = (/ 011,051,204,205,033,034,001,059,087,084,007,083,208,011,065,066 /) !parameter
  pds6 = (/ 109,109,001,001,109,109,001,001,001,001,109,001,001,001,001,001 /) !level
  pds7 = (/ 001,001,000,000,001,001,000,000,000,000,001,000,000,000,000,000 /) !height
! index 10 indicates instantaneous, 003 indicates time average
  pds16 = (/010,010,010,010,010,010,010,003,010,003,010,010,010,010,010,010 /) !3-ave; 10-fcst

  ngdasT1534 = (gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold)

  varfield = 0 
  errorcode = 1
  var_status = .false. 

  inquire (file=fname, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(fname),'r',iret)
     if(iret.ne.0) then 
        write(LIS_logunit,*) 'ERR: Could not open file: ',trim(fname)
        errorcode = 0
        return
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LIS_verify(iret, 'error in grib_count_in_file in read_gdasT1534')

     allocate(lb(gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold))
     allocate(f(gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold))
     
     do kk=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LIS_warning(iret, 'error in grib_new_from_file in read_gdasT1534')
        if(iret.ne.0) then 
           write(LIS_logunit,*) &

                'Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        ! Trap the old "Could not find correct forcing parameter in file"
        ! error from LIS 6.  This error occurred right before a GDAST1534
        ! grid change.  LIS would try to read ahead, but the new data
        ! would be on the new grid, so LIS would misread it resulting in
        ! the above error message.  The LIS would roll back to the previous
        ! day for GDAST1534 forcing.
        ! Trap this by checking the number of values in one of the
        ! GRIB fields.
        call grib_get_size(igrib,'values',grid_size)
        if ( grid_size /= ngdasT1534 ) then
           write(LIS_logunit,*) &
              'ERR: Number of values does not match expected', trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_gdasT1534')
        call grib_get(igrib,'indicatorOfTypeOfLevel',pds6_val,rc)
        call LIS_verify(rc,'error in grib_get: indicatorOfTypeOfLevel in read_gdasT1534')
        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in read_gdasT1534')

        call grib_get(igrib,'timeRangeIndicator',pds16_val,rc)
        call LIS_verify(rc, 'error in grib_get: timeRangeIndicator in read_gdasT1534')

        var_found = .false. 
        do iv=1,gdasT1534_struc(n)%nmif
           if((pds5_val.eq.pds5(iv)).and.&
                (pds6_val.eq.pds6(iv)).and.&
                (pds7_val.eq.pds7(iv)).and.&
                (pds16_val.eq.pds16(iv))) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo
        f = LIS_rc%udef
        call grib_get(igrib,'values',f,rc)
        call LIS_warning(rc, 'error in grib_get:values in read_gdasT1534')

        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)

           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in read_gdasT1534')

        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in read_gdasT1534')
        
        if(var_found) then 
           lb = .false.
           do t=1,ngdasT1534
              if(f(t).ne.missingValue) lb(t) = .true. 
           enddo
           
           pcp_flag = .false. 
           if(var_index.eq.8) pcp_flag = .true.
!           if(var_index.eq.1) pcp_flag = .true.
           
           call interp_gdasT1534(n, findex,pcp_flag,ngdasT1534,f,&
                lb,LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield)
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    glbdata(var_index,LIS_domain(n)%gindex(c,r)) =&
                         varfield(c,r)
                 endif
              enddo
           enddo
        endif

     enddo
     call grib_close_file(ftn)

     deallocate(lb)
     deallocate(f)     
         
     do kk=1,gdasT1534_struc(n)%nmif
        if(.not.var_status(kk)) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           return
        endif
     enddo
  else
     write(LIS_logunit,*) &
          'Could not find file: ',trim(fname)
     errorcode = 0
  endif
#endif     
end subroutine retrieve_gdasT1534_variables
!BOP
! !ROUTINE: interp_gdasT1534
! \label{interp_gdasT1534}
!
! !INTERFACE:
subroutine interp_gdasT1534(n, findex, pcp_flag, ngdasT1534,f,lb,lis_gds,nc,nr, &
     varfield)
! !USES:
  use LIS_coreMod
  use gdasT1534_forcingMod,   only : gdasT1534_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  logical, intent(in) :: pcp_flag
  integer, intent(in) :: ngdasT1534
  real, intent(out)   :: f(ngdasT1534)
  logical*1           :: lb(ngdasT1534)
  real                :: lis_gds(50)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given GDAST1534 field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[ngdasT1534]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: count1,i,j,mo,c,r

  real, dimension(nc*nr) :: lis1d

  logical*1 :: lo(nc*nr)

!=== End variable declarations
!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------  
  if(LIS_rc%met_interp(findex).eq."bilinear") then 
     call bilinear_interp(lis_gds,lb,f,lo,lis1d,gdasT1534_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
          gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221,&
          gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
          gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,LIS_rc%udef, iret)
  elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 
     if (pcp_flag) then 
        call conserv_interp(lis_gds,lb,f,lo,lis1d,gdasT1534_struc(n)%mi,mo, & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdasT1534_struc(n)%w112,gdasT1534_struc(n)%w122,&
             gdasT1534_struc(n)%w212,gdasT1534_struc(n)%w222,&
             gdasT1534_struc(n)%n112,gdasT1534_struc(n)%n122,&
             gdasT1534_struc(n)%n212,gdasT1534_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,lb,f,lo,lis1d,gdasT1534_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
             gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221,&
             gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
             gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,LIS_rc%udef,iret)
     endif
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
     call neighbor_interp(lis_gds,lb,f,&
          lo,lis1d,gdasT1534_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdasT1534_struc(n)%n113,LIS_rc%udef,iret)     
  endif

!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GDAST1534 & LDAS. For LDAS land 
! points not included in GDAST1534 geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0 
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gdasT1534
