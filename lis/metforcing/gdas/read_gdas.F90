!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_gdas
! \label{read_gdas}
!
! !REVISION HISTORY:
!  14 Dec 2000: Urszula Jambor; Rewrote geteta.f to use GDAS forcing in GLDAS
!  15 Mar 2001: Jon Gottschalck; Added additional parameters and octets in 
!               which to search in GRIB files
!  01 Jun 2001: Urszula Jambor; Added option to get forcing from different 
!               files (F00 instantaneous and F06 six hour means)
!  29 Jan 2003: Urszula Jambor; Rewrote code, uses GETGB call to replace
!               ungribgdas.  Interpolation now occurs in interp_gdas.  
!               Using GETGB avoids problems with the Oct2002 GDAS 
!               grid update.
!  12 Nov 2003: Matt Rodell; Check to make sure input file exists before
!		opening and thereby creating a new, empty file.
!  14 Nov 2003: Matt Rodell; Ensure lugb varies in call to baopen
!  05 Feb 2004: James Geiger; Added GrADS-DODS Server functionality
!  29 Apr 2010: Sujay Kumar: Fixed the mixing of instantaneous and time 
!               averaged fields.
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!
! !INTERFACE:
subroutine read_gdas( order, n, findex, &
     name00, name03, name06, F06flag, ferror,try )
! !USES:  
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_get_nstep, LIS_date2time
  use LIS_metforcingMod,  only : LIS_forc
  use gdas_forcingMod,    only : gdas_struc
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_surfaceModelDataMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)          :: order    
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: name00
  character(len=*), intent(in) :: name03
  character(len=*), intent(in) :: name06
  logical, intent(in)          :: F06flag
  integer, intent(out)         :: ferror 
  integer, intent(inout)       :: try
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GDAS forecast datasets, transforms into 9 LIS forcing 
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
!    name of the 3 hour GDAS forecast file
!  \item[name06]
!    name of the 6 hour GDAS forecast file
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
!  \item[interp\_gdas](\ref{interp_gdas}) \newline
!    spatially interpolates a GDAS variable
!  \end{description}
!EOP
!==== Local Variables=======================
  
  character(len=80) :: fname
  integer :: iv, c,r,t
  integer :: ferror1, ferror2, ferror3
  integer :: nforce
  integer :: ngdas
  real    :: glbdata_i(10,LIS_rc%ngrid(n))
  real    :: glbdata_a(10,LIS_rc%ngrid(n))
  real    :: glbdata_a_f06(10,LIS_rc%ngrid(n))
  integer :: nstep
  logical :: dataStrucflag
!=== End Variable Definition =======================
  glbdata_i = LIS_rc%udef
  glbdata_a = LIS_rc%udef
  glbdata_a_f06 = LIS_rc%udef
  ngdas = (gdas_struc(n)%ncold*gdas_struc(n)%nrold)
  dataStrucflag = .false.
!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
  nstep = LIS_get_nstep(LIS_rc,n)

  nforce = gdas_struc(n)%nmif

!--------------------------------------------------------------------------
! Check model timestep and output interval and stop if the timestep
! and/or output interval will cause incorrect output to be generated.
!--------------------------------------------------------------------------
  if(LIS_rc%ts .ge. 10800 .AND. LIS_rc%ts .lt. 21600) then
     write(LIS_logunit,*) '[WARN] Model timestep is greater than or equal'
     write(LIS_logunit,*) '[WARN]   to 3 hours, which will cause errors in the '
     write(LIS_logunit,*) '[WARN]   GDAS reader. Change the model timestep to a '
     write(LIS_logunit,*) '[WARN]   value less than 3 hours (1hr is suggested.)'
     call LIS_endrun()
 endif
  if(LIS_sfmodel_struc(n)%outInterval .ge. 21600 .AND. LIS_rc%ts .eq. 21600) then
     write(LIS_logunit,*) '[WARN] Model timestep is 6hr and output interval is greater than'
     write(LIS_logunit,*) '[WARN]   or equal to 6hr. This setup can cause issues in the reader'
     write(LIS_logunit,*) '[WARN]   where the output is not the true 6hr average and/or the '
     write(LIS_logunit,*) '[WARN]   data being written is shifted by one timestep.'
     write(LIS_logunit,*) '[WARN] It is suggested that the model timestep be changed to 1hr or less.'
     call LIS_endrun()
  endif

!--------------------------------------------------------------------------
! if there's a problem then ferror is set to zero
! read instantaneous fields
!--------------------------------------------------------------------------

  iv = 0

!--------------------------------------------------------------------------
! Set up to open file and retrieve specified field 
!--------------------------------------------------------------------------
  fname = name00
  if(gdas_struc(n)%dstrucchange1 .AND.  gdas_struc(n)%gdastime1 .ge. gdas_struc(n)%datastructime1) then
    dataStrucflag = .true.  !HKB Use special routine for f00 files following 2019 Jun 12 12Z GDAS upgrades
  endif
  call retrieve_gdas_variables(n, findex, fname, dataStrucflag, glbdata_i, ferror1)
  dataStrucflag = .false. !Reset flag since f03 and f06 files are not affected by 2019 Jun 12 upgrade

!--------------------------------------------------------------------------
! read 3hr forecast for time averaged fields
!--------------------------------------------------------------------------
  fname = name03
  call retrieve_gdas_variables(n, findex, fname, dataStrucflag, glbdata_a, ferror2)

!--------------------------------------------------------------------------
! read 6hr forecast for time averaged fields, if required. 
!--------------------------------------------------------------------------

  if(F06flag) then 
     fname = name06
     call retrieve_gdas_variables(n, findex, fname, dataStrucflag, glbdata_a_f06, ferror3)
  end if
  
  ferror = 1
  if(ferror1.eq.0.or.ferror2.eq.0) then 
     ferror = 0
  endif
  if(F06flag) then
     if(ferror.eq.0.or.ferror3.eq.0) then 
        ferror = 0
     endif
  endif
!--------------------------------------------------------------------------
! Place the interpolated data into the LIS arrays
!--------------------------------------------------------------------------

  do iv=1,9
     do t=1,LIS_rc%ngrid(n)
        if(F06flag) then 
           if(iv.eq.3.or.iv.eq.4.or.iv.eq.8.or.iv.eq.9) then! these are time avgd fields
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = 2*glbdata_a_F06(iv,t)-glbdata_a(iv,t)
              else
                 gdas_struc(n)%metdata2(iv,t) = 2*glbdata_a_F06(iv,t)-glbdata_a(iv,t)
              endif
           else ! these are instantaneous
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = glbdata_i(iv,t)
              else
                 gdas_struc(n)%metdata2(iv,t) = glbdata_i(iv,t)
              endif
           endif
        else
           if(iv.eq.3.or.iv.eq.4.or.iv.eq.8.or.iv.eq.9) then! these are time avgd fields
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = glbdata_a(iv,t)                         
              else
                 gdas_struc(n)%metdata2(iv,t) = glbdata_a(iv,t)                         
              endif
           else
              if(order.eq.1) then 
                 gdas_struc(n)%metdata1(iv,t) = glbdata_i(iv,t)                         
              else
                 gdas_struc(n)%metdata2(iv,t) = glbdata_i(iv,t)                         
              endif
           endif
        endif
     enddo
  enddo

  return

end subroutine read_gdas

!BOP
! 
! !ROUTINE: retrieve_gdas_variables
! \label{retrieve_gdas_variables}
! 
! !INTERFACE: 
subroutine retrieve_gdas_variables(n, findex, fname, dataStrucflag, glbdata, errorcode)
! !USES: 
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit,LIS_getNextUnitNumber,& 
       LIS_releaseUnitNumber, LIS_verify, LIS_warning
  use gdas_forcingMod,    only : gdas_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer               :: n 
  integer               :: findex
  character(len=*)      :: fname
  logical               :: dataStrucflag
  real                  :: glbdata(10,LIS_rc%ngrid(n))
  integer               :: errorcode
! 
! !DESCRIPTION: 
!   This subroutine retrieves GDAS forcing variables, and interpolates
!  them to the LIS grid. 
! 
!EOP

  integer               :: ngdas
  real, allocatable :: f(:)
  real, dimension(LIS_rc%lnc(n), LIS_rc%lnr(n)) :: varfield
  integer :: igrib
  integer :: iv,ivmax,c,r,t
  real    :: missingValue 
  integer :: iret
  integer :: ftn 
  integer, dimension(gdas_struc(n)%nmif) :: pds5, pds7, pds6,pds16
  integer :: pds5_val, pds7_val, pds16_val
  logical*1, allocatable :: lb(:)
  logical :: file_exists
  integer :: kk
  integer :: var_index
  integer :: nvars
  integer :: rc
  logical :: var_status(gdas_struc(n)%nmif)
  logical :: pcp_flag, var_found
  integer :: grid_size

#if(defined USE_GRIBAPI) 
!  pds5 = (/ 011,051,204,205,033,034,001,059,214,084 /) !parameter
!  pds6 = -1
!  pds7 = (/ 002,002,000,000,010,010,000,000,000,000 /) !htlev2
!! index 10 indicates instantaneous, 003 indicates time average
!  pds16 = (/010,010,003,003,010,010,010,003,003,003 /) 

! EMK...Drop last entry (for albedo); not used, and it exceeds array dimension
  pds5 = (/ 011,051,204,205,033,034,001,059,214 /) !parameter
  pds6 = -1
  pds7 = (/ 002,002,000,000,010,010,000,000,000 /) !htlev2

  if(dataStrucflag) then
    ! HKB...All instantaneous fields in f00 files
    pds16 = (/010,010,010,010,010,010,010,010,010 /)
    ivmax = 7
  else
    ! index 10 indicates instantaneous, 003 indicates time average
    pds16 = (/010,010,003,003,010,010,010,003,003 /) 
    ivmax = 9
  endif

  ngdas = (gdas_struc(n)%ncold*gdas_struc(n)%nrold)

  varfield = 0 
  errorcode = 1
  var_status = .false. 

  inquire (file=fname, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(fname),'r',iret)
     if(iret.ne.0) then 
        write(LIS_logunit,*) '[ERR] Could not open file: ',trim(fname)
        errorcode = 0
        return
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LIS_warning(iret, 'error in grib_count_in_file in read_gdas')
     if(iret.ne.0) then 
        errorcode = 0
        return 
     endif

     allocate(lb(gdas_struc(n)%ncold*gdas_struc(n)%nrold))
     allocate(f(gdas_struc(n)%ncold*gdas_struc(n)%nrold))
     
     do kk=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LIS_warning(iret, 'error in grib_new_from_file in read_gdas')
        if(iret.ne.0) then 
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        ! Trap the old "Could not find correct forcing parameter in file"
        ! error from LIS 6.  This error occurred right before a GDAS
        ! grid change.  LIS would try to read ahead, but the new data
        ! would be on the new grid, so LIS would misread it resulting in
        ! the above error message.  The LIS would roll back to the previous
        ! day for GDAS forcing.
        ! Trap this by checking the number of values in one of the
        ! GRIB fields.
        call grib_get_size(igrib,'values',grid_size)
        if ( grid_size /= ngdas ) then
           write(LIS_logunit,*) &
              '[ERR] Number of values does not match expected', trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_gdas')

        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in read_gdas')

        call grib_get(igrib,'timeRangeIndicator',pds16_val,rc)
        call LIS_verify(rc, 'error in grib_get: timeRangeIndicator in read_gdas')

        var_found = .false. 
        do iv=1,ivmax
           if((pds5_val.eq.pds5(iv)).and.&
                (pds7_val.eq.pds7(iv)).and.&
                (pds16_val.eq.pds16(iv))) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo
        f = -9999.0
        call grib_get(igrib,'values',f,rc)
        call LIS_warning(rc, 'error in grib_get:values in read_gdas')

        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in read_gdas')

        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in read_gdas')
        
        if(var_found) then 
           lb = .false.
           do t=1,ngdas
              if(f(t).ne.missingValue) lb(t) = .true. 
           enddo
           
           pcp_flag = .false. 
           if(var_index.eq.8.or.var_index.eq.9) pcp_flag = .true. 
           
           call interp_gdas(n, findex,pcp_flag,ngdas,f,&
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
         
     do kk=1,ivmax
        if(.not.var_status(kk)) then 
           write(LIS_logunit,*) &
                '[ERR] Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           
           return
        endif
     enddo
  else
     write(LIS_logunit,*) &
          '[ERR] Could not find file: ',trim(fname)
     errorcode = 0
  endif
#endif     
end subroutine retrieve_gdas_variables
!BOP
! !ROUTINE: interp_gdas
! \label{interp_gdas}
!
! !INTERFACE:
subroutine interp_gdas(n, findex, pcp_flag, ngdas,f,lb,lis_gds,nc,nr, &
     varfield)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use gdas_forcingMod,   only : gdas_struc
  use LIS_logMod,        only : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  logical, intent(in) :: pcp_flag
  integer, intent(in) :: ngdas
  real, intent(out)   :: f(ngdas)
  logical*1           :: lb(ngdas)
  real                :: lis_gds(50)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given GDAS field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[ngdas]
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
!   \item[upscaleByAveraging](\ref{upscaleByAveraging}) \newline
!     upscales finer resolution forcing data to coarser resolution running
!     domain by averaging
! \end{description}
!EOP
  integer :: iret
  integer :: count1,i,j,mo

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
  if ( gdas_struc(n)%met_interp == "bilinear" ) then
     call bilinear_interp(lis_gds,lb,f,lo,lis1d,gdas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gdas_struc(n)%w111,gdas_struc(n)%w121,&
          gdas_struc(n)%w211,gdas_struc(n)%w221,&
          gdas_struc(n)%n111,gdas_struc(n)%n121,&
          gdas_struc(n)%n211,gdas_struc(n)%n221,LIS_rc%udef, iret)
  elseif ( gdas_struc(n)%met_interp == "budget-bilinear" ) then
     if (pcp_flag) then 
        call conserv_interp(lis_gds,lb,f,lo,lis1d,gdas_struc(n)%mi,mo, & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdas_struc(n)%w112,gdas_struc(n)%w122,&
             gdas_struc(n)%w212,gdas_struc(n)%w222,&
             gdas_struc(n)%n112,gdas_struc(n)%n122,&
             gdas_struc(n)%n212,gdas_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,lb,f,lo,lis1d,gdas_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gdas_struc(n)%w111,gdas_struc(n)%w121,&
             gdas_struc(n)%w211,gdas_struc(n)%w221,&
             gdas_struc(n)%n111,gdas_struc(n)%n121,&
             gdas_struc(n)%n211,gdas_struc(n)%n221,LIS_rc%udef,iret)
     endif
  elseif ( gdas_struc(n)%met_interp == "average" ) then
    call upscaleByAveraging(gdas_struc(n)%mi, mo, LIS_rc%udef, &
         gdas_struc(n)%n111, lb, f, lo, lis1d)
  else
     write(LIS_logunit,*) '[ERR] The specified spatial interpolation option '
     write(LIS_logunit,*) '[ERR] is not supported for GDAS.'
     write(LIS_logunit,*) '[ERR] LIS is stopping.'
     call LIS_endrun()
  endif

!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GDAS & LDAS. For LDAS land 
! points not included in GDAS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gdas

