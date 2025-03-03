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
! !ROUTINE: read_gfs
! \label{read_gfs}
!
! !REVISION HISTORY:
!  16 Mar 2008: Sujay Kumar; Initial specification
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!
! !INTERFACE:
subroutine read_gfs( order, n, findex, name00, name03, name06, F06flag, ferror,try )

! !USES:  
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_timeMgrMod,    only : LDT_get_nstep, LDT_date2time
  use LDT_metforcingMod, only : LDT_forc
  use LDT_logMod,        only : LDT_logunit
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use gfs_forcingMod,    only : gfs_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: order    
  integer, intent(in)           :: n
  integer, intent(in)           :: findex
  character(len=*), intent(in)  :: name00
  character(len=*), intent(in)  :: name03
  character(len=*), intent(in)  :: name06
  logical, intent(in)           :: F06flag
  integer, intent(out)          :: ferror 
  integer, intent(inout)        :: try
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GFS forecast datasets, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain.
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
!    name of the 3 hour GFS forecast file
!  \item[name06]
!    name of the 6 hour GFS forecast file
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
!  \item[interp\_gfs](\ref{interp_gfs}) \newline
!    spatially interpolates a GFS variable
!  \end{description}
!EOP
!==== Local Variables=======================

  character(len=LDT_CONST_PATH_LEN) :: fname
  integer :: iv, c,r,t
  integer :: ferror1, ferror2, ferror3
  integer :: nforce
  integer :: ngfs
  real    :: glbdata_i(10,LDT_rc%ngrid(n))
  real    :: glbdata_a(10,LDT_rc%ngrid(n))
  real    :: glbdata_a_f06(10,LDT_rc%ngrid(n))
  integer :: nstep

!=== End Variable Definition =======================

  glbdata_a = LDT_rc%udef
  glbdata_a_f06 = LDT_rc%udef
  ngfs = (gfs_struc(n)%nc*gfs_struc(n)%nr)
!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
  nstep = LDT_get_nstep(LDT_rc,n)

  nforce = gfs_struc(n)%nmif

!--------------------------------------------------------------------------
! if there's a problem then ferror is set to zero
! read instantaneous fields
!--------------------------------------------------------------------------
  iv = 0

!--------------------------------------------------------------------------
! Set up to open file and retrieve specified field 
!--------------------------------------------------------------------------
  fname = name00
  call retrieve_gfs_variables(n, findex, fname,glbdata_i, ferror1)

!--------------------------------------------------------------------------
! read 3hr forecast for time averaged fields
!--------------------------------------------------------------------------
  fname = name03
  call retrieve_gfs_variables(n, findex, fname,glbdata_a, ferror2)

!--------------------------------------------------------------------------
! read 6hr forecast for time averaged fields, if required. 
!--------------------------------------------------------------------------
  if(F06flag) then 
     fname = name06
     call retrieve_gfs_variables(n, findex, fname,glbdata_a_f06, ferror3)
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
! Place the interpolated data into the LDT arrays
!--------------------------------------------------------------------------

  do iv=1,9
     do t=1,LDT_rc%ngrid(n)
        if(F06flag) then 
           if(iv.eq.3.or.iv.eq.4.or.iv.eq.8.or.iv.eq.9) then! these are time avgd fields
              if(order.eq.1) then 
                 LDT_forc(n,findex)%metdata1(iv,t) = 2*glbdata_a_F06(iv,t)-glbdata_a(iv,t)
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = 2*glbdata_a_F06(iv,t)-glbdata_a(iv,t)
              endif
           else ! these are instantaneous
              if(order.eq.1) then 
                 LDT_forc(n,findex)%metdata1(iv,t) = glbdata_i(iv,t)
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = glbdata_i(iv,t)
              endif
           endif
        else
           if(iv.eq.3.or.iv.eq.4.or.iv.eq.8.or.iv.eq.9) then! these are time avgd fields
              if(order.eq.1) then 
                 LDT_forc(n,findex)%metdata1(iv,t) = glbdata_a(iv,t)                         
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = glbdata_a(iv,t)                         
              endif
           else
              if(order.eq.1) then 
                 LDT_forc(n,findex)%metdata1(iv,t) = glbdata_i(iv,t)                         
              else
                 LDT_forc(n,findex)%metdata2(iv,t) = glbdata_i(iv,t)                         
              endif
           endif
        endif
     enddo
  enddo

  return

end subroutine read_gfs

!BOP
! 
! !ROUTINE: retrieve_gfs_variables
! \label{retrieve_gfs_variables}
! 
! !INTERFACE: 
subroutine retrieve_gfs_variables(n, findex, fname, glbdata, errorcode)
! !USES: 
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit,LDT_getNextUnitNumber,& 
       LDT_releaseUnitNumber, LDT_verify, LDT_warning
  use gfs_forcingMod,    only : gfs_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer               :: n 
  integer               :: findex
  character(len=*)      :: fname
  real                  :: glbdata(10,LDT_rc%ngrid(n))
  integer               :: errorcode
! 
! !DESCRIPTION: 
!   This subroutine retrieves GFS forcing variables, and interpolates
!  them to the LDT grid. 
! 
!EOP

  integer               :: ngfs
  real, allocatable :: f(:)
  real, dimension(LDT_rc%lnc(n), LDT_rc%lnr(n)) :: varfield
  integer :: igrib
  integer :: iv,c,r,t
  real    :: missingValue 
  integer :: iret
  integer :: ftn 
  integer, dimension(gfs_struc(n)%nmif) :: pds5, pds7, pds6,pds16
  integer :: pds5_val, pds7_val, pds16_val
  logical*1, allocatable :: lb(:)
  logical :: file_exists
  integer :: kk
  integer :: var_index
  integer :: nvars
  integer :: rc
  logical :: var_status(gfs_struc(n)%nmif)
  logical :: pcp_flag, var_found

#if(defined USE_GRIBAPI) 
  pds5 = (/ 011,051,204,205,033,034,001,059,214/) !parameter
  pds6 = -1
  pds7 = (/ 002,002,000,000,10,10,000,000,000 /) !htlev2
! index 10 indicates instantaneous, 003 indicates time average
  pds16 = (/010,010,003,003,010,010,010,003,003 /) 
  

  ngfs = (gfs_struc(n)%nc*gfs_struc(n)%nr)
  
  varfield = 0 
  errorcode = 1
  var_status = .false. 

  inquire (file=fname, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(fname),'r',iret)
     if(iret.ne.0) then 
        write(LDT_logunit,*) &
             'Could not open file: ',trim(fname)
        errorcode = 0
        return
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LDT_verify(iret, 'error in grib_count_in_file in read_gfs')

     allocate(lb(gfs_struc(n)%nc*gfs_struc(n)%nr))
     allocate(f(gfs_struc(n)%nc*gfs_struc(n)%nr))
     
     do kk=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LDT_warning(iret, 'error in grib_new_from_file in read_gfs')
        if(iret.ne.0) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LDT_verify(rc, 'error in grib_get: indicatorOfParameter in read_gfs')

        call grib_get(igrib,'level',pds7_val,rc)
        call LDT_verify(rc, 'error in grib_get: level in read_gfs')

        call grib_get(igrib,'timeRangeIndicator',pds16_val,rc)
        call LDT_verify(rc, 'error in grib_get: timeRangeIndicator in read_gfs')

        var_found = .false. 
        do iv=1,9
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
        call LDT_warning(rc, 'error in grib_get:values in read_gfs')

        if(rc.ne.0) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LDT_verify(rc, 'error in grib_get:missingValue in read_gfs')

        call grib_release(igrib,rc)
        call LDT_verify(rc, 'error in grib_release in read_gfs')
        
        if(var_found) then 
           lb = .false.
           do t=1,ngfs
              if(f(t).ne.missingValue) lb(t) = .true. 
           enddo
           
           pcp_flag = .false. 
           if(var_index.eq.8.or.var_index.eq.9) pcp_flag = .true. 
           
           call interp_gfs(n, findex,pcp_flag,ngfs,f,&
                lb,LDT_rc%gridDesc(n,:), &
                LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)
           
           do r=1,LDT_rc%lnr(n)
              do c=1,LDT_rc%lnc(n)
                 if(LDT_domain(n)%gindex(c,r).ne.-1) then 
                    glbdata(var_index,LDT_domain(n)%gindex(c,r)) =&
                         varfield(c,r)
                 endif
              enddo
           enddo
        endif

     enddo
     call grib_close_file(ftn)

     deallocate(lb)
     deallocate(f)     
         
     do kk=1,9
        if(.not.var_status(kk)) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           errorcode = 0
           
           return
        endif
     enddo
  else
     write(LDT_logunit,*) &
          'Could not find file: ',trim(fname)
     errorcode = 0
  endif
#endif     
end subroutine retrieve_gfs_variables


!BOP
! !ROUTINE: interp_gfs
! \label{interp_gfs}
!
! !INTERFACE:
subroutine interp_gfs(n, findex, pcp_flag,ngfs,f,lb,lis_gds,nc,nr, &
     varfield)
! !USES:
  use LDT_coreMod,     only : LDT_rc, LDT_domain
  use gfs_forcingMod,  only : gfs_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  logical, intent(in) :: pcp_flag
  integer, intent(in) :: ngfs
  real, intent(out)   :: f(ngfs)
  logical*1           :: lb(ngfs)
  real                :: lis_gds(20)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given GFS field 
!   to the LDT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[pcp\_flag]
!  flag indicating precipitation variables
! \item[ngfs]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LDT grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LDT grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LDT grid
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
  integer :: count1,i,j,mo

  real, dimension(nc*nr) :: lis1d

  logical*1 :: lo(nc*nr)

!=== End variable declarations
  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LDT grid
!-----------------------------------------------------------------------  
  if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 
     call bilinear_interp(lis_gds,lb,f,lo,lis1d,gfs_struc(n)%mi,mo,&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          gfs_struc(n)%w111,gfs_struc(n)%w121,&
          gfs_struc(n)%w211,gfs_struc(n)%w221,&
          gfs_struc(n)%n111,gfs_struc(n)%n121,&
          gfs_struc(n)%n211,gfs_struc(n)%n221,LDT_rc%udef, iret)
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
     if (pcp_flag)then     
        call conserv_interp(lis_gds,lb,f,lo,lis1d,gfs_struc(n)%mi,mo, & 
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             gfs_struc(n)%w112,gfs_struc(n)%w122,&
             gfs_struc(n)%w212,gfs_struc(n)%w222,&
             gfs_struc(n)%n112,gfs_struc(n)%n122,&
             gfs_struc(n)%n212,gfs_struc(n)%n222,LDT_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,lb,f,lo,lis1d,gfs_struc(n)%mi,mo,&
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             gfs_struc(n)%w111,gfs_struc(n)%w121,&
             gfs_struc(n)%w211,gfs_struc(n)%w221,&
             gfs_struc(n)%n111,gfs_struc(n)%n121,&
             gfs_struc(n)%n211,gfs_struc(n)%n221,LDT_rc%udef,iret)
     endif
  endif

!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GFS & LDAS. For LDAS land 
! points not included in GFS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gfs
