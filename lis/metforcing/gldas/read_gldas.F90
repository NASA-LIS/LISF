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
! !ROUTINE: read_gldas
! \label{read_gldas}
!
! !REVISION HISTORY:
!  19 Sept 2008: Sujay Kumar: Initial Implementation
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!
! !INTERFACE:
subroutine read_gldas( order, n, findex, name, ferror, try )
! !USES:  
  use LIS_coreMod,         only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,      only : LIS_get_nstep
  use LIS_metforcingMod,   only : LIS_forc
  use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_warning
  use gldas_forcingMod,    only : gldas_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: order    
  integer, intent(in)           :: n
  integer, intent(in)           :: findex
  character(len=*),  intent(in) :: name
  integer, intent(out)          :: ferror 
  integer, intent(inout)        :: try
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GLDAS data and spatially interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 3 hour GLDAS forecast file
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \item[try]
!    index of the tries (in case of missing data)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gldas](\ref{interp_gldas}) \newline
!    spatially interpolates a GLDAS variable
!  \item[fillgaps\_gldas](\ref{fillgaps_gldas}) \newline
!    fills the data gaps due to mismatches in landmask of LIS
!    domain and GLDAS mask.
!  \end{description}
!EOP
!==== Local Variables=======================
  integer                    :: ftn
  integer                    :: c,r,t
  integer                    :: nforce, ngldas
  logical                    :: file_exists
  integer                    :: var_index
  integer                    :: iv, iv_total
  logical                    :: var_found,pcp_flag
  logical                    :: var_status(10)
  logical*1, allocatable         :: lb(:)
  real,  allocatable             :: f(:)
  real                       :: missingValue
  real                       :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                    :: kk,nvars
  integer                    :: igrib
  integer                    :: pds7_1(10), pds6_1(10)
  integer                    :: pds5(10), pds7(10), pds6(10)
  integer                    :: pds5_val, pds7_val, pds6_val
  integer                    :: jpds5, jpds7, jpds6
  integer                    :: rc,status,iret
!=== End Variable Definition =======================

#if(defined USE_GRIBAPI) 
!--------------------------------------------------------------------------
! Set the GRIB parameter specifiers
!--------------------------------------------------------------------------
  pds5 = (/ 011,051,204,205,033,034,001,132,131,63 /) !parameter
  pds6 = (/ 001,001,001,001,001,001,001,001,001,001 /)
  pds7 = (/ 000,000,000,000,000,000,000,000,000,000 /) !htlev2

  pds6_1 = (/ 001,001,001,001,105,105,001,001,001,001 /)
  pds7_1 = (/ 000,000,000,000,010,010,000,000,000,000 /) !htlev2

  ngldas = (gldas_struc(n)%ncold*gldas_struc(n)%nrold)
  nforce = LIS_rc%met_nf(findex)
  ferror = 1  

  iv_total = 9
  inquire (file=name, exist=file_exists)
  if (file_exists) then      

     call grib_open_file(ftn,trim(name),'r',iret)
     if(iret.ne.0) then 
        write(LIS_logunit,*) &
             'Could not open file: ',trim(name)
        ferror = 0
        return
     endif

     call grib_count_in_file(ftn,nvars,iret)
     call LIS_verify(iret, 'error in grib_count_in_file in read_gldas')

     allocate(lb(gldas_struc(n)%ncold*gldas_struc(n)%nrold))
     allocate(f(gldas_struc(n)%ncold*gldas_struc(n)%nrold))
     
     do kk=1,nvars
        call grib_new_from_file(ftn, igrib, iret)
        call LIS_warning(iret, 'error in grib_new_from_file in read_gldas')
        if(iret.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve message in file: ',trim(name)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_gldas')

        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in read_gldas')

        call grib_get(igrib,'indicatorOfTypeOfLevel',pds6_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfTypeOfLevel in read_gldas')
        
        var_found = .false. 
        do iv=1,iv_total
           jpds5  = pds5(iv)
           jpds7  = pds7(iv)
           jpds6 = pds6(iv)
           if(LIS_rc%yr.ge.2008.and.LIS_rc%mo.gt.4) then 
              jpds6 = pds6_1(iv)
              jpds7 = pds7_1(iv)
           endif
   !change back on 2008/06        
           if(LIS_rc%yr.ge.2008.and.LIS_rc%mo.gt.5) then 
              jpds6 = pds6(iv)
              jpds7 = pds7(iv)
           endif

           if((pds5_val.eq.jpds5).and.&
                (pds7_val.eq.jpds7).and.&
                (pds6_val.eq.jpds6)) then
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo

        f = -9999.0
        call grib_get(igrib,'values',f,rc)
        call LIS_warning(rc, 'error in grib_get:values in read_gldas')

        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve values in file: ',trim(name)
           ferror = 0
           deallocate(lb)
           deallocate(f)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in read_gldas')

        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in read_gldas')
        
        if (var_found) then 
           lb = .false.
           do t=1,ngldas
              if(f(t).ne.missingValue) lb(t) = .true. 
           enddo
           where ( f == missingValue )
              f = LIS_rc%udef
           endwhere
           
           pcp_flag = .false. 
           if(var_index.eq.8.or.var_index.eq.9) pcp_flag = .true. 
           
           call interp_gldas(n, findex, pcp_flag,ngldas,&
                f,lb,LIS_rc%gridDesc(n,:), &
                LIS_rc%lnc(n),LIS_rc%lnr(n),varfield)

           call fillgaps_gldas(n,1,varfield)

           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gldas_struc(n)%metdata1(var_index,&
                            LIS_domain(n)%gindex(c,r)) =&
                            varfield(c,r)
                    else
                       gldas_struc(n)%metdata2(var_index,&
                            LIS_domain(n)%gindex(c,r)) = &
                            varfield(c,r)
                    endif
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
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(name)
           ferror = 0
           
           return
        endif
     enddo
  else
     write(LIS_logunit,*) &
          'Could not find file: ',trim(name)
     ferror = 0
  endif

#endif


end subroutine read_gldas


!BOP
! !ROUTINE: interp_gldas
! \label{interp_gldas}
!
! !INTERFACE:
subroutine interp_gldas(n,findex,pcp_flag,ngldas,f,lb,lis_gds,nc,nr, &
                        varfield)
! !USES:
  use LIS_coreMod,      only : LIS_rc, LIS_domain
  use gldas_forcingMod, only : gldas_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  logical, intent(in) :: pcp_flag
  integer, intent(in) :: ngldas
  real, intent(out)   :: f(ngldas)
  logical*1           :: lb(ngldas)
  real                :: lis_gds(50)
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real, intent(out)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine spatially interpolates a given GLDAS field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[findex]
!  index of the forcing source
! \item[pcp\_flag]
!  flag indicating if precip variables are being interpolated
! \item[ngldas]
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
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
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
! Interpolate to LIS grid
!-----------------------------------------------------------------------  
  if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
     call bilinear_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gldas_struc(n)%w111,gldas_struc(n)%w121,&
          gldas_struc(n)%w211,gldas_struc(n)%w221,&
          gldas_struc(n)%n111,gldas_struc(n)%n121,&
          gldas_struc(n)%n211,gldas_struc(n)%n221,LIS_rc%udef, iret)
  elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
     if (pcp_flag)then     
        call conserv_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo, & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gldas_struc(n)%w112,gldas_struc(n)%w122,&
             gldas_struc(n)%w212,gldas_struc(n)%w222,&
             gldas_struc(n)%n112,gldas_struc(n)%n122,&
             gldas_struc(n)%n212,gldas_struc(n)%n222,LIS_rc%udef,iret)
     else 
        call bilinear_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gldas_struc(n)%w111,gldas_struc(n)%w121,&
             gldas_struc(n)%w211,gldas_struc(n)%w221,&
             gldas_struc(n)%n111,gldas_struc(n)%n121,&
             gldas_struc(n)%n211,gldas_struc(n)%n221,LIS_rc%udef,iret)
     endif
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
     call neighbor_interp(lis_gds,lb,f,lo,lis1d,gldas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gldas_struc(n)%n113,LIS_rc%udef,iret)
  endif
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GLDAS & LDAS. For LDAS land 
! points not included in GLDAS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gldas
