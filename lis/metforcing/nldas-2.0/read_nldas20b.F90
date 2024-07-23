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
! !ROUTINE: read_nldas20b
!  \label{read_nldas20b}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from read_nldas2b.F90)
!
! !INTERFACE:
subroutine read_nldas20b(n,kk,findex,order,name,ferror)
! !USES:
  use LIS_coreMod, only        : LIS_rc,LIS_domain
  use LIS_logMod, only         : LIS_logunit,LIS_verify,LIS_warning
  use LIS_metforcingMod, only  : LIS_forc
  use nldas20_forcingMod, only : nldas20_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)          :: n
  integer, intent(in)          :: kk     ! Forecast member index
  integer, intent(in)          :: findex ! Forcing index
  integer, intent(in)          :: order
  character(len=*), intent(in) :: name
  integer, intent(out)         :: ferror
!
! !DESCRIPTION:
!  For the given time, reads values from NLDAS-2 netCDF-4 FORB data,
!  transforms into 10 variables, and interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[kk]
!    forecast member index
!  \item[findex]
!    forcing index
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    hourly instance, order=2, read the next hourly instance)
!  \item[name]
!    name of the hourly NLDAS-2 forecast file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[interp\_nldas20](\ref{interp_nldas20}) \newline
!    spatially interpolates a NLDAS-2 variable
!  \end{description}
!
!EOP
  integer                :: iv,ftn
  integer                :: nldas20,paramid
  integer                :: k,t,c,r,iret,rc
  real, parameter        :: missingValue = -9999.0
  integer, parameter     :: nvars = 10
  logical                :: pcp_flag
  logical                :: file_exists
  character(len=20)      :: input_varname(nvars)
  logical*1, allocatable :: lb(:)
  real, allocatable      :: f(:)
  real, allocatable      :: nldas20_forcing(:,:)
  real                   :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real    :: dummy(nldas20_struc(n)%ncold,nldas20_struc(n)%nrold)

  ferror = 1
  iv = 0

  input_varname(1)  = "SWdown"
  input_varname(2)  = "Rainf"
  input_varname(3)  = "CRainf"
  input_varname(4)  = "ACond"
  input_varname(5)  = "Tair"
  input_varname(6)  = "Qair"
  input_varname(7)  = "PSurf"
  input_varname(8)  = "Wind_E"
  input_varname(9)  = "Wind_N"
  input_varname(10) = "PhiS"

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  nldas20 = (nldas20_struc(n)%ncold*nldas20_struc(n)%nrold)

  allocate(nldas20_forcing(nldas20_struc(n)%ncold*nldas20_struc(n)%nrold,nvars))

  varfield = 0
  ferror = 1

  inquire(file=name,exist=file_exists)
  if (file_exists) then
     iret = nf90_open(path=name,mode=NF90_NOWRITE,ncid=ftn)
     if (iret.ne.0) then
        write(LIS_logunit,*) "[WARN] Could not open file: ",trim(name)
        ferror = 0
        return
     endif

     allocate(lb(nldas20_struc(n)%ncold*nldas20_struc(n)%nrold))
     allocate(f(nldas20_struc(n)%ncold*nldas20_struc(n)%nrold))

     do k = 1,nvars
        iret = nf90_inq_varid(ftn,trim(input_varname(k)),paramid)
        call LIS_verify(iret,trim(input_varname(k))//              &
             " field not found in the hourly file")
        iret = nf90_get_var(ftn,paramid,dummy)
        call LIS_verify(iret,"Error in nf90_get_var")
        ! write(LIS_logunit,*) "[INFO] Read field: ",input_varname(k)

        f = LIS_rc%udef     ! Initialize forcing
        t = 0
        do r = 1,nldas20_struc(n)%nrold
           do c = 1,nldas20_struc(n)%ncold
              t = t + 1
              f(t) = dummy(c,r)
           enddo
        enddo

        lb = .false.
        do t = 1,nldas20
           if (f(t).ne.missingValue) then
              nldas20_forcing(t,k) = f(t)
              lb(t) = .true.
           else
              nldas20_forcing(t,k) = LIS_rc%udef
           endif
        enddo
     enddo
     deallocate(f)

     iret = nf90_close(ftn)

     do iv = 1,nvars
        pcp_flag = .false.
        if ((iv.eq.2).or.(iv.eq.3)) pcp_flag = .true.

        call interp_nldas20(n,findex,LIS_rc%mo,pcp_flag,nldas20,   &
             nldas20_forcing(:,iv),                                &
             lb,LIS_rc%gridDesc(n,:),                              &
             LIS_rc%lnc(n),LIS_rc%lnr(n),varfield)

        do r = 1,LIS_rc%lnr(n)
           do c = 1,LIS_rc%lnc(n)
              if (LIS_domain(n)%gindex(c,r).ne.-1) then

! MODEL LEVEL TAIR CASE
                 if (iv.eq.5) then
                    if (nldas20_struc(n)%model_level_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,1,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,1,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL LEVEL SPFH CASE
                 if (iv.eq.6) then
                    if (nldas20_struc(n)%model_level_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,2,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,2,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! USE MODEL BASED DSWRF CASE
                 if (iv.eq.1) then
                    if (nldas20_struc(n)%model_dswrf_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,3,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,3,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL LEVEL USE AERODYNAMIC CONDUCTANCE
                 if (iv.eq.4) then
                    if (nldas20_struc(n)%model_level_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,13,         &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,13,         &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL LEVEL UWIND CASE
                 if (iv.eq.8) then
                    if (nldas20_struc(n)%model_level_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,5,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,5,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL LEVEL VWIND CASE
                 if (iv.eq.9) then
                    if (nldas20_struc(n)%model_level_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,6,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,6,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL LEVEL PRESSURE CASE
                 if (iv.eq.7) then
                    if (nldas20_struc(n)%model_level_press.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,7,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,7,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL BASED PRECIP CASE
                 if (iv.eq.2) then
                    if (nldas20_struc(n)%model_pcp_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,8,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,8,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL BASED CONV. PRECIP CASE
                 if (iv.eq.3) then
                    if (nldas20_struc(n)%model_pcp_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,9,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,9,          &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

! MODEL FORCING HEIGHT CASE
                 if (iv.eq.10) then
                    if (nldas20_struc(n)%model_level_data.gt.0) then
                       if (order.eq.1) then
                          nldas20_struc(n)%metdata1(kk,12,         &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       elseif (order.eq.2) then
                          nldas20_struc(n)%metdata2(kk,12,         &
                               LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                       endif
                    endif
                 endif

              endif
           enddo
        enddo
     enddo
     deallocate(lb)
  else
     write(LIS_logunit,*) "[WARN] Could not find file: ",trim(name)
     ferror = 0
  endif

  deallocate(nldas20_forcing)
#endif

end subroutine read_nldas20b

