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
! !ROUTINE: read_nldas20a
!  \label{read_nldas20a}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from read_nldas2a.F90)
!
! !INTERFACE:
subroutine read_nldas20a(n,kk,findex,order,month,name,ferror)
! !USES:
  use LIS_coreMod
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
  integer, intent(out)         :: month
  character(len=*), intent(in) :: name
  integer, intent(out)         :: ferror
!
! !DESCRIPTION:
!  For the given time, reads values from NLDAS-2 netCDF-4 FORA data,
!  transforms into 11 variables, and interpolates to the LIS domain.
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
!  \item[month]
!    current month
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
  integer, parameter     :: nvars = 11
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

  input_varname(1)  = "Tair"
  input_varname(2)  = "Qair"
  input_varname(3)  = "SWdown"
  input_varname(4)  = "LWdown"
  input_varname(5)  = "Wind_E"
  input_varname(6)  = "Wind_N"
  input_varname(7)  = "PSurf"
  input_varname(8)  = "Rainf"
  input_varname(9)  = "CRainf_frac"
  input_varname(10) = "PotEvap"
  input_varname(11) = "CAPE"

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

     if (LIS_rc%met_ecor(findex).ne."none") then
        do t = 1,nldas20
           if (lb(t)) then
              call nldas20_ec_removal(n,t,nldas20_forcing(t,1),    &
                   nldas20_forcing(t,2),nldas20_forcing(t,4),      &
                   nldas20_forcing(t,7))
           endif
        enddo
     endif

     do iv = 1,nvars
        pcp_flag = .false.
        if ((iv.eq.8).or.(iv.eq.9)) pcp_flag = .true.

        call interp_nldas20(n,findex,month,pcp_flag,nldas20,       &
             nldas20_forcing(:,iv),                                &
             lb,LIS_rc%gridDesc(n,:),                              &
             LIS_rc%lnc(n),LIS_rc%lnr(n),varfield)

        do r = 1,LIS_rc%lnr(n)
           do c = 1,LIS_rc%lnc(n)
              if (LIS_domain(n)%gindex(c,r).ne.-1) then
                 if (order.eq.1) then
                    nldas20_struc(n)%metdata1(kk,iv,               &
                         LIS_domain(n)%gindex(c,r)) = varfield(c,r)
                 elseif (order.eq.2) then
                    nldas20_struc(n)%metdata2(kk,iv,               &
                         LIS_domain(n)%gindex(c,r)) = varfield(c,r)
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

end subroutine read_nldas20a

!BOP
! !ROUTINE: interp_nldas20
! \label{interp_nldas20}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from read_nldas2a.F90)
!
! !INTERFACE:
subroutine interp_nldas20(n,findex,month,pcp_flag,input_size,    &
     input_data,input_bitmap,lis_gds,nc,nr,output_2d)
! !USES:
  use LIS_coreMod, only        : LIS_rc,LIS_domain
  use nldas20_forcingMod, only : nldas20_struc
  use LIS_spatialDownscalingMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)   :: n
  integer, intent(in)   :: findex
  integer, intent(in)   :: month
  logical, intent(in)   :: pcp_flag
  integer, intent(in)   :: input_size
  real, intent(in)      :: input_data(input_size)
  logical*1, intent(in) :: input_bitmap(input_size)
  real, intent(in)      :: lis_gds(50)
  integer, intent(in)   :: nc
  integer, intent(in)   :: nr
  real, intent(inout)   :: output_2d(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given NLDAS-2 field to the LIS grid.
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing index
!  \item[month]
!    current month
!  \item[pcp_flag]
!    flag indicating if it is a precip field
!  \item[input\_size]
!    number of elements in the input grid
!  \item[input\_bitmap]
!    input bitmap
!  \item[lis\_gds]
!    array description of the LIS grid
!  \item[nc]
!    number of columns (in the east-west dimension) in the LIS grid
!  \item[nr]
!    number of rows (in the north-south dimension) in the LIS grid
!  \item[output\_2d]
!    output interpolated field
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using nearest neighbor interpolation
! \end{description}
!
!EOP
  integer :: iret
  integer :: mo
  integer :: count1,i,j
  real, dimension(nc*nr) :: output_data
  logical*1 :: output_bitmap(nc*nr)

!=== End variable declarations

  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap.
!-----------------------------------------------------------------------
  output_bitmap = .true.

!-----------------------------------------------------------------------
! Apply downscaling
!-----------------------------------------------------------------------

!Bailing Li
!  if LIS_rc%pcp_downscale(findex).eq.1: spatial downscaling and scaling factors are calculated in LIS
!  if LIS_rc%pcp_downscale(findex).eq.2: spatial downscaling and bias-correction using ratios of two
!     PCP climatologies which are calculated in LDT and stored in lis input file

  if (pcp_flag) then
     if (LIS_rc%pcp_downscale(findex).eq.1) then
!input_data becomes the ratio field.
        call LIS_generatePcpClimoRatioField(n,findex,"NLDAS2",     &
             month,input_size,input_data,input_bitmap)
     elseif (pcp_flag.and.(LIS_rc%pcp_downscale(findex).eq.2)) then
        call LIS_readPcpClimoRatioField(n,findex,month)
     endif
  endif

!-----------------------------------------------------------------------
! Interpolate to LIS grid
!-----------------------------------------------------------------------
  select case(LIS_rc%met_interp(findex))

  case("bilinear")
     call bilinear_interp(lis_gds,input_bitmap,input_data,      &
          output_bitmap,output_data,nldas20_struc(n)%mi,mo,     &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                  &
          nldas20_struc(n)%w111,nldas20_struc(n)%w121,          &
          nldas20_struc(n)%w211,nldas20_struc(n)%w221,          &
          nldas20_struc(n)%n111,nldas20_struc(n)%n121,          &
          nldas20_struc(n)%n211,nldas20_struc(n)%n221,          &
          LIS_rc%udef,iret)

  case("budget-bilinear")
     if (pcp_flag) then
        call conserv_interp(lis_gds,input_bitmap,input_data,    &
             output_bitmap,output_data,nldas20_struc(n)%mi,mo,  &
             LIS_domain(n)%lat,LIS_domain(n)%lon,               &
             nldas20_struc(n)%w112,nldas20_struc(n)%w122,       &
             nldas20_struc(n)%w212,nldas20_struc(n)%w222,       &
             nldas20_struc(n)%n112,nldas20_struc(n)%n122,       &
             nldas20_struc(n)%n212,nldas20_struc(n)%n222,       &
             LIS_rc%udef,iret)

     else
        call bilinear_interp(lis_gds,input_bitmap,input_data,   &
             output_bitmap,output_data,nldas20_struc(n)%mi,mo,  &
             LIS_domain(n)%lat,LIS_domain(n)%lon,               &
             nldas20_struc(n)%w111,nldas20_struc(n)%w121,       &
             nldas20_struc(n)%w211,nldas20_struc(n)%w221,       &
             nldas20_struc(n)%n111,nldas20_struc(n)%n121,       &
             nldas20_struc(n)%n211,nldas20_struc(n)%n221,       &
             LIS_rc%udef,iret)
     endif

  case("neighbor")
     call neighbor_interp(lis_gds,input_bitmap,input_data,      &
          output_bitmap,output_data,nldas20_struc(n)%mi,mo,     &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                  &
          nldas20_struc(n)%n113,LIS_rc%udef,iret)

  end select

  if (pcp_flag.and.(LIS_rc%pcp_downscale(findex).ne.0)) then
     call LIS_pcpClimoDownscaling(n,findex,month,nc*nr,         &
          output_data,output_bitmap)
  endif

!-----------------------------------------------------------------------
! convert the interpolated data to 2d.
!-----------------------------------------------------------------------
  count1 = 0
  do j = 1,nr
     do i = 1,nc
        output_2d(i,j) = output_data(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_nldas20

