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
! !ROUTINE: read_merra2bias
! \label{read_merra2bias}
!
! !REVISION HISTORY:
! 29 Jun 2026: Kristen Whitney, initial code (based on geos-itbias)
!
! !INTERFACE:
subroutine read_merra2bias(n,order,month,findex,                     &
     merra2name,                                                       &
     lapseratefname,                                                 &
     merra2biasforc,ferror)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc
  use merra2bias_forcingMod, only : merra2bias_struc
  use LIS_constantsMod, only : LIS_CONST_LAPSE_RATE
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)          :: n
  integer, intent(in)          :: order
  integer, intent(in)          :: month
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: merra2name
  character(len=*), intent(in) :: lapseratefname
  real, intent(inout)          :: &
       merra2biasforc(merra2bias_struc(n)%nvars,&
       LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer, intent(out)         :: ferror

! !DESCRIPTION:
!  For the given time, reads parameters from MERRA2bias data, transforms
!  into 4 LIS forcing parameters and interpolates to the LIS domain. \newline
!
! MERRA2bias FORCING VARIABLES: \newline
!  1. tair    Temperature of air at 2-m, 10-m, or LML [$K$] \newline
!  2. qair    Specific humidity of air at 2-m, 10-m, or LML [$kg/kg$] \newline
!  3. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  4. lwgab   Downward longwave radiation at the ground [$W/m^2$] \newline
!
!  The arguments are:
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 1 hour MERRA2bias analysis file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP
  integer   :: ftn_merra2
  integer   :: tmpId,qId,psId
  integer   :: lwgabId
  integer   :: nr_index,nc_index
  logical   :: file_exists
  integer   :: mo
  real      :: tair(merra2bias_struc(n)%ncold,merra2bias_struc(n)%nrold)
  real      :: qair(merra2bias_struc(n)%ncold,merra2bias_struc(n)%nrold)
  real      :: ps(merra2bias_struc(n)%ncold,merra2bias_struc(n)%nrold)
  real      :: lwgab(merra2bias_struc(n)%ncold,merra2bias_struc(n)%nrold)
  integer                           :: ftn_drate
  integer                           :: lapserateid
  real, allocatable                 :: lapse_rate_in(:,:)
  integer :: nc_sub, nr_sub
  real, allocatable :: lapse_rate_sub(:,:)
  real, allocatable :: lat_sub(:), lon_sub(:)
  real, allocatable :: lat_full(:), lon_full(:)
  integer :: i0, j0
  integer :: latid, lonid

  external :: interp_merra2bias_var
  external :: interp_lapserate_merra2bias
! __________________________________________________________________________

#if (defined USE_NETCDF3)
  write(LIS_logunit,*) "[ERR] MERRA2bias reader requires NetCDF4"
  call LIS_endrun()
#endif

#if (defined USE_NETCDF4)
  ferror = 0
  nr_index = merra2bias_struc(n)%nrold
  nc_index = merra2bias_struc(n)%ncold
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

! Read required fields from MERRA2bias input file:
  inquire(file=merra2name,exist=file_exists)
  if(file_exists) then
     write(LIS_logunit,*)                                           &
          '[INFO] Reading MERRA2bias file: ',trim(merra2name)
     call LIS_verify(nf90_open(path=trim(merra2name),                 &
          mode=NF90_NOWRITE,ncid=ftn_merra2),                         &
          'nf90_open failed for merra2file in read_merra2bias')

     call LIS_verify(nf90_inq_varid(ftn_merra2,'PS',psId),            &
          'nf90_inq_varid failed for ps in read_merra2bias')
     call LIS_verify(nf90_get_var(ftn_merra2,psId,ps),                &
          'nf90_get_var failed for ps in read_merra2bias')

     call LIS_verify(nf90_inq_varid(ftn_merra2,'LWGAB',lwgabId),      &
          'nf90_inq_varid failed for LWGAB in read_merra2bias')
     call LIS_verify(nf90_get_var(ftn_merra2,lwgabId,lwgab),          &
          'nf90_get_var failed for LWGAB in read_merra2bias')

     call LIS_verify(nf90_inq_varid(ftn_merra2,'QV2M',qId),           &
          'nf90_inq_varid failed for QV2M in read_merra2bias')
     call LIS_verify(nf90_get_var(ftn_merra2,qId,qair),               &
          'nf90_get_var failed for QV2M in read_merra2bias')

     call LIS_verify(nf90_inq_varid(ftn_merra2,'T2M',tmpId),          &
          'nf90_inq_varid failed for T2M in read_merra2bias')
     call LIS_verify(nf90_get_var(ftn_merra2,tmpId,Tair),             &
          'nf90_get_var failed for T2M in read_merra2bias')

     call LIS_verify(nf90_close(ftn_merra2),                          &
          '[WARN] Failed to close slvfile in read_merra2bias.')

     call interp_merra2bias_var(n,findex,month,tair,  1,.false.,    &
          merra2biasforc)
     call interp_merra2bias_var(n,findex,month,qair,  2,.false.,    &
          merra2biasforc)
     call interp_merra2bias_var(n,findex,month,ps,    3,.false.,    &
          merra2biasforc)
     call interp_merra2bias_var(n,findex,month,lwgab, 4,.false.,     &
          merra2biasforc)

    ! Obtain lapse rates
     if ((merra2bias_struc(n)%usedynlapserate.eq.1).and.            &
          ((LIS_rc%met_ecor(findex).eq."lapse-rate").or.               &
          (LIS_rc%met_ecor(findex).eq."lapse-rate and slope-aspect").or.&
          (LIS_rc%met_ecor(findex).eq."micromet"))) then
        inquire(file=trim(lapseratefname),exist=file_exists)
        if (file_exists) then
           write(LIS_logunit,*) 'Reading lapse rate file: ', &
                trim(lapseratefname)

           ! open file and get ID for lapse rate
           call LIS_verify(nf90_open(path=trim(lapseratefname),     &
                mode=NF90_NOWRITE,ncid=ftn_drate),                  &
                'nf90_open failed for '//trim(lapseratefname))
           call LIS_verify(nf90_inq_varid(ftn_drate,'lapse_rate',   &
                lapserateid),                                       &
                'nf90_inq_varid failed for lapse_rate')

           ! Get subdomain dimensions directly from dimids 1=lat, 2=lon
           call LIS_verify(nf90_inquire_dimension(ftn_drate,1, &
                len=nr_sub), 'inq lat dim failed (sub)')
           call LIS_verify(nf90_inquire_dimension(ftn_drate,2, &
                len=nc_sub), 'inq lon dim failed (sub)')

           ! allocate full-size lapse rate
           allocate(lapse_rate_in(merra2bias_struc(n)%ncold,   &
                merra2bias_struc(n)%nrold))
           lapse_rate_in = 1.e+15
           ! Case 1: full global lapse-rate file
           if (nc_sub == merra2bias_struc(n)%ncold .and.       &
                nr_sub == merra2bias_struc(n)%nrold) then
              call LIS_verify(nf90_get_var(ftn_drate,lapserateid, &
                   lapse_rate_in), 'get lapse_rate failed')

              ! Case 2: subdomain lapse-rate file
           else if (nc_sub <= merra2bias_struc(n)%ncold .and.  &
                nr_sub <= merra2bias_struc(n)%nrold) then
              allocate(lapse_rate_sub(nc_sub,nr_sub))
              allocate(lat_sub(nr_sub), lon_sub(nc_sub))
              allocate(lat_full(merra2bias_struc(n)%nrold), &
                   lon_full(merra2bias_struc(n)%ncold))

              ! Read subdomain values
              call LIS_verify(nf90_get_var(ftn_drate,lapserateid,    &
                   lapse_rate_sub), 'get lapse_rate_sub failed')
              call LIS_verify(nf90_inq_varid(ftn_drate,'lat',latid), &
                   'inq lat varid failed')
              call LIS_verify(nf90_inq_varid(ftn_drate,'lon',lonid), &
                   'inq lon varid failed')
              call LIS_verify(nf90_get_var(ftn_drate,latid,lat_sub), &
                   'get lat_sub failed')
              call LIS_verify(nf90_get_var(ftn_drate,lonid,lon_sub), &
                   'get lon_sub failed')

              ! Read full-domain coordinates from slv file
              call LIS_verify(nf90_open(path=trim(merra2name),         &
                   mode=NF90_NOWRITE,ncid=ftn_merra2),                 &
                   'nf90_open failed for merra2file in read_merra2bias')
              call LIS_verify(nf90_inq_varid(ftn_merra2,'lat',latid),  &
                   'inq lat varid failed (full)')
              call LIS_verify(nf90_inq_varid(ftn_merra2,'lon',lonid),  &
                   'inq lon varid failed (full)')
              call LIS_verify(nf90_get_var(ftn_merra2,latid,lat_full), &
                   'get lat_full failed')
              call LIS_verify(nf90_get_var(ftn_merra2,lonid,lon_full), &
                   'get lon_full failed')
              call LIS_verify(nf90_close(ftn_merra2), &
                   '[WARN] Failed to close slvfile in read_merra2bias.')

              ! Find paste offsets
              call find_index(lon_full, lon_sub(1), i0)
              call find_index(lat_full, lat_sub(1), j0)

              ! Paste subdomain into full-domain buffer
              lapse_rate_in(i0:i0+nc_sub-1, j0:j0+nr_sub-1) = &
                   lapse_rate_sub(:,:)

              ! Cleanup
              deallocate(lapse_rate_sub, lat_sub, lon_sub, lat_full, &
                   lon_full)

              ! Case 3: unexpected
           else
              write(LIS_logunit,*) '[ERR] Unexpected lapse-rate dims: ', &
                   nc_sub, ' x ', nr_sub
              call LIS_endrun()
           end if

           call LIS_verify(nf90_close(ftn_drate), &
                'close lapse-rate file failed')
           call interp_lapserate_merra2bias(n,order,lapse_rate_in)

           deallocate(lapse_rate_in)

        else
           write(LIS_logunit,*) '[WARN] Could not find ', &
                trim(lapseratefname)
           write(LIS_logunit,*) '[WARN] Using static lapse rate.'
           merra2bias_struc(n)%lapserate1 = LIS_CONST_LAPSE_RATE
           merra2bias_struc(n)%lapserate2 = LIS_CONST_LAPSE_RATE
        endif
     endif
  else
     call LIS_endrun()
  endif

#endif
contains
  subroutine find_index(x,val,idx)
    real, intent(in) :: x(:), val
    integer, intent(out) :: idx
    real :: dmin
    integer :: k
    dmin = huge(1.0); idx = -1
    do k=1,size(x)
       if (abs(x(k)-val) < dmin) then
          dmin = abs(x(k)-val)
          idx = k
       end if
    enddo
  end subroutine find_index
end subroutine read_merra2bias

!BOP
!
! !ROUTINE: interp_merra2bias_var
! \label{interp_merra2bias_var}
!
! !INTERFACE:
subroutine interp_merra2bias_var(n,findex,month,input_var,var_index, &
     pcp_flag,merra2biasforc)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use merra2bias_forcingMod, only : merra2bias_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(merra2bias_struc(n)%ncold,       &
       merra2bias_struc(n)%nrold)
  integer, intent(in)    :: var_index
  logical, intent(in)    :: pcp_flag
  real,    intent(inout) :: merra2biasforc(merra2bias_struc(n)%nvars,  &
       LIS_rc%lnc(n)*LIS_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine spatially interpolates a MERRA2bias field
!  to the LIS running domain
!
!EOP
  integer   :: c,r,k,iret
  real      :: f (merra2bias_struc(n)%ncold*merra2bias_struc(n)%nrold)
  logical*1 :: lb(merra2bias_struc(n)%ncold*merra2bias_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size

  external :: conserv_interp
  external :: bilinear_interp
  external :: neighbor_interp

! _____________________________________________________________

  input_size = merra2bias_struc(n)%ncold*merra2bias_struc(n)%nrold

!-----------------------------------------------------------------------
! Apply downscaling
!-----------------------------------------------------------------------
  lb = .true.
  do r = 1,merra2bias_struc(n)%nrold
     do c = 1,merra2bias_struc(n)%ncold
        k = c+(r-1)*merra2bias_struc(n)%ncold
        f(k) = input_var(c,r)
        if (f(k).eq.1.e+15) then
           f(k)  = LIS_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo

  if (pcp_flag.and.(LIS_rc%pcp_downscale(findex).ne.0)) then
     ! input_data becomes the ratio field.
     call LIS_generatePcpClimoRatioField(n,findex,"MERRA2bias",       &
          month,input_size,f,lb)
  endif

  if (pcp_flag.and.                                                    &
       (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")) then

     call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,                 &
          merra2biasforc(var_index,:),                                 &
          merra2bias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),          &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          merra2bias_struc(n)%w112,merra2bias_struc(n)%w122,           &
          merra2bias_struc(n)%w212,merra2bias_struc(n)%w222,           &
          merra2bias_struc(n)%n112,merra2bias_struc(n)%n122,           &
          merra2bias_struc(n)%n212,merra2bias_struc(n)%n222,           &
          LIS_rc%udef,iret)

  elseif ((trim(LIS_rc%met_interp(findex)).eq."bilinear").or.          &
       (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")) then
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,                &
          merra2biasforc(var_index,:),                                 &
          merra2bias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),          &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          merra2bias_struc(n)%w111,merra2bias_struc(n)%w121,           &
          merra2bias_struc(n)%w211,merra2bias_struc(n)%w221,           &
          merra2bias_struc(n)%n111,merra2bias_struc(n)%n121,           &
          merra2bias_struc(n)%n211,merra2bias_struc(n)%n221,           &
          LIS_rc%udef,iret)

  elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,                &
          merra2biasforc(var_index,:),                                 &
          merra2bias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),          &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          merra2bias_struc(n)%n113,LIS_rc%udef,iret)

  else
     write(LIS_logunit,*) '[ERR] Spatial interpolation option '//      &
          trim(LIS_rc%met_interp(findex))//                            &
          ' not supported for MERRA2bias.'
     call LIS_endrun()
  endif

  if (pcp_flag.and.(LIS_rc%pcp_downscale(findex).ne.0)) then
     call LIS_pcpClimoDownscaling(n,findex,month,                      &
          LIS_rc%lnc(n)*LIS_rc%lnr(n),merra2biasforc(var_index,:),lo)
  endif

end subroutine interp_merra2bias_var

!BOP
!
! !ROUTINE: interp_lapserate_merra2bias
! \label{interp_lapserate_merra2bias}
!
! !INTERFACE:
subroutine interp_lapserate_merra2bias(n,order,input_var)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use merra2bias_forcingMod, only : merra2bias_struc
  use LIS_constantsMod, only : LIS_CONST_LAPSE_RATE
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n, order
  real,    intent(in)    :: input_var(merra2bias_struc(n)%ncold, &
       merra2bias_struc(n)%nrold)

!
! !DESCRIPTION:
!  This subroutine spatially interpolates a MERRA2 Bias-corrected field
!  to the LIS running domain
!
!EOP
  integer   :: c,r,k,iret
  integer   :: gid
  real      :: f (merra2bias_struc(n)%ncold*merra2bias_struc(n)%nrold)
  logical*1 :: lb(merra2bias_struc(n)%ncold*merra2bias_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size
  real      :: output_var(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  external :: bilinear_interp
! _____________________________________________________________

  input_size = merra2bias_struc(n)%ncold*merra2bias_struc(n)%nrold

  lb = .true.
  do r=1,merra2bias_struc(n)%nrold
     do c=1,merra2bias_struc(n)%ncold
        k= c+(r-1)*merra2bias_struc(n)%ncold
        f(k) = input_var(c,r)
        if ( f(k) == 1.e+15 ) then
           f(k)  = LIS_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo

  call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,             &
       output_var(:),                                            &
       merra2bias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),       &
       LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
       merra2bias_struc(n)%w111,merra2bias_struc(n)%w121,        &
       merra2bias_struc(n)%w211,merra2bias_struc(n)%w221,        &
       merra2bias_struc(n)%n111,merra2bias_struc(n)%n121,        &
       merra2bias_struc(n)%n211,merra2bias_struc(n)%n221,        &
       LIS_rc%udef, iret)

  if (order.eq.1) then
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              gid = LIS_domain(n)%gindex(c,r)
              if(.not.isNaN(output_var(c+(r-1)*LIS_rc%lnc(n)))) then
                 merra2bias_struc(n)%lapserate1(gid) = &
                      output_var(c+(r-1)*LIS_rc%lnc(n))/1000.0
              else
                 merra2bias_struc(n)%lapserate1(gid) = &
                      LIS_CONST_LAPSE_RATE
              endif
           endif
        enddo
     enddo
  endif
  if (order.eq.2) then
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              gid = LIS_domain(n)%gindex(c,r)
              if(.not.isNaN(output_var(c+(r-1)*LIS_rc%lnc(n)))) then
                 merra2bias_struc(n)%lapserate2(gid) = &
                      output_var(c+(r-1)*LIS_rc%lnc(n))/1000.0
              else
                 merra2bias_struc(n)%lapserate2(gid) = &
                      LIS_CONST_LAPSE_RATE
              endif
           endif
        enddo
     enddo
  endif

end subroutine interp_lapserate_merra2bias
