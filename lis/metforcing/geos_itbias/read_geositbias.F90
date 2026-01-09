!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_geositbias
! \label{read_geositbias}
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it)
! 07 Jan 2026: Kristen Whitney, initial code for using dynamic lapse rate
! 09 Jan 2026: Eric Kemp, reformatted.
!
! !INTERFACE:
subroutine read_geositbias(n,order,month,findex,                     &
     geosname,                                                       &
     lapseratefname,                                                 &
     geositbiasforc,ferror)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc
  use geositbias_forcingMod, only : geositbias_struc
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
  character(len=*), intent(in) :: geosname
  character(len=*), intent(in) :: lapseratefname
  real, intent(inout)          :: &
       geositbiasforc(geositbias_struc(n)%nvars,&
       LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer, intent(out)         :: ferror

! !DESCRIPTION:
!  For the given time, reads parameters from GEOS-ITbias data, transforms
!  into 9 LIS forcing parameters and interpolates to the LIS domain. \newline
!
! GEOS-ITbias FORCING VARIABLES: \newline
!  1. tair    Temperature of air at 2-m, 10-m, or LML [$K$] \newline
!  2. qair    Specific humidity of air at 2-m, 10-m, or LML [$kg/kg$] \newline
!  3. swgdn   Downward shortwave radiation at the ground [$W/m^2$] \newline
!  4. lwgab   Downward longwave radiation at the ground [$W/m^2$] \newline
!  5. uwind   Zonal wind at 2-m, 10-m, or LML [$m/s$] \newline
!  6. vwind   Meridional wind at 2-m, 10-m, or LML [$m/s$] \newline
!  7. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  8. pretot  Total precipitation [$mm/s$] \newline
!  9. precon  Convective precipitation [$mm/s$] \newline
! 10. presno  Precipitation falling as snow [$mm/s$] \newline
! 11. swland  Net shortwave radiation at the ground [$W/m^2$] \newline
! 12. pardr   Surface downward PAR beam flux [$W/m^2$] \newline
! 13. pardf   Surface downward PAR diffuse flux [$W/m^2$] \newline
! 14. hlml    Height of center of lowest model layer (LML) [$m$] \newline
!
!  The arguments are:
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 1 hour GEOS-ITbias analysis file
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

  integer   :: ftn_geos
  integer   :: tmpId,qId,psId
  integer   :: lwgabId
  integer   :: nr_index,nc_index
  logical   :: file_exists
  integer   :: mo
  real      :: tair(geositbias_struc(n)%ncold,geositbias_struc(n)%nrold)
  real      :: qair(geositbias_struc(n)%ncold,geositbias_struc(n)%nrold)
  real      :: ps(geositbias_struc(n)%ncold,geositbias_struc(n)%nrold)
  real      :: lwgab(geositbias_struc(n)%ncold,geositbias_struc(n)%nrold)
  integer                           :: ftn_drate
  integer                           :: lapserateid
  real, allocatable                 :: lapse_rate_in(:,:)
  integer :: nc_sub, nr_sub
  real, allocatable :: lapse_rate_sub(:,:)
  real, allocatable :: lat_sub(:), lon_sub(:)
  real, allocatable :: lat_full(:), lon_full(:)
  integer :: i0, j0
  integer :: latid, lonid

  external :: interp_geositbias_var
  external :: interp_lapserate_geositbias
! __________________________________________________________________________

#if (defined USE_NETCDF3)
  write(LIS_logunit,*) "[ERR] GEOS-ITbias reader requires NetCDF4"
  call LIS_endrun()
#endif

#if (defined USE_NETCDF4)
  ferror = 0
  nr_index = geositbias_struc(n)%nrold
  nc_index = geositbias_struc(n)%ncold
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

! Read single layer file (*slv) fields:
  inquire(file=geosname,exist=file_exists)
  if(file_exists) then
     write(LIS_logunit,*)                                           &
          '[INFO] Reading GEOS-ITbias file: ',trim(geosname)
     call LIS_verify(nf90_open(path=trim(geosname),                 &
          mode=NF90_NOWRITE,ncid=ftn_geos),                         &
          'nf90_open failed for geosfile in read_geositbias')

     call LIS_verify(nf90_inq_varid(ftn_geos,'PS',psId),            &
          'nf90_inq_varid failed for ps in read_geositbias')
     call LIS_verify(nf90_get_var(ftn_geos,psId,ps),                &
          'nf90_get_var failed for ps in read_geositbias')

     call LIS_verify(nf90_inq_varid(ftn_geos,'LWGAB',lwgabId),      &
          'nf90_inq_varid failed for LWGAB in read_geositbias')
     call LIS_verify(nf90_get_var(ftn_geos,lwgabId,lwgab),          &
          'nf90_get_var failed for LWGAB in read_geositbias')

     call LIS_verify(nf90_inq_varid(ftn_geos,'QV2M',qId),           &
          'nf90_inq_varid failed for QV2M in read_geositbias')
     call LIS_verify(nf90_get_var(ftn_geos,qId,qair),               &
          'nf90_get_var failed for QV2M in read_geositbias')

     call LIS_verify(nf90_inq_varid(ftn_geos,'T2M',tmpId),          &
          'nf90_inq_varid failed for T2M in read_geositbias')
     call LIS_verify(nf90_get_var(ftn_geos,tmpId,Tair),             &
          'nf90_get_var failed for T2M in read_geositbias')

     call LIS_verify(nf90_close(ftn_geos),                          &
          '[WARN] Failed to close slvfile in read_geositbias.')

     call interp_geositbias_var(n,findex,month,tair,  1,.false.,    &
          geositbiasforc)
     call interp_geositbias_var(n,findex,month,qair,  2,.false.,    &
          geositbiasforc)
     call interp_geositbias_var(n,findex,month,ps,    7,.false.,    &
          geositbiasforc)
     call interp_geositbias_var(n,findex,month,lwgab,4,.false.,     &
          geositbiasforc)

    ! Obtain lapse rates
     if ((geositbias_struc(n)%usedynlapserate.eq.1).and.            &
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
           allocate(lapse_rate_in(geositbias_struc(n)%ncold,   &
                geositbias_struc(n)%nrold))
           lapse_rate_in = 1.e+15
           ! Case 1: full global lapse-rate file
           if (nc_sub == geositbias_struc(n)%ncold .and.       &
                nr_sub == geositbias_struc(n)%nrold) then
              call LIS_verify(nf90_get_var(ftn_drate,lapserateid, &
                   lapse_rate_in), 'get lapse_rate failed')

              ! Case 2: subdomain lapse-rate file
           else if (nc_sub <= geositbias_struc(n)%ncold .and.  &
                nr_sub <= geositbias_struc(n)%nrold) then
              allocate(lapse_rate_sub(nc_sub,nr_sub))
              allocate(lat_sub(nr_sub), lon_sub(nc_sub))
              allocate(lat_full(geositbias_struc(n)%nrold), &
                   lon_full(geositbias_struc(n)%ncold))

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
              call LIS_verify(nf90_open(path=trim(geosname),         &
                   mode=NF90_NOWRITE,ncid=ftn_geos),                 &
                   'nf90_open failed for geosfile in read_geositbias')
              call LIS_verify(nf90_inq_varid(ftn_geos,'lat',latid),  &
                   'inq lat varid failed (full)')
              call LIS_verify(nf90_inq_varid(ftn_geos,'lon',lonid),  &
                   'inq lon varid failed (full)')
              call LIS_verify(nf90_get_var(ftn_geos,latid,lat_full), &
                   'get lat_full failed')
              call LIS_verify(nf90_get_var(ftn_geos,lonid,lon_full), &
                   'get lon_full failed')
              call LIS_verify(nf90_close(ftn_geos), &
                   '[WARN] Failed to close slvfile in read_geositbias.')

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
           call interp_lapserate_geositbias(n,order,lapse_rate_in)

           deallocate(lapse_rate_in)

        else
           write(LIS_logunit,*) '[WARN] Could not find ', &
                trim(lapseratefname)
           write(LIS_logunit,*) '[WARN] Using static lapse rate.'
           geositbias_struc(n)%lapserate1 = LIS_CONST_LAPSE_RATE
           geositbias_struc(n)%lapserate2 = LIS_CONST_LAPSE_RATE
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
end subroutine read_geositbias

!BOP
!
! !ROUTINE: interp_geositbias_var
! \label{interp_geositbias_var}
!
! !INTERFACE:
subroutine interp_geositbias_var(n,findex,month,input_var,var_index, &
     pcp_flag,geositbiasforc)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use geositbias_forcingMod, only : geositbias_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(geositbias_struc(n)%ncold,       &
       geositbias_struc(n)%nrold)
  integer, intent(in)    :: var_index
  logical, intent(in)    :: pcp_flag
  real,    intent(inout) :: geositbiasforc(geositbias_struc(n)%nvars,  &
       LIS_rc%lnc(n)*LIS_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine spatially interpolates a GEOS-ITbias field
!  to the LIS running domain
!
!EOP
  integer   :: c,r,k,iret
  real      :: f (geositbias_struc(n)%ncold*geositbias_struc(n)%nrold)
  logical*1 :: lb(geositbias_struc(n)%ncold*geositbias_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size

  external :: conserv_interp
  external :: bilinear_interp
  external :: neighbor_interp

! _____________________________________________________________

  input_size = geositbias_struc(n)%ncold*geositbias_struc(n)%nrold

!-----------------------------------------------------------------------
! Apply downscaling
!-----------------------------------------------------------------------
  lb = .true.
  do r = 1,geositbias_struc(n)%nrold
     do c = 1,geositbias_struc(n)%ncold
        k = c+(r-1)*geositbias_struc(n)%ncold
        f(k) = input_var(c,r)
        if (f(k).eq.1.e+15) then
           f(k)  = LIS_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo

  if (pcp_flag.and.(LIS_rc%pcp_downscale(findex).ne.0)) then
     ! input_data becomes the ratio field.
     call LIS_generatePcpClimoRatioField(n,findex,"GEOS-ITbias",       &
          month,input_size,f,lb)
  endif

  if (pcp_flag.and.                                                    &
       (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")) then

     call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,                 &
          geositbiasforc(var_index,:),                                 &
          geositbias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),          &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          geositbias_struc(n)%w112,geositbias_struc(n)%w122,           &
          geositbias_struc(n)%w212,geositbias_struc(n)%w222,           &
          geositbias_struc(n)%n112,geositbias_struc(n)%n122,           &
          geositbias_struc(n)%n212,geositbias_struc(n)%n222,           &
          LIS_rc%udef,iret)

  elseif ((trim(LIS_rc%met_interp(findex)).eq."bilinear").or.          &
       (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")) then
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,                &
          geositbiasforc(var_index,:),                                 &
          geositbias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),          &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          geositbias_struc(n)%w111,geositbias_struc(n)%w121,           &
          geositbias_struc(n)%w211,geositbias_struc(n)%w221,           &
          geositbias_struc(n)%n111,geositbias_struc(n)%n121,           &
          geositbias_struc(n)%n211,geositbias_struc(n)%n221,           &
          LIS_rc%udef,iret)

  elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,                &
          geositbiasforc(var_index,:),                                 &
          geositbias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),          &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          geositbias_struc(n)%n113,LIS_rc%udef,iret)

  else
     write(LIS_logunit,*) '[ERR] Spatial interpolation option '//      &
          trim(LIS_rc%met_interp(findex))//                            &
          ' not supported for GEOS-ITbias.'
     call LIS_endrun()
  endif

  if (pcp_flag.and.(LIS_rc%pcp_downscale(findex).ne.0)) then
     call LIS_pcpClimoDownscaling(n,findex,month,                      &
          LIS_rc%lnc(n)*LIS_rc%lnr(n),geositbiasforc(var_index,:),lo)
  endif

end subroutine interp_geositbias_var

!BOP
!
! !ROUTINE: interp_lapserate_geositbias
! \label{interp_lapserate_geositbias}
!
! !INTERFACE:
subroutine interp_lapserate_geositbias(n,order,input_var)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use geositbias_forcingMod, only : geositbias_struc
  use LIS_constantsMod, only : LIS_CONST_LAPSE_RATE
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n, order
  real,    intent(in)    :: input_var(geositbias_struc(n)%ncold, &
       geositbias_struc(n)%nrold)

!
! !DESCRIPTION:
!  This subroutine spatially interpolates a GEOS-IT Bias-corrected field
!  to the LIS running domain
!
!EOP

  integer   :: c,r,k,iret
  integer   :: gid
  real      :: f (geositbias_struc(n)%ncold*geositbias_struc(n)%nrold)
  logical*1 :: lb(geositbias_struc(n)%ncold*geositbias_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size
  real      :: output_var(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  external :: bilinear_interp
! _____________________________________________________________

  input_size = geositbias_struc(n)%ncold*geositbias_struc(n)%nrold

  lb = .true.
  do r=1,geositbias_struc(n)%nrold
     do c=1,geositbias_struc(n)%ncold
        k= c+(r-1)*geositbias_struc(n)%ncold
        f(k) = input_var(c,r)
        if ( f(k) == 1.e+15 ) then
           f(k)  = LIS_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo

  call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,             &
       output_var(:),                                            &
       geositbias_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),       &
       LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
       geositbias_struc(n)%w111,geositbias_struc(n)%w121,        &
       geositbias_struc(n)%w211,geositbias_struc(n)%w221,        &
       geositbias_struc(n)%n111,geositbias_struc(n)%n121,        &
       geositbias_struc(n)%n211,geositbias_struc(n)%n221,        &
       LIS_rc%udef, iret)

  if (order.eq.1) then
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then
              gid = LIS_domain(n)%gindex(c,r)
              if(.not.isNaN(output_var(c+(r-1)*LIS_rc%lnc(n)))) then
                 geositbias_struc(n)%lapserate1(gid) = &
                      output_var(c+(r-1)*LIS_rc%lnc(n))/1000.0
              else
                 geositbias_struc(n)%lapserate1(gid) = &
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
                 geositbias_struc(n)%lapserate2(gid) = &
                      output_var(c+(r-1)*LIS_rc%lnc(n))/1000.0
              else
                 geositbias_struc(n)%lapserate2(gid) = &
                      LIS_CONST_LAPSE_RATE
              endif
           endif
        enddo
     enddo
  endif

end subroutine interp_lapserate_geositbias
