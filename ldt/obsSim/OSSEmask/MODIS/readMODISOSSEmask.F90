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
!
! !ROUTINE: readMODISOSSEmask
! \label{readMODISOSSEmask}
!
! !INTERFACE:
subroutine readMODISOSSEmask(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
#if (defined USE_GRIBAPI)
  use grib_api
#endif
  use LDT_coreMod
  use LDT_obsSimMod
  use LDT_historyMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use MODISOSSEmask_Mod
!
! !DESCRIPTION:
!  This routine reads the MODIS orbital masks and interpolates to the
!  LDT domain
!EOP
  implicit none

  integer,   intent(in) :: n
  integer,   parameter  :: ndays = 16
  character(len=LDT_CONST_PATH_LEN)         :: fname
  integer               :: c,r,k
  integer               :: maskid
  logical               :: file_exists
  integer               :: ftn
  integer               :: iret
  integer               :: tindex
  integer               :: lat_off,lon_off
  real                  :: mask(MODISOSSEmaskData%nr,&
       MODISOSSEmaskData%nc,ndays)
  real                  :: mask1d(MODISOSSEmaskData%nc*MODISOSSEmaskData%nr)
  real                  :: mask_2d(LDT_rc%lnc(n),LDT_rc%lnr(n))

  ! create MODIS filename
  call create_MODIS_ossemask_filename(&
       MODISOSSEmaskData%odir,&
       MODISOSSEmaskData%datares,&
       fname)

! read, subset and interpolate the data

  inquire(file=trim(fname),exist=file_exists)

  if(file_exists) then
     write(LDT_logunit,*) '[INFO] reading MODIS mask ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

     iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
     call LDT_verify(iret, 'Error opening file '//trim(fname))

     call LDT_verify(nf90_inq_varid(ftn, "BinaryViewMask", maskid),&
          'nf90_inq_varid failed in readMODISOSSEmask')

     lat_off = nint((MODISOSSEmaskData%gridDesci(4)+90)/0.01)+1
     lon_off = nint((MODISOSSEmaskData%gridDesci(5)+180)/0.01)+1

     call LDT_verify(nf90_get_var(ftn,maskid,mask,&
          start=(/lat_off,lon_off,1/),&
          count=(/MODISOSSEmaskData%nr,MODISOSSEmaskData%nc,ndays/)),&
          'nf90_get_Var failed in readMODISOSSEmask')

     !subset data to the current day
     !index 1 in the data (october 1) is doy 274. The data repeats
     !every 16 days.
     !ignoring leap years for now

     if(LDT_rc%doy.ge.274) then
        tindex = mod((LDT_rc%doy - 274),7)+1
     else
        tindex = mod(LDT_rc%doy +94, 7)+1
     endif

     mask1d(:) = 0.0

     do r=1,MODISOSSEmaskData%nr
        do c=1, MODISOSSEmaskData%nc
           mask1d(c+(r-1)*MODISOSSEmaskData%nc) = mask(r,c,tindex)
        enddo
     enddo

     call convertMODISOSSEmaskToLDTgrid(n,mask1d(:),mask_2d(:,:))

#endif
  else
     write(LDT_logunit,*) '[WARN] file '//trim(fname)
     write(LDT_logunit,*) '[WARN] not found ...'
     mask_2d = LDT_rc%udef
  endif

  LDT_obsSim_struc%datamask = mask_2d

end subroutine readMODISOSSEmask

!BOP
!
! !ROUTINE: create_MODIS_ossemask_filename
! \label{create_MODIS_ossemask_filename}
!
! !INTERFACE:
subroutine create_MODIS_ossemask_filename(odir,datares, fname)
  !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod

   implicit none
! !ARGUMENTS:

   character(len=*)             :: odir
   real                         :: datares
   character(len=*)             :: fname
!
! !DESCRIPTION:
!  Create the file name for the MODIS mask files
!
!
!EOP

   if(datares.eq.0.01) then
      fname = trim(odir)//'/MODIS_sw2330km_res1km.nc'
   elseif(datares.eq.0.05) then
      fname = trim(odir)//'/MODIS_sw2330km_res5km.nc'
   elseif(datares.eq.0.25) then
      fname = trim(odir)//'/MODIS_sw2330km_res25km.nc'
   endif
 end subroutine create_MODIS_ossemask_filename

!BOP
!
! !ROUTINE: convertMODISOSSEmaskToLDTgrid
! \label{convertMODISOSSEmaskToLDTgrid}
!
! !INTERFACE:
 subroutine convertMODISOSSEmaskToLDTgrid(n, var_inp, var_out)
! !USES:
   use LDT_coreMod
   use MODISOSSEmask_Mod

   implicit none
! !ARGUMENTS:
   integer         :: n
   real            :: var_inp(MODISOSSEmaskData%nc*MODISOSSEmaskData%nr)
   real            :: var_out(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This routine interpolates or upscales the input data to
!  the LDT grid. If the input data is finer than the LDT
!  grid, the input data is upscaled. If the input data is
!  coarser, then it is interpolated to the LDT grid.
!
!EOP
   integer         :: ios
   integer         :: c,r
   logical*1       :: var_data_b(MODISOSSEmaskData%nc*MODISOSSEmaskData%nr)
   real            :: varobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   logical*1       :: varobs_b_ip(MODISOSSEmaskData%nc*MODISOSSEmaskData%nr)

   do r=1,MODISOSSEmaskData%nr
      do c=1, MODISOSSEmaskData%nc
         if(var_inp(c+(r-1)*MODISOSSEmaskData%nc).ne.LDT_rc%udef) then
            var_data_b(c+(r-1)*MODISOSSEmaskData%nc) = .true.
         else
            var_data_b(c+(r-1)*MODISOSSEmaskData%nc) = .false.
         endif
      enddo
   enddo

   if(LDT_isLDTatAfinerResolution(n,MODISOSSEmaskData%datares)) then

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!--------------------------------------------------------------------------
      call bilinear_interp(LDT_rc%gridDesc(n,:),&
           var_data_b, var_inp, varobs_b_ip, varobs_ip, &
           MODISOSSEmaskData%nc*MODISOSSEmaskData%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           MODISOSSEmaskData%w11, MODISOSSEmaskData%w12, &
           MODISOSSEmaskData%w21, MODISOSSEmaskData%w22, &
           MODISOSSEmaskData%n11, MODISOSSEmaskData%n12, &
           MODISOSSEmaskData%n21, MODISOSSEmaskData%n22, &
           LDT_rc%udef, ios)

      call neighbor_interp(LDT_rc%gridDesc(n,:),&
           var_data_b, var_inp, varobs_b_ip, varobs_ip, &
           MODISOSSEmaskData%nc*MODISOSSEmaskData%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           MODISOSSEmaskData%n11, LDT_rc%udef, ios)
   else
      call upscaleByAveraging(&
           MODISOSSEmaskData%nc*MODISOSSEmaskData%nr,&
           LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
           MODISOSSEmaskData%n11,var_data_b, var_inp, varobs_b_ip,varobs_ip)

   endif

   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if(varobs_b_ip(c+(r-1)*LDT_rc%lnc(n))) then
            var_out(c,r) = varobs_ip(c+(r-1)*LDT_rc%lnc(n))
         else
            var_out(c,r) = LDT_rc%udef
         endif
      enddo
   enddo

 end subroutine convertMODISOSSEmaskToLDTgrid
