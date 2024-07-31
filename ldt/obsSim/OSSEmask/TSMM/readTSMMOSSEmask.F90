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
! !ROUTINE: readTSMMOSSEmask
! \label{readTSMMOSSEmask}
!
! !INTERFACE: 
subroutine readTSMMOSSEmask(n)
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
  use TSMMOSSEmask_Mod
!
! !DESCRIPTION: 
!  This routine reads the TSMM orbital masks and interpolates to the
!  LDT domain
!EOP
  implicit none

  integer,   intent(in) :: n
  integer,   parameter  :: ndays = 7
  character(len=LDT_CONST_PATH_LEN)         :: fname
  integer               :: c,r,k
  integer               :: maskid
  logical               :: file_exists
  integer               :: ftn
  integer               :: iret
  integer               :: tindex
  integer               :: lat_off,lon_off
  real                  :: mask(TSMMOSSEmaskData%nr,&
       TSMMOSSEmaskData%nc,ndays)
  real                  :: mask1d(TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr)
  real                  :: mask_2d(LDT_rc%lnc(n),LDT_rc%lnr(n))

  ! create TSMM filename
  call create_TSMM_ossemask_filename(&
       TSMMOSSEmaskData%odir,&
       TSMMOSSEmaskData%datares,&
       fname)

! read, subset and interpolate the data
  
  inquire(file=trim(fname),exist=file_exists)

  if(file_exists) then 
     write(LDT_logunit,*) '[INFO] reading TSMM mask ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
        
     iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
     call LDT_verify(iret, 'Error opening file '//trim(fname))
     
     call LDT_verify(nf90_inq_varid(ftn, "BinaryViewMask", maskid),&
          'nf90_inq_varid failed in readTSMMOSSEmask')

     lat_off = nint((TSMMOSSEmaskData%gridDesci(4)+90)/0.01)+1
     lon_off = nint((TSMMOSSEmaskData%gridDesci(5)+180)/0.01)+1
     
     call LDT_verify(nf90_get_var(ftn,maskid,mask,&
          start=(/lat_off,lon_off,1/),&
          count=(/TSMMOSSEmaskData%nr,TSMMOSSEmaskData%nc,ndays/)),&
          'nf90_get_Var failed in readTSMMOSSEmask')
     
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
     
     do r=1,TSMMOSSEmaskData%nr
        do c=1, TSMMOSSEmaskData%nc
           mask1d(c+(r-1)*TSMMOSSEmaskData%nc) = mask(r,c,tindex)
        enddo
     enddo

     call convertTSMMOSSEmaskToLDTgrid(n,mask1d(:),mask_2d(:,:))
     
#endif
  else
     write(LDT_logunit,*) '[WARN] file '//trim(fname)
     write(LDT_logunit,*) '[WARN] not found ...'
     mask_2d = LDT_rc%udef
  endif

  LDT_obsSim_struc%datamask = mask_2d

end subroutine readTSMMOSSEmask

!BOP
!
! !ROUTINE: create_TSMM_ossemask_filename
! \label{create_TSMM_ossemask_filename}
!
! !INTERFACE:
subroutine create_TSMM_ossemask_filename(odir,datares, fname)
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
!  Create the file name for the TSMM mask files
!
!
!EOP

   if(datares.eq.0.01) then
      fname = trim(odir)//'/TSMM_sw250km_res1km.nc'
   elseif(datares.eq.0.05) then
      fname = trim(odir)//'/TSMM_sw250km_res5km.nc'
   elseif(datares.eq.0.25) then
      fname = trim(odir)//'/TSMM_sw250km_res25km.nc'
   endif
 end subroutine create_TSMM_ossemask_filename

!BOP
! 
! !ROUTINE: convertTSMMOSSEmaskToLDTgrid
! \label{convertTSMMOSSEmaskToLDTgrid}
!
! !INTERFACE: 
 subroutine convertTSMMOSSEmaskToLDTgrid(n, var_inp, var_out)
! !USES:    
   use LDT_coreMod
   use TSMMOSSEmask_Mod

   implicit none
! !ARGUMENTS: 
   integer         :: n 
   real            :: var_inp(TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr)
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
   logical*1       :: var_data_b(TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr)
   real            :: varobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   logical*1       :: varobs_b_ip(TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr)

   do r=1,TSMMOSSEmaskData%nr
      do c=1, TSMMOSSEmaskData%nc
         if(var_inp(c+(r-1)*TSMMOSSEmaskData%nc).ne.LDT_rc%udef) then 
            var_data_b(c+(r-1)*TSMMOSSEmaskData%nc) = .true. 
         else
            var_data_b(c+(r-1)*TSMMOSSEmaskData%nc) = .false.
         endif
      enddo
   enddo

   if(LDT_isLDTatAfinerResolution(n,TSMMOSSEmaskData%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
      call bilinear_interp(LDT_rc%gridDesc(n,:),&
           var_data_b, var_inp, varobs_b_ip, varobs_ip, &
           TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           TSMMOSSEmaskData%w11, TSMMOSSEmaskData%w12, &
           TSMMOSSEmaskData%w21, TSMMOSSEmaskData%w22, &
           TSMMOSSEmaskData%n11, TSMMOSSEmaskData%n12, &
           TSMMOSSEmaskData%n21, TSMMOSSEmaskData%n22, &
           LDT_rc%udef, ios)

      call neighbor_interp(LDT_rc%gridDesc(n,:),&
           var_data_b, var_inp, varobs_b_ip, varobs_ip, &
           TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           TSMMOSSEmaskData%n11, LDT_rc%udef, ios)
   else
      call upscaleByAveraging(&
           TSMMOSSEmaskData%nc*TSMMOSSEmaskData%nr,&
           LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
           TSMMOSSEmaskData%n11,var_data_b, var_inp, varobs_b_ip,varobs_ip)
      
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
         
 end subroutine convertTSMMOSSEmaskToLDTgrid

 
