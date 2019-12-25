!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_era5
! \label{read_era5}
! 
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:
subroutine read_era5(n, kk,order, month, findex,          &
     fname, ferror)
! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod
  use LDT_FORC_AttributesMod
  use LDT_metforcingMod, only : LDT_forc
  use era5_forcingMod, only : era5_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: kk
  integer, intent(in)          :: order
  integer, intent(in)          :: month
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: fname
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  ERA5 data, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain. \newline
!
!  ZS: surface orography (m)
!  UREF: reference height for the wind (m)
!  ZREF: reference height for the temperature (m)
!
!  Tair: Atmospheric temperature (K) (2m)
!  Qair: Atmospheric humidity (kg/kg)
!  PSurf: Atmospheric pressure (Pa)
!  Rainf: Rain (kg/m2/s)
!  Snowf: Snow (kg/m2/s)
!  Wind: Wind speed (m/s) (10m)
!  Wind_DIR: Wind direction (degrees from N, clockwise)
!  LWdown: Long-wave radiation (W/m2)
!  DIR_SWdown: direct short-wave radiation (W/m2)
!  SCA_SWdown: diffuse short-wave radiation (W/m2)
!  CO2air: near surface CO2 concentration (kg/m3)
!   
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 1 hour ERA5 analysis file
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
  
  integer   :: ftn
  integer   :: tmpId, qId, windId, lwdId,psId,rainfId,snowfID,dirSWId,difSWId
  integer   :: tindex
  logical   :: file_exists
  integer   :: c,r,t,k,iret
  integer   :: mo
  logical   :: read_lnd

  real      :: tair(era5_struc(n)%npts)
  real      :: qair(era5_struc(n)%npts)
  real      :: wind(era5_struc(n)%npts)
  real      :: ps(era5_struc(n)%npts)
  real      :: rainf(era5_struc(n)%npts)
  real      :: snowf(era5_struc(n)%npts)
  real      :: dirswd(era5_struc(n)%npts)
  real      :: difswd(era5_struc(n)%npts)
  real      :: swd(era5_struc(n)%npts)
  real      :: lwd(era5_struc(n)%npts)
! __________________________________________________________________________

#if (defined USE_NETCDF4) 
  ferror = 1 

  mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

! Read single layer file (*slv) fields:
  inquire(file=fname,exist=file_exists) 
  if(file_exists) then 
     write(LDT_logunit,*)'[INFO] Reading ERA5 file (bookend,', order,' -',trim(fname), ')'
     call LDT_verify(nf90_open(path=trim(fname), mode=NF90_NOWRITE, &
          ncid=ftn), 'nf90_open failed in read_era5')

     call LDT_verify(nf90_inq_varid(ftn,'PSurf',psId), &
          'nf90_inq_varid failed for psurf in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
          'nf90_inq_varid failed for Tair in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'Qair',qId), &
          'nf90_inq_varid failed for Qair in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'Wind',windId), &
          'nf90_inq_varid failed for Wind in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
          'nf90_inq_varid failed for Rainf in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'Snowf',snowfId), &
          'nf90_inq_varid failed for Snowf in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'DIR_SWdown',dirSWId), &
          'nf90_inq_varid failed for DIR_SWdown in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'SCA_SWdown',difSWId), &
          'nf90_inq_varid failed for SCA_SWdown in read_era5')
     call LDT_verify(nf90_inq_varid(ftn,'LWdown',lwdId), &
          'nf90_inq_varid failed for LWdown in read_era5')

     if(LDT_rc%da.gt.1) then 
        tindex = (LDT_rc%da - 1)*24 + LDT_rc%hr + 1
     else
        tindex = LDT_rc%hr + 1
     endif

     call LDT_verify(nf90_get_var(ftn,psId, ps, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for ps in read_era5') 
     call LDT_verify(nf90_get_var(ftn,tmpId, tair, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for tair in read_era5') 
     call LDT_verify(nf90_get_var(ftn,qId, qair, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for qair in read_era5') 
     call LDT_verify(nf90_get_var(ftn,windId, wind, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for wind in read_era5') 
     call LDT_verify(nf90_get_var(ftn,rainfId, rainf, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for rainf in read_era5') 
     call LDT_verify(nf90_get_var(ftn,snowfId, snowf, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for snowf in read_era5') 
     call LDT_verify(nf90_get_var(ftn,dirSWId, dirswd, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for dirSW in read_era5') 
     call LDT_verify(nf90_get_var(ftn,difSWId, difswd, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for difSW in read_era5') 
     call LDT_verify(nf90_get_var(ftn,lwdId, lwd, &
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for lwd in read_era5') 

     call LDT_verify(nf90_close(ftn), &
          'failed to close file in read_era5')

     do t=1,era5_struc(n)%npts
        swd(t) = dirswd(t) + difswd(t)
     enddo

     call interp_era5_var(n,kk,order,findex,month,tair,  1, .false.)
     call interp_era5_var(n,kk,order,findex,month,qair,  2, .false.)
     call interp_era5_var(n,kk,order,findex,month,swd,   3, .false.)
     call interp_era5_var(n,kk,order,findex,month,lwd,   4, .false.)
     call interp_era5_var(n,kk,order,findex,month,wind,  5, .false.)
     call interp_era5_var(n,kk,order,findex,month,ps,    6, .false.)
     call interp_era5_var(n,kk,order,findex,month,rainf, 7, .true.)
     call interp_era5_var(n,kk,order,findex,month,snowf, 8, .true.)

  else
     write(LDT_logunit,*) '[ERR] ',trim(fname)//' does not exist'
     call LDT_endrun()
  endif

  
#endif
end subroutine read_era5

!BOP
! 
! !ROUTINE: interp_era5_var
! \label{interp_era5_var}
! 
! !INTERFACE: 
subroutine interp_era5_var(n,kk,order,findex, month, input_var,  var_index, &
     pcp_flag)

! !USES: 
  use LDT_coreMod
  use LDT_logMod
  use era5_forcingMod, only : era5_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: kk
  integer, intent(in)    :: order
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(era5_struc(n)%npts)
  integer, intent(in)    :: var_index
  logical, intent(in)    :: pcp_flag

!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a ERA5 field
!  to the LDT running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  integer   :: doy
  integer   :: ftn
  integer   :: pcp1Id, pcp2Id, pcp3Id, pcp4Id,pcp5Id, pcp6Id
  real      :: f (era5_struc(n)%ncold*era5_struc(n)%nrold)
  logical*1 :: lb(era5_struc(n)%ncold*era5_struc(n)%nrold)
  logical*1 :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer   :: input_size
  real      :: era5forc(LDT_rc%lnc(n)*LDT_rc%lnr(n))
! _____________________________________________________________

  input_size = era5_struc(n)%ncold*era5_struc(n)%nrold

!-----------------------------------------------------------------------    
! Apply corrections
!-----------------------------------------------------------------------  
     
  lb = .false.
  do r=1,era5_struc(n)%nrold
     do c=1,era5_struc(n)%ncold           
        k= c+(r-1)*era5_struc(n)%ncold
        if(era5_struc(n)%g2p(c,r).gt.0) then 
           f(k) = input_var(era5_struc(n)%g2p(c,r))
           lb(k) = .true.
        else
           f(k) = LDT_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo
  
  if(pcp_flag.and.&
       trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
     
     call conserv_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          era5forc(:), &
          era5_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),& 
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5_struc(n)%w112,era5_struc(n)%w122,&
          era5_struc(n)%w212,era5_struc(n)%w222,&
          era5_struc(n)%n112,era5_struc(n)%n122,&
          era5_struc(n)%n212,era5_struc(n)%n222,&
          LDT_rc%udef, iret)
     
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear".or.&
       trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
     call bilinear_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          era5forc(:), &
          era5_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n), & 
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5_struc(n)%w111,era5_struc(n)%w121,&
          era5_struc(n)%w211,era5_struc(n)%w221,&
          era5_struc(n)%n111,era5_struc(n)%n121,&
          era5_struc(n)%n211,era5_struc(n)%n221,&
          LDT_rc%udef, iret)
     
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then 
     call neighbor_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          era5forc(:),era5_struc(n)%mi,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n),&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5_struc(n)%n113,LDT_rc%udef,iret)
  else
     write(LDT_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(LDT_rc%met_gridtransform(findex))//&
          ' not supported for ERA5'
     call LDT_endrun()
  endif
  
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(LDT_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              era5_struc(n)%metdata1(var_index,&
                   LDT_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LDT_rc%lnc(n))
           elseif(order.eq.2) then 
              era5_struc(n)%metdata2(var_index,&
                   LDT_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LDT_rc%lnc(n))
           endif
        endif
     enddo
  enddo

end subroutine interp_era5_var




