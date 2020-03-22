!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_era5
! \label{read_era5}
! 
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:
subroutine read_era5(n, kk,order, year, month, day, hour, findex,          &
     fname, ferror)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc
  use era5_forcingMod, only : era5_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: kk
  integer, intent(in)          :: order
  integer, intent(in)          :: year
  integer, intent(in)          :: month
  integer, intent(in)          :: day
  integer, intent(in)          :: hour
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: fname
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  ERA5 data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain. \newline
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
!EOP
  
  integer   :: ftn
  integer   :: tmpId, qId, windId, lwdId,psId,rainfId,snowfID,dirSWId,difSWId
  integer   :: tindex
  logical   :: file_exists
  integer   :: c,r,t,k,l,iret
  integer   :: mo,rec_size
  logical   :: read_lnd
  logical   :: read_flag
  
  real      :: ps(era5_struc(n)%npts)
  real      :: tair(era5_struc(n)%npts)
  real      :: qair(era5_struc(n)%npts)
  real      :: swd(era5_struc(n)%npts)
  real      :: lwd(era5_struc(n)%npts)
  real      :: wind(era5_struc(n)%npts)
  real      :: rainf(era5_struc(n)%npts)
  real      :: snowf(era5_struc(n)%npts)
  real      :: dirsw(era5_struc(n)%npts)
  real      :: difsw(era5_struc(n)%npts)

  real      :: varfield(LIS_rc%lnc(n)*LIS_rc%lnr(n))

! __________________________________________________________________________

  ferror = 1

#if (defined USE_NETCDF4) 

  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

  inquire(file=fname,exist=file_exists) 
  if(file_exists) then 
     
     write(LIS_logunit,*)'[INFO] Reading ERA5 file ',trim(fname)

     tindex = (day - 1)*24 + hour + 1

     call LIS_verify(nf90_open(path=trim(fname), mode=NF90_NOWRITE, &
          ncid=ftn), 'nf90_open failed in read_era5')
     
     call LIS_verify(nf90_inq_varid(ftn,'PSurf',psId), &
          'nf90_inq_varid failed for psurf in read_era5')
     call LIS_verify(nf90_get_var(ftn,psId, ps,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for ps in read_era5') 
     
     if(order.eq.1) then 
        call interp_era5_var(n,findex,month,ps, .false.,&
             varfield)
     else
        call interp_era5_var(n,findex,month,ps, .false.,&
             varfield)
     endif
     call assign_processed_era5forc(n,kk,order,6,&
          varfield)

     call LIS_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
          'nf90_inq_varid failed for Tair in read_era5')
     call LIS_verify(nf90_get_var(ftn,tmpId, tair,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for tair in read_era5') 
     
     if(order.eq.1) then 
        call interp_era5_var(n,findex,month,tair, .false.,&
             varfield)
     else
        call interp_era5_var(n,findex,month,tair,  .false.,&
             varfield)
     endif
     call assign_processed_era5forc(n,kk,order,1,&
          varfield)

     call LIS_verify(nf90_inq_varid(ftn,'Qair',qId), &
          'nf90_inq_varid failed for Qair in read_era5')
     call LIS_verify(nf90_get_var(ftn,qId, qair,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for qair in read_era5') 
     
     if(order.eq.1) then 
        call interp_era5_var(n,findex,month,qair, .false.,&
             varfield)
     else
        call interp_era5_var(n,findex,month,qair,  .false.,&
             varfield)
     endif
     call assign_processed_era5forc(n,kk,order,2,&
          varfield)

     call LIS_verify(nf90_inq_varid(ftn,'Wind',windId), &
          'nf90_inq_varid failed for Wind in read_era5')
     call LIS_verify(nf90_get_var(ftn,windId, wind,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for wind in read_era5') 
     if(order.eq.1) then 
        call interp_era5_var(n,findex,month,wind, .false.,&
             varfield)
     else
        call interp_era5_var(n,findex,month,wind,  .false.,&
             varfield)
     endif
     call assign_processed_era5forc(n,kk,order,5,&
          varfield)

     call LIS_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
          'nf90_inq_varid failed for Rainf in read_era5')
     call LIS_verify(nf90_get_var(ftn,rainfId, rainf,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for rainf in read_era5') 
     call LIS_verify(nf90_inq_varid(ftn,'Snowf',snowfId), &
          'nf90_inq_varid failed for Snowf in read_era5')
     call LIS_verify(nf90_get_var(ftn,snowfId, snowf,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for snowf in read_era5') 

     do t=1,era5_struc(n)%npts
        if(rainf(t).ne.LIS_rc%udef.and.&
             snowf(t).ne.LIS_rc%udef) then 
           rainf(t) = rainf(t) + snowf(t)
        else
           rainf(t) = LIS_rc%udef
        endif
     enddo

     if(order.eq.1) then 
        call interp_era5_var(n,findex,month,rainf,.true.,&
             varfield)
     else
        call interp_era5_var(n,findex,month,rainf, .true.,&
             varfield)
     endif

     call assign_processed_era5forc(n,kk,order,7,&
          varfield)

     call LIS_verify(nf90_inq_varid(ftn,'DIR_SWdown',dirSWId), &
          'nf90_inq_varid failed for DIR_SWdown in read_era5')
     call LIS_verify(nf90_get_var(ftn,dirSWId, swd,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for dirSW in read_era5') 
     call LIS_verify(nf90_inq_varid(ftn,'SCA_SWdown',difSWId), &
          'nf90_inq_varid failed for SCA_SWdown in read_era5')
     call LIS_verify(nf90_get_var(ftn,difSWId, difsw,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for difSW in read_era5') 
     
     do t=1,era5_struc(n)%npts
        if(swd(t).ne.LIS_rc%udef.and.&
             difsw(t).ne.LIS_rc%udef) then 
           swd(t) = swd(t) + difsw(t)
        else
           swd(t) = LIS_rc%udef
        endif
     enddo

     if(order.eq.1) then 
        call interp_era5_var(n,findex,month,swd,  .false.,&
             varfield)
     else
        call interp_era5_var(n,findex,month,swd,   .false.,&
             varfield)
     endif

     call assign_processed_era5forc(n,kk,order,3,&
          varfield)

     call LIS_verify(nf90_inq_varid(ftn,'LWdown',lwdId), &
          'nf90_inq_varid failed for LWdown in read_era5')
     call LIS_verify(nf90_get_var(ftn,lwdId, lwd,&
          start=(/1,tindex/),count=(/era5_struc(n)%npts,1/)),&
          'nf90_get_var failed for lwd in read_era5')
     
     if(order.eq.1) then 
        call interp_era5_var(n,findex,month,lwd,  .false.,&
             varfield)
     else
        call interp_era5_var(n,findex,month,lwd,   .false.,&
             varfield)
     endif

     call assign_processed_era5forc(n,kk,order,4,&
          varfield)
        
     call LIS_verify(nf90_close(ftn), &
          'failed to close file in read_era5')
                
  else
     write(LIS_logunit,*) '[ERR] ',trim(fname)//' does not exist'
     call LIS_endrun()        
  endif
#endif

end subroutine read_era5


!BOP
! 
! !ROUTINE: interp_era5_var
! \label{interp_era5_var}
! 
! !INTERFACE: 
subroutine interp_era5_var(n,findex, month, input_var, &
     pcp_flag, output_var)

! !USES: 
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use era5_forcingMod, only : era5_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(era5_struc(n)%npts)
  real,    intent(out)   :: output_var(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical, intent(in)    :: pcp_flag

!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a single ERA5 field
!  to the LIS running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  integer   :: doy
  integer   :: ftn
  integer   :: pcp1Id, pcp2Id, pcp3Id, pcp4Id,pcp5Id, pcp6Id
  real      :: f (era5_struc(n)%ncold*era5_struc(n)%nrold)
  logical*1 :: lb(era5_struc(n)%ncold*era5_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size
! _____________________________________________________________

  input_size = era5_struc(n)%ncold*era5_struc(n)%nrold
  output_var = LIS_rc%udef

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
           f(k) = LIS_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo
  
!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    
     
  if(pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0) then 
     call LIS_generatePcpClimoRatioField(n,findex,"ERA5",&
          month, & 
          input_size, &
          f, &
          lb)     
  endif
          
  if(pcp_flag.and.&
       trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
     
     call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          era5_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5_struc(n)%w112,era5_struc(n)%w122,&
          era5_struc(n)%w212,era5_struc(n)%w222,&
          era5_struc(n)%n112,era5_struc(n)%n122,&
          era5_struc(n)%n212,era5_struc(n)%n222,&
          LIS_rc%udef, iret)
     
  elseif(trim(LIS_rc%met_interp(findex)).eq."bilinear".or.&
       trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          era5_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5_struc(n)%w111,era5_struc(n)%w121,&
          era5_struc(n)%w211,era5_struc(n)%w221,&
          era5_struc(n)%n111,era5_struc(n)%n121,&
          era5_struc(n)%n211,era5_struc(n)%n221,&
          LIS_rc%udef, iret)
     
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:),era5_struc(n)%mi,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          era5_struc(n)%n113,LIS_rc%udef,iret)
  else
     write(LIS_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(LIS_rc%met_interp(findex))//&
          ' not supported for ERA5'
     call LIS_endrun()
  endif
  
  if( pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0 ) then 
     
     call LIS_pcpClimoDownscaling(n, findex, month,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n), output_var(:), lo)
     
  endif

end subroutine interp_era5_var

!BOP
! 
! !ROUTINE: assign_processed_era5forc
! \label{assign_processed_era5forc}
! 
! !INTERFACE: 
subroutine assign_processed_era5forc(n,kk,order,var_index,era5forc)
! !USES: 
  use LIS_coreMod
  use era5_forcingMod, only : era5_struc
! 
! !DESCRIPTION: 
! 
!  This subroutine assigns the processed ERA5 forcing variables to 
!  the module datastructures (for use in temporal interpolation)
!
!EOP
  implicit none

  integer :: n
  integer :: kk
  integer :: order
  integer :: var_index
  real    :: era5forc(LIS_rc%lnc(n)*LIS_rc%lnc(n))
  

  integer :: c,r

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              era5_struc(n)%metdata1(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LIS_rc%lnc(n))
           elseif(order.eq.2) then 
              era5_struc(n)%metdata2(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   era5forc(c+(r-1)*LIS_rc%lnc(n))
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_era5forc
