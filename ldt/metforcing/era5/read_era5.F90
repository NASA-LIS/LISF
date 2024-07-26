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
  integer   :: c,r,t,k,l,iret
  integer   :: mo,rec_size
  logical   :: read_lnd
  logical   :: read_flag
  
  real, allocatable      :: tair(:,:)
  real, allocatable      :: qair(:,:)
  real, allocatable      :: swd(:,:)
  real, allocatable      :: lwd(:,:)
  real, allocatable      :: wind(:,:)
  real, allocatable      :: ps(:,:)
  real, allocatable      :: rainf(:,:)
  real, allocatable      :: snowf(:,:)
  real, allocatable      :: dirsw(:,:)
  real, allocatable      :: difsw(:,:)

  integer :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/
! __________________________________________________________________________

  read_flag = .false.
  ! Check if it is the switch of a month
  if(order.eq.1) then
     if(era5_struc(n)%mo1.ne.month) then
        era5_struc(n)%mo1 = month
        read_flag = .true.
     endif
  else
     if(era5_struc(n)%mo2.ne.month) then
        era5_struc(n)%mo2 = month
        read_flag = .true.
     endif
  endif
  ferror = 1

  if(read_flag) then 
#if (defined USE_NETCDF4) 

     if((mod(year,4) .eq. 0 .and. mod(year, 100).ne.0) &!leap year
          .or.(mod(year,400) .eq.0)) then 
        days(2) = 29
     else 
        days(2) = 28
     endif
     
     mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

! Read single layer file (*slv) fields:
     inquire(file=fname,exist=file_exists) 
     if(file_exists) then 
        if(order.eq.1) then
           era5_struc(n)%ps1     = LDT_rc%udef
           era5_struc(n)%tair1   = LDT_rc%udef
           era5_struc(n)%qair1   = LDT_rc%udef
           era5_struc(n)%wind1   = LDT_rc%udef
           era5_struc(n)%rainf1  = LDT_rc%udef
           era5_struc(n)%swd1    = LDT_rc%udef
           era5_struc(n)%lwd1    = LDT_rc%udef
        else
           era5_struc(n)%ps2     = LDT_rc%udef
           era5_struc(n)%tair2   = LDT_rc%udef
           era5_struc(n)%qair2   = LDT_rc%udef
           era5_struc(n)%wind2   = LDT_rc%udef
           era5_struc(n)%rainf2  = LDT_rc%udef
           era5_struc(n)%swd2    = LDT_rc%udef
           era5_struc(n)%lwd2    = LDT_rc%udef
        endif

        write(LDT_logunit,*)'[INFO] Reading ERA5 file (bookend,', order,' -',trim(fname), ')'

        rec_size = days(month)*24 + 1

        call LDT_verify(nf90_open(path=trim(fname), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_era5')

        allocate(ps(era5_struc(n)%npts,rec_size))
        
        call LDT_verify(nf90_inq_varid(ftn,'PSurf',psId), &
             'nf90_inq_varid failed for psurf in read_era5')
        call LDT_verify(nf90_get_var(ftn,psId, ps),&
             'nf90_get_var failed for ps in read_era5') 

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_era5_var(n,findex,month,ps(:,l),   .false.,&
                   era5_struc(n)%ps1(:,l))              
           enddo
        else
           do l=1,rec_size
              call interp_era5_var(n,findex,month,ps(:,l),    .false.,&
                   era5_struc(n)%ps2(:,l))
           enddo
        endif
        deallocate(ps)

        allocate(tair(era5_struc(n)%npts,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
             'nf90_inq_varid failed for Tair in read_era5')
        call LDT_verify(nf90_get_var(ftn,tmpId, tair),&
             'nf90_get_var failed for tair in read_era5') 

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_era5_var(n,findex,month,tair(:,l), .false.,&
                   era5_struc(n)%tair1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_era5_var(n,findex,month,tair(:,l),  .false.,&
                   era5_struc(n)%tair2(:,l))
           enddo
        endif
        deallocate(tair)

        allocate(qair(era5_struc(n)%npts,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Qair',qId), &
             'nf90_inq_varid failed for Qair in read_era5')
        call LDT_verify(nf90_get_var(ftn,qId, qair),&
             'nf90_get_var failed for qair in read_era5') 

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_era5_var(n,findex,month,qair(:,l), .false.,&
                   era5_struc(n)%qair1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_era5_var(n,findex,month,qair(:,l),  .false.,&
                   era5_struc(n)%qair2(:,l))
           enddo
        endif
        deallocate(qair)

        allocate(wind(era5_struc(n)%npts,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Wind',windId), &
             'nf90_inq_varid failed for Wind in read_era5')
        call LDT_verify(nf90_get_var(ftn,windId, wind),&
             'nf90_get_var failed for wind in read_era5') 
        if(order.eq.1) then 
           do l=1,rec_size
              call interp_era5_var(n,findex,month,wind(:,l), .false.,&
                   era5_struc(n)%wind1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_era5_var(n,findex,month,wind(:,l),  .false.,&
                   era5_struc(n)%wind2(:,l))
           enddo
        endif
        deallocate(wind)
        
        allocate(rainf(era5_struc(n)%npts,rec_size))
        allocate(snowf(era5_struc(n)%npts,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
             'nf90_inq_varid failed for Rainf in read_era5')
        call LDT_verify(nf90_get_var(ftn,rainfId, rainf),&
             'nf90_get_var failed for rainf in read_era5') 
        call LDT_verify(nf90_inq_varid(ftn,'Snowf',snowfId), &
             'nf90_inq_varid failed for Snowf in read_era5')
        call LDT_verify(nf90_get_var(ftn,snowfId, snowf),&
             'nf90_get_var failed for snowf in read_era5') 

        do t=1,era5_struc(n)%npts
           do l=1,rec_size
              if(rainf(t,l).ne.LDT_rc%udef.and.&
                   snowf(t,l).ne.LDT_rc%udef) then 
                 rainf(t,l) = rainf(t,l) + snowf(t,l)
              else
                 rainf(t,l) = LDT_rc%udef
              endif
           enddo
        enddo
        deallocate(snowf)
        if(order.eq.1) then 
           do l=1,rec_size
              call interp_era5_var(n,findex,month,rainf(:,l),.true.,&
                   era5_struc(n)%rainf1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_era5_var(n,findex,month,rainf(:,l), .true.,&
                   era5_struc(n)%rainf2(:,l))
           enddo
        endif
        deallocate(rainf)

        allocate(difSW(era5_struc(n)%npts,rec_size))
        allocate(swd(era5_struc(n)%npts,rec_size))

        call LDT_verify(nf90_inq_varid(ftn,'DIR_SWdown',dirSWId), &
             'nf90_inq_varid failed for DIR_SWdown in read_era5')
        call LDT_verify(nf90_get_var(ftn,dirSWId, swd),&
             'nf90_get_var failed for dirSW in read_era5') 
        iret = nf90_inq_varid(ftn,'SCA_SWdown',difSWId)
        !Assumes that diffuse field doesn't exist
        if(iret.ne.0) then 
           difsw = 0 
        endif

        do t=1,era5_struc(n)%npts
           do l=1,rec_size
              if(swd(t,l).ne.LDT_rc%udef.and.&
                   difsw(t,l).ne.LDT_rc%udef) then 
                 swd(t,l) = swd(t,l) + difsw(t,l)
              else
                 swd(t,l) = LDT_rc%udef
              endif
           enddo
        enddo
        deallocate(difsw)

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_era5_var(n,findex,month,swd(:,l),  .false.,&
                   era5_struc(n)%swd1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_era5_var(n,findex,month,swd(:,l),   .false.,&
                   era5_struc(n)%swd2(:,l))
           enddo
        endif
        deallocate(swd)

        allocate(lwd(era5_struc(n)%npts,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'LWdown',lwdId), &
             'nf90_inq_varid failed for LWdown in read_era5')
        call LDT_verify(nf90_get_var(ftn,lwdId, lwd),&
             'nf90_get_var failed for lwd in read_era5')

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_era5_var(n,findex,month,lwd(:,l),  .false.,&
                   era5_struc(n)%lwd1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_era5_var(n,findex,month,lwd(:,l),   .false.,&
                   era5_struc(n)%lwd2(:,l))
           enddo
        endif
        deallocate(lwd)
        
        call LDT_verify(nf90_close(ftn), &
             'failed to close file in read_era5')
                
     else
        write(LDT_logunit,*) '[ERR] ',trim(fname)//' does not exist'
        call LDT_endrun()
        
     endif
  endif

  tindex = (day - 1)*24 + hour + 1
  
  if(order.eq.1) then 
     call assign_processed_era5forc(n,kk,order,1,era5_struc(n)%tair1(:,tindex))
     call assign_processed_era5forc(n,kk,order,2,era5_struc(n)%qair1(:,tindex))
     call assign_processed_era5forc(n,kk,order,3,era5_struc(n)%swd1(:,tindex))
     call assign_processed_era5forc(n,kk,order,4,era5_struc(n)%lwd1(:,tindex))
     call assign_processed_era5forc(n,kk,order,5,era5_struc(n)%wind1(:,tindex))
     call assign_processed_era5forc(n,kk,order,6,era5_struc(n)%ps1(:,tindex))
     call assign_processed_era5forc(n,kk,order,7,era5_struc(n)%rainf1(:,tindex))
  else
     call assign_processed_era5forc(n,kk,order,1,era5_struc(n)%tair2(:,tindex))
     call assign_processed_era5forc(n,kk,order,2,era5_struc(n)%qair2(:,tindex))
     call assign_processed_era5forc(n,kk,order,3,era5_struc(n)%swd2(:,tindex))
     call assign_processed_era5forc(n,kk,order,4,era5_struc(n)%lwd2(:,tindex))
     call assign_processed_era5forc(n,kk,order,5,era5_struc(n)%wind2(:,tindex))
     call assign_processed_era5forc(n,kk,order,6,era5_struc(n)%ps2(:,tindex))
     call assign_processed_era5forc(n,kk,order,7,era5_struc(n)%rainf2(:,tindex))
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
  use LDT_coreMod
  use LDT_logMod
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
  real,    intent(out)   :: output_var(LDT_rc%lnc(n)*LDT_rc%lnr(n))
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
! _____________________________________________________________

  input_size = era5_struc(n)%ncold*era5_struc(n)%nrold
  output_var = LDT_rc%udef

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
  
!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    
     

  if(pcp_flag.and.&
       trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
     
     call conserv_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
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
          output_var(:), &
          era5_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n), & 
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5_struc(n)%w111,era5_struc(n)%w121,&
          era5_struc(n)%w211,era5_struc(n)%w221,&
          era5_struc(n)%n111,era5_struc(n)%n121,&
          era5_struc(n)%n211,era5_struc(n)%n221,&
          LDT_rc%udef, iret)
     
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then 
     call neighbor_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:),era5_struc(n)%mi,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n),&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5_struc(n)%n113,LDT_rc%udef,iret)
  else
     write(LDT_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(LDT_rc%met_gridtransform(findex))//&
          ' not supported for ERA5'
     call LDT_endrun()
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
  use LDT_coreMod
  use era5_forcingMod, only : era5_struc
!
! !DESCRIPTION: 
!  This routine assigns the interpolated ERA5 forcing data
!  to the module data structures to be used later for 
!  time interpolation 
!
!EOP
  implicit none

  integer :: n
  integer :: kk
  integer :: order
  integer :: var_index
  real    :: era5forc(LDT_rc%lnc(n)*LDT_rc%lnc(n))
  

  integer :: c,r

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
end subroutine assign_processed_era5forc






