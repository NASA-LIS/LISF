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
! !ROUTINE: read_era5cds
! \label{read_era5cds}
! 
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code 
! 15 apr 2025: Hiroko Beaudoing, adopted ERA5 routines for the public CDS
!                               data format
!
! !INTERFACE:
subroutine read_era5cds(n, kk,order, year, month, day, hour, findex,          &
     instfile, avgfile, ferror)
! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod
  use LDT_FORC_AttributesMod
  use LDT_metforcingMod, only : LDT_forc
  use era5cds_forcingMod, only : era5cds_struc
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
  character(len=*), intent(in) :: instfile, avgfile
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  ERA5 data, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain. \newline
!
!  ERA5 model output variables used to force LIS are provided in 2 files:
!  1) inst, Instantaneous values, available every hour
!  2) accum, Time integrated values (accumulations) over the past hour.
!
!  ERA5 FORCING VARIABLES:
!  1. T         inst,  Atmospheric temperature (K) (2m)
!  2. D         inst,  Dewpoint temperature (K) (2m)
!  3. SSRD      accum, surface solar radiation downwards [J m**-2]
!  4. STRD      accum, surface thermal radiation downwards [J m**-2]
!  5. U         inst, zonal wind [m/s] (10m)
!  6. V         inst, meridional wind [m/s] (10m)
!  7. SP        inst, surface pressure [Pa]
!  8. LSP       accum, large scale precipitation [m]
!  9. CP        accum, convective precipitation [m]
!
!  also inst file includes:
!  Z: Geopotential (m**2 s**-2)
!
!  NOTE 1: be aware that ECMWF outputs large-scale and convective precipitation
!  separately.  For total precipitation, need to sum the two fields,
!  LSP+CP=TP.
!  NOTE 2: only time2 SW & LW flux accumulations used in interpolation
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
  logical   :: file_exists
  integer   :: c,r,t,iv
  integer   :: fret,ferror1,ferror2
  integer   :: timestep
  real      :: glbdata_i(9,LDT_rc%ngrid(n))
  real      :: glbdata_a(9,LDT_rc%ngrid(n))
! __________________________________________________________________________

  glbdata_i = LDT_rc%udef
  glbdata_a = LDT_rc%udef

  ferror1 = 1
  ferror2 = 1    ! 1 success
  fret = 0       ! 0 success

  timestep = hour + 1
  if ( timestep .ge. 24 ) timestep = 1

  call retrieve_inst_era5cdsvars(n, findex, timestep, instfile, glbdata_i, fret)
  if (fret.ne.0) then
     ferror1 = 0
  endif

  if ( order.eq.2 ) then
  call retrieve_accum_era5cdsvars(n, findex, avgfile, glbdata_a, &
                                timestep, fret)
   if (fret.ne.0) then
      ferror2 = 0
   endif

  endif  ! order==2

  if ( ferror1 == 1 .and. ferror2 == 1) then
   ferror = 1     ! success
  else
   ferror = 0     ! problem, roll back one day
  endif
  if ( ferror == 1 ) then  ! only proceed if retrieve calls were successful
   do iv=1,9
     do t=1,LDT_rc%ngrid(n)
        ! these are time avgd fields
        if ( iv.eq.3 .or. iv.eq.4 .or. iv.eq.8 .or. iv.eq.9 ) then
           if(order.eq.1) then
              era5cds_struc(n)%metdata1(iv,t) = glbdata_a(iv,t)
           else
              era5cds_struc(n)%metdata2(iv,t) = glbdata_a(iv,t)
           endif
        ! these are instantaneous
        else
           if(order.eq.1) then
              era5cds_struc(n)%metdata1(iv,t) = glbdata_i(iv,t)
           else
              era5cds_struc(n)%metdata2(iv,t) = glbdata_i(iv,t)
           endif
        endif
     enddo
   enddo
  endif    ! ferror == 1

end subroutine read_era5cds

!BOP
! !ROUTINE: retrieve_inst_era5cdsvars
! \label{retrieve_inst_era5cdsvars}
!
! !INTERFACE:
subroutine retrieve_inst_era5cdsvars(n, findex, timestep, month, instfile, glbdata, fret)
  use LDT_coreMod,      only : LDT_rc, LDT_domain
  use LDT_logMod,       only : LDT_logunit, LDT_verify
  use era5cds_forcingMod, only : era5cds_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none
! !ARGUMENTS:
  integer,   intent(in)        :: n
  integer,   intent(in)        :: findex
  integer,   intent(in)        :: timestep
  integer,   intent(in)        :: month
  character(len=*), intent(in) :: instfile
  real                         :: glbdata(9,LDT_rc%ngrid(n))
  integer, intent(inout)       :: fret
!
! !DESCRIPTION:
! This routine opens the corresponding ERA5CDS data file to extract
! the specified variable, which represents an instantaneous value.
! Should be used for near-surface temperature, specific humidity,
! winds, and surface pressure.
!
!EOP
  integer,parameter  :: N_INST_VARS=5
  integer            :: iret
  integer            :: c,r,iv,v,i,j
  integer            :: varid
  real               :: missingValue
  real               :: e,es,td,ta
  logical            :: file_exists
  character(len=3)   :: varname(N_INST_VARS)
  real,  allocatable     :: f(:)
  logical*1, allocatable :: lb(:)
  logical            :: pcp_flag
  integer            :: rel_index(N_INST_VARS)
  real, allocatable, dimension(:,:)  :: datain
  real, dimension(LDT_rc%lnc(n)*LDT_rc%lnr(n))  :: varfield
  !=== set variable name
  varname = (/ "t2m","d2m","u10","v10","sp " /)
  rel_index = (/ 1, 2, 5, 6, 7 /) ! index of variable w.r.t. the
                                  ! full list of 9 forcing variables
  inquire(file=instfile,exist=file_exists)
  if ( file_exists ) then
    allocate(datain(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold))
    allocate(f(era5cds_struc(n)%ncold*era5cds_struc(n)%nrold))
    allocate(lb(era5cds_struc(n)%ncold*era5cds_struc(n)%nrold))

    call LDT_verify(nf90_open(path=trim(instfile), mode=NF90_NOWRITE, &
             ncid=iret), 'nf90_open failed in read_era5cds')
    ! Forcing variable loop:
    do v = 1, N_INST_VARS
       iv = rel_index(v)
       call LDT_verify(nf90_inq_varid(iret, trim(varname(v)),varId), &
             'nf90_inq_varid failed in retrieve_inst')
       call LDT_verify(nf90_get_var(iret,varId, datain, &
             start=(/1,1,timestep/), &
             count=(/era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,1/)),&
             'nf90_get_var failed in retrieve_inst')
       call LDT_verify(nf90_get_att(iret,varId, "_FillValue", missingValue), &
             'nf90_get_att failed in retrieve_inst for FillValue')

     !=== Transferring current data to 1-D array for interpolation
     c=0
     do j=1,era5cds_struc(n)%nrold
       do i=1,era5cds_struc(n)%ncold
          c = c + 1
          f(c) = datain(i,j)
          if ( f(c) .ne. missingValue ) then
             lb(c) = .true.
          endif
       enddo
     enddo
     pcp_flag = .false.
     call interp_era5cds_var(n,findex,month,f, &
                             lb,pcp_flag,varfield)

     if ( iv .eq. 2 ) then   ! convert from dew point to qair
       do c=1,LDT_rc%ngrid(n)
         td = varfield(c) - 273.15   ! in celsius
         ta = glbdata(1,c) - 273.15   ! in celsius
         e = 6.11* 10**(7.5*td / (237.3+td)) ! actual vapor pressure
         es = 6.11*10**(7.5*ta / (237.3+ta)) ! saturation vapor pressure
         glbdata(iv,c) = e / es
       end do
     else
       glbdata(iv,:) = varfield
     end if

    end do   ! v

    call LDT_verify(nf90_close(iret), 'failed to close file in retrieve_inst')
    deallocate(f)
    deallocate(lb)
    deallocate(datain)
    fret = 0
 else
    fret = -1
 endif

end subroutine retrieve_inst_era5cdsvars
!BOP
! !ROUTINE: retrieve_accum_ecmwfvars
! \label{retrieve_accum_ecmwfvars}
!
! !INTERFACE:
subroutine retrieve_accum_era5cdsvars(n, findex, month, avgfile, glbdata1, &
                                    timestep, fret)
  use LDT_coreMod,      only : LDT_rc, LDT_domain
  use LDT_logMod,       only : LDT_logunit, LDT_verify
  use era5cds_forcingMod, only : era5cds_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  implicit none
! !ARGUMENTS:
  integer,   intent(in)        :: n
  integer,   intent(in)        :: findex
  integer,   intent(in)        :: month
  character(len=*), intent(in) :: avgfile
  real, intent(inout)          :: glbdata1(9,LDT_rc%ngrid(n))
  integer                      :: timestep
  integer, intent(inout)       :: fret
!
! !DESCRIPTION:
! This routine opens the corresponding ERA5CDS data file to extract
! the specified variable, which represents an accumulated value.
! Should be used for shortwave, longwave, lsp, and cp.
!
!EOP
  integer,parameter  :: N_ACCUM_VARS=4
  integer            :: nera5cds
  integer            :: iret,gbret
  integer            :: c,r,iv,ftn,i,j,v
  integer            :: varid
  real               :: missingValue
  logical            :: file_exists
  character(len=4)   :: varname(N_ACCUM_VARS)
  real,  allocatable     :: f(:)
  logical*1, allocatable :: lb(:)
  real, allocatable, dimension(:,:)  :: datain
  logical            :: pcp_flag
  integer            :: rel_index(N_ACCUM_VARS)
  real,dimension(LDT_rc%lnc(n)*LDT_rc%lnr(n)) :: varfield

  !=== set GRIB shortName
  varname = (/ "ssrd","strd","lsp ","cp  " /)
  rel_index = (/ 3, 4, 8, 9 /) ! index of variable w.r.t. the
                               ! full list of 9 forcing variables

  nera5cds = era5cds_struc(n)%ncold*era5cds_struc(n)%nrold

  inquire(file=avgfile,exist=file_exists)
  if ( file_exists ) then
    allocate(datain(era5cds_struc(n)%ncold,era5cds_struc(n)%nrold))
    allocate(f(nera5cds))
    allocate(lb(nera5cds))

    call LDT_verify(nf90_open(path=trim(avgfile), mode=NF90_NOWRITE, &
             ncid=iret), 'nf90_open failed in retrieve_avg_era5')
    ! Forcing variable loop:
    do v = 1, N_ACCUM_VARS
       iv = rel_index(v)
       call LDT_verify(nf90_inq_varid(iret, trim(varname(v)),varId), &
             'nf90_inq_varid failed in retrieve_avg')
       call LDT_verify(nf90_get_var(iret,varId, datain, &
             start=(/1,1,timestep/), &
             count=(/era5cds_struc(n)%ncold,era5cds_struc(n)%nrold,1/)),&
             'nf90_get_var failed in retrieve_avg')
       call LDT_verify(nf90_get_att(iret,varId, "_FillValue", missingValue), &
             'nf90_get_att failed in retrieve_avg for FillValue')

       !=== Transferring current data to 1-D array for interpolation
       c=0
       do j=1,era5cds_struc(n)%nrold
         do i=1,era5cds_struc(n)%ncold
            c = c + 1
            f(c) = datain(i,j)
            if ( f(c) .ne. missingValue ) then
               lb(c) = .true.
            endif
         enddo
       enddo
       pcp_flag = .false.
       if(iv.eq.8.or.iv.eq.9) pcp_flag = .true.

       call interp_era5cds_var(n,findex,month,f, &
                             lb,pcp_flag,varfield)

       if ( iv == 3 .or. iv == 4) then ! swd/lwd - convert accumulated field to rate
          do c=1,LDT_rc%ngrid(n)
             glbdata1(iv,c) = varfield(c) / (1.0*60*60)
          enddo
       elseif ( iv == 8 .or. iv == 9 ) then ! lsp or cp [m] -> mm/s
          do c=1,LDT_rc%ngrid(n)
             glbdata1(iv,c) = (varfield(c) * 1000.0)/(1.0*60*60)
          enddo
       end if
    end do   ! v

    do c=1,LDT_rc%ngrid(n)
       glbdata1(8,c) = glbdata1(8,c) + glbdata1(9,c)
    enddo

    call LDT_verify(nf90_close(iret), 'failed to close file in retrieve_avg_era5')
    deallocate(f)
    deallocate(lb)
    deallocate(datain)
    fret = 0
  else
    fret = -1
  end if

end subroutine retrieve_accum_era5cdsvars
!BOP
! 
! !ROUTINE: interp_era5cds_var
! \label{interp_era5cds_var}
! 
! !INTERFACE: 
subroutine interp_era5cds_var(n,findex, month, input_var, &
     pcp_flag, output_var)

! !USES: 
  use LDT_coreMod
  use LDT_logMod
  use era5cds_forcingMod, only : era5cds_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(era5cds_struc(n)%npts)
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
  real      :: f (era5cds_struc(n)%ncold*era5cds_struc(n)%nrold)
  logical*1 :: lb(era5cds_struc(n)%ncold*era5cds_struc(n)%nrold)
  logical*1 :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer   :: input_size
! _____________________________________________________________

  input_size = era5cds_struc(n)%ncold*era5cds_struc(n)%nrold
  output_var = LDT_rc%udef
  lo = .true.

!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    

  if(pcp_flag.and.&
       trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
     
     call conserv_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          era5cds_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),& 
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5cds_struc(n)%w112,era5cds_struc(n)%w122,&
          era5cds_struc(n)%w212,era5cds_struc(n)%w222,&
          era5cds_struc(n)%n112,era5cds_struc(n)%n122,&
          era5cds_struc(n)%n212,era5cds_struc(n)%n222,&
          LDT_rc%udef, iret)
     
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear".or.&
       trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 

     call bilinear_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          era5cds_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n), & 
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5cds_struc(n)%w111,era5cds_struc(n)%w121,&
          era5cds_struc(n)%w211,era5cds_struc(n)%w221,&
          era5cds_struc(n)%n111,era5cds_struc(n)%n121,&
          era5cds_struc(n)%n211,era5cds_struc(n)%n221,&
          LDT_rc%udef, iret)
     
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then 
     call neighbor_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:),era5cds_struc(n)%mi,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n),&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          era5cds_struc(n)%n113,LDT_rc%udef,iret)
  else
     write(LDT_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(LDT_rc%met_gridtransform(findex))//&
          ' not supported for ERA5CDS'
     call LDT_endrun()
  endif
  

end subroutine interp_era5cds_var
