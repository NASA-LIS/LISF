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
! !ROUTINE: read_gefs_operational
!  \label{read_gefs_operational}
!
! !REVISION HISTORY:
! 28 Jan 2021: Sujay Kumar; Updated for GEFS operational data
!
! !INTERFACE:
subroutine read_gefs_operational(n, m, findex, order, filename, ferror)

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_metforcingMod,  only : LIS_forc
  use gefs_forcingMod,    only : gefs_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)       :: findex
  integer, intent(in)       :: n
  integer, intent(in)       :: m
  integer, intent(in)       :: order
  character*140, intent(in) :: filename
  integer, intent(out)      :: ferror
!
! !DESCRIPTION:
!  For the given time, reads the forcing data from the 
!  GEFS file, transforms into LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[filename]
!    name of the file to be read
!  \item[ferror]
!    return error flag (0-fail, 1-success)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gefs](\ref{interp_gefs}) \newline
!    Performs spatial interpolation of GEFS forecast data to the LIS grid
!  \end{description}

!EOP
  integer     :: ftn, rc
  logical     :: file_exists  
  integer     :: c,r,i,k
  integer     :: nvars

  ! Grib specific dimenstions:
  integer     :: igrib
  integer     :: nfcsttimes
  integer     :: numpts
  integer     :: grib_gridsize

  ! Eccodes - Keys:
  integer      :: level
  character*50     :: shortName
  ! For more Grib Eccodes keys, see web document:
  !  https://confluence.ecmwf.int/download/attachments/97363968/eccodes-keys-2018.pdf?api=v2
  !  Slides: 7 to 20, forecast keys on slide 14

  real, allocatable  :: gefs_grib_data(:)              ! Read-in data
  real        :: eval, es
  real        :: relhum(LIS_rc%lnc(n),LIS_rc%lnr(n)) 
  real        :: sp(LIS_rc%lnc(n),LIS_rc%lnr(n)) 
  real        :: tair(LIS_rc%lnc(n),LIS_rc%lnr(n)) 
  real        :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n)) ! Interp field
  logical     :: pcp_flag                              ! Precip flag for spatial interp

! ______________________________________________________________________________

  varfield = 0 
  ferror = 1
  numpts = gefs_struc(n)%nc*gefs_struc(n)%nr
  allocate( gefs_grib_data(numpts) )
  gefs_grib_data = LIS_rc%udef

  pcp_flag = .false.

#if(defined USE_GRIBAPI)

! Check if file exists: ----------------------------------

  inquire( file=filename, exist=file_exists )
  if (file_exists) then      

     call grib_open_file(ftn,trim(filename),'r',rc)
     call LIS_verify(rc,'grib_open_file error in read_gefs_operational')

     call grib_count_in_file(ftn,nvars,rc)
     call LIS_verify(rc,'grib_count_in_file error in read_gefs_operational')
     
     do k=1,nvars
        call grib_new_from_file(ftn,igrib,rc)
        call LIS_verify(rc,'error in grib_new_from_file for data1 in read_gefs_operational')
         
        call grib_get(igrib,'shortName',shortName,rc)
        call LIS_verify(rc,'error in grib_get: shortName in read_gefs_operational')

        call grib_get(igrib,'level',level,rc)
        call LIS_verify(rc,'error in grib_get: level in read_gefs_operational')
                

        if(shortName.eq."2t") then !air_temperature
           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
        
           tair=varfield

           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gefs_struc(n)%metdata1(1,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    elseif(order.eq.2) then 
                       gefs_struc(n)%metdata2(1,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    endif
                 endif
              end do
           enddo
        elseif(shortName.eq."2r".and.level.eq.2) then !relative_humidity

           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
        
           relhum = varfield
          

        elseif(shortName.eq."dswrf") then !shortwave radiation
           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
                   
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gefs_struc(n)%metdata1(3,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    elseif(order.eq.2) then 
                       gefs_struc(n)%metdata2(3,m,&
                            LIS_domain(n)%gindex(c,r))&
                            = varfield(c,r)
                    endif
                 endif
              end do
           enddo
        elseif(shortName.eq."dlwrf") then !shortwave radiation
           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
        
          
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gefs_struc(n)%metdata1(4,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    elseif(order.eq.2) then 
                       gefs_struc(n)%metdata2(4,m,&
                            LIS_domain(n)%gindex(c,r))&
                            = varfield(c,r)
                    endif
                 endif
              end do
           enddo
        elseif(shortname.eq."10u".and.level.eq.10) then !Uwind

           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
                  

           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gefs_struc(n)%metdata1(5,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    elseif(order.eq.2) then 
                       gefs_struc(n)%metdata2(5,m,&
                            LIS_domain(n)%gindex(c,r))&
                            = varfield(c,r)
                    endif
                 endif
              end do
           enddo

        elseif(shortName.eq."10v".and.level.eq.10) then !vwind
           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
                  
          
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then
                    if(order.eq.1) then 
                       gefs_struc(n)%metdata1(6,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    elseif(order.eq.2) then 
                       gefs_struc(n)%metdata2(6,m,&
                            LIS_domain(n)%gindex(c,r))&
                            = varfield(c,r)
                    endif
                 endif
              end do
           enddo
        elseif(shortName.eq."sp") then !surface pressure
           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
        
           sp = varfield

          
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gefs_struc(n)%metdata1(7,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    elseif(order.eq.2) then 
                       gefs_struc(n)%metdata2(7,m,&
                            LIS_domain(n)%gindex(c,r))&
                            = varfield(c,r)
                    endif
                 endif
              end do
           enddo

        elseif(shortName.eq."tp") then !total precip
           gefs_grib_data = LIS_rc%udef
        
           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_operational')
        
           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                numpts, gefs_grib_data )
           
           pcp_flag = .true.
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
                  
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    if(order.eq.1) then 
                       gefs_struc(n)%metdata1(8,m,&
                            LIS_domain(n)%gindex(c,r)) &
                            = varfield(c,r)
                    elseif(order.eq.2) then 
                       gefs_struc(n)%metdata2(8,m,&
                            LIS_domain(n)%gindex(c,r))&
                            = varfield(c,r)
                    endif
                 endif
              end do
           enddo

        endif
        call grib_release(igrib,rc)
        call LIS_verify(rc,'error in grib_release in read_gefs_operational') 
     enddo

!  convert relative humidity to specific humidity using Bolton (1980)
!
!   e = 6.112*exp((17.67*Td)/(Td + 243.5)); 
!   q = (0.622 * e)/(p - (0.378 * e)); 
!    where: 
!       e = vapor pressure in mb; 
!       Td = dew point in deg C; 
!       p = surface pressure in mb; 
!       q = specific humidity in kg/kg. 
!
!   Td dew point temperature is
!   es = 6.112 * exp((17.67 * T)/(T + 243.5)); 
!   e = es * (RH/100.0); 
!   Td = log(e/6.112)*243.5/(17.67-log(e/6.112)); 
!     where: 
!       T = temperature in deg C; 
!       es = saturation vapor pressure in mb; 
!       e = vapor pressure in mb; 
!       RH = Relative Humidity in percent; 
!       Td = dew point in deg C 

 
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 tair(c,r) =tair(c,r)-273.15 !to celcius
                 es = 6.112*exp((17.67+tair(c,r))/(tair(c,r)+243.5))
                 eval = es*(relhum(c,r)/100.0)
                 gefs_struc(n)%metdata1(2,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = (0.622*eval)/(sp(c,r)*0.01-0.378*eval)

              elseif(order.eq.2) then 
                 tair(c,r) =tair(c,r)-273.15 !to celcius
                 es = 6.112*exp((17.67+tair(c,r))/(tair(c,r)+243.5))
                 eval = es*(relhum(c,r)/100.0)
                 
                 gefs_struc(n)%metdata2(2,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = (0.622*eval)/(sp(c,r)*0.01-0.378*eval)

              endif
           endif
        end do
     enddo
     

     call grib_close_file(ftn)
  
! GEFS file not found:
  else
     write(LIS_logunit,*) &
       '[ERR] Could not find file: ',trim(filename)
     ferror = 0
  endif

#endif
  deallocate( gefs_grib_data)
end subroutine read_gefs_operational
