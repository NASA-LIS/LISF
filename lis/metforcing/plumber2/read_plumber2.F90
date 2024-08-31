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
! !ROUTINE: read_plumber2
! \label{read_plumber2}
! 
! !REVISION HISTORY:
! 16 Sep 2021: Mark Beauharnois, initial version based on PUMET reader
!
! !INTERFACE:      
subroutine read_plumber2(n,order,findex,&
     filename, metforc,ferror)
! !USES:
  use LIS_coreMod, only : LIS_rc,LIS_domain, LIS_masterproc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use plumber2_forcingMod, only : plumber2_struc

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: order
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: filename
  real, intent(inout) :: metforc(LIS_rc%met_nf(findex),LIS_rc%ngrid(n))
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  PLUMBER2 data file, transforms into LIS_met_nf(findex) LIS forcing 
!  parameters. Spatial interpolation for PLUMBER2 removed. 

!  Index order:
    !  1. prec rain (total precip)
    !  2. psurf
    !  3. qair
    !  4. tair
    !  5. swdown
    !  6. lwdown
    !  7. wind u
    !  8. wind v
    !  9. lai    currently not used

!  lwdown surface_downwelling_longwave_flux_in_air [W/m^2]
!  precip total precipitation_flux [kg/m^2/s]
!  psurf surface air_pressure [Pa]
!  spfh specific_humidity [ ]
!  swdown surface_downwelling_shortwave_flux_in_air [W/m^2]
!  tair 2m air_temperature [K]
!  wind wind_speed [m/s] 
!  lai leaf_area_index[-] 
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    instance, order=2, read the next instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the PLUMBER2 Met Forcing file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!EOP
  
  integer   :: rIndex
  integer   :: ftn
  integer   :: tairId, qairId, windId, psurfId
  integer   :: precId, swdId, lwdId
  logical   :: file_exists
  integer   :: c,r,t,k
  real      :: tairData,qairData,windData,precData
  real      :: swdData,lwdData,psurfData
  real*8    :: FVal

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 

  ferror = 0 

  rIndex = plumber2_struc(n)%read_index

  inquire(file=filename,exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*)'[INFO] Reading PLUMBER2 file (bookend,', order,' ... ',trim(filename)
     call LIS_verify(nf90_open(path=trim(filename), mode=NF90_SHARE, &
          ncid=ftn), 'nf90_open failed for PLUMBER2 forcing file in read_plumber2')

!    -----------------------READ lwdown--------------------------------------------------

     call LIS_verify(nf90_inq_varid(ftn,'LWdown',lwdId), &
          'nf90_inq_varid failed for LWdown in read_plumber2')
     call LIS_verify(nf90_get_var(ftn,lwdId,lwdData, &
           start=(/1,1,rIndex/)), &
          'nf90_get_var failed for LWdown in read_plumber2')
     call LIS_verify(nf90_get_att(ftn, lwdId,'_FillValue', FVal), &
          'nf90_get_att failed for LWdown:_FillValue in read_plumber2') 

!    -------------------------READ Precip Flux-----------------------------------------------

     call LIS_verify(nf90_inq_varid(ftn,'Precip',precId), &
          'nf90_inq_varid failed for precipitation flux in read_plumber2')
     call LIS_verify(nf90_get_var(ftn,precId,precData, &
           start=(/1,1,rIndex/)), &
          'nf90_get_var failed for precipitation flux in read_plumber2')
     call LIS_verify(nf90_get_att(ftn, precId,'_FillValue', FVal), &
          'nf90_get_att failed for Precip:_FillValue in read_plumber2') 

!    -------------------------READ psurf------------------------------------------------

     call LIS_verify(nf90_inq_varid(ftn,'Psurf',psurfId), &
          'nf90_inq_varid failed for SURFACE PRESSURE (psurf) in read_plumber2')
     call LIS_verify(nf90_get_var(ftn,psurfId, psurfData, &
          start=(/1,1,rIndex/)), &
          'nf90_get_var failed for SURFACE PRESSURE (psurf) in read_plumber2') 
     call LIS_verify(nf90_get_att(ftn, psurfId,'_FillValue', FVal), &
          'nf90_get_att failed for psurf:_FillValue in read_plumber2') 

!    --------------------------READ spfh-----------------------------------------------

     call LIS_verify(nf90_inq_varid(ftn,'Qair',qairId), &
          'nf90_inq_varid failed for Sp. Humidity (spfh) in read_plumber2')
     call LIS_verify(nf90_get_var(ftn,qairId, qairData, &
          start=(/1,1,rIndex/)), &
          'nf90_get_var failed for Sp. Humidity (spfh) in read_plumber2')
     call LIS_verify(nf90_get_att(ftn, qairId,'_FillValue', FVal), &
          'nf90_get_att failed for spfh:_FillValue in read_plumber2') 

!    --------------------------READ swdown-----------------------------------------------

     call LIS_verify(nf90_inq_varid(ftn,'SWdown',swdId), &
          'nf90_inq_varid failed for SWdown in read_plumber2')
     call LIS_verify(nf90_get_var(ftn,swdId,swdData, &
          start=(/1,1,rIndex/)), &
          'nf90_get_var failed for SWdown in read_plumber2')
     call LIS_verify(nf90_get_att(ftn,swdId,'_FillValue', FVal), &
          'nf90_get_att failed for SWdown:_FillValue in read_plumber2') 

!    -------------------------READ tair------------------------------------------------

     call LIS_verify(nf90_inq_varid(ftn,'Tair',tairId), &
          'nf90_inq_varid failed for surface air temp (Tair) in read_plumber2')
     call LIS_verify(nf90_get_var(ftn,tairId,tairData, &
        start=(/1,1,rIndex/)), &
       'nf90_get_var failed for surface air temperature (Tair) in read_plumber2')
     call LIS_verify(nf90_get_att(ftn,tairId,'_FillValue', FVal), &
          'nf90_get_att failed for Tair:_FillValue in read_plumber2') 

!    --------------------------READ wind-----------------------------------------------

     call LIS_verify(nf90_inq_varid(ftn,'Wind',windId), &
          'nf90_inq_varid failed for scalar wind in read_plumber2')
     call LIS_verify(nf90_get_var(ftn,windId,windData, &
          start=(/1,1,rIndex/)), &
          'nf90_get_var failed for scalar wind in read_plumber2')
     call LIS_verify(nf90_get_att(ftn,windId,'_FillValue', FVal), &
          'nf90_get_att failed for Wind:_FillValue in read_plumber2') 

!    --------------------------READ LAI-----------------------------------------------

!    call LIS_verify(nf90_inq_varid(ftn,'LAI',laiId), &
!         'nf90_inq_varid failed for LAI in read_plumber2')
!    call LIS_verify(nf90_get_var(ftn,laiId,laiData), &
!         'nf90_get_var failed for LAI in read_plumber2')
!    call LIS_verify(nf90_get_att(ftn,laiId,'_FillValue', FVal), &
!         'nf90_get_att failed for LAI:_FillValue in read_plumber2') 

!    =========================================================================
!
    !  1. prec rain (total precip)
    !  2. psurf
    !  3. qair
    !  4. tair
    !  5. swdown
    !  6. lwdown
    !  7. wind u
    !  8. wind v
    !  9. lai   currently not used
  
    metforc(1, LIS_rc%ngrid(n)) = precData
    metforc(2, LIS_rc%ngrid(n)) = psurfData
    metforc(3, LIS_rc%ngrid(n)) = qairData
    metforc(4, LIS_rc%ngrid(n)) = tairData
    metforc(5, LIS_rc%ngrid(n)) = swdData
    metforc(6, LIS_rc%ngrid(n)) = lwdData
    metforc(7, LIS_rc%ngrid(n)) = windData   ! u wind
    metforc(8, LIS_rc%ngrid(n)) = 0.0        ! v wind
!   metforc(9, LIS_rc%ngrid(n)) = laiData

  else
     write(LIS_logunit,*) trim(filename)//' does not exist'
     call LIS_endrun()
  endif
  
  call LIS_verify(nf90_close(ftn),&
       'nf90_close failed for PLUMBER2 forcing file in read_plumber2')

#endif
end subroutine read_plumber2
